# R/10_forecast_fpca_VAR_1_7_holiday_regime.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)
library(lubridate)
library(readr)

dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

K <- 2
NBASIS <- 12
NORDER <- 4
TEST_START <- as.Date("2023-01-01")
MIN_TRAIN_DAYS <- 365
MIN_LAG_DAYS <- 28
DEBUG_N <- Inf

MODEL_NAME_BASE <- "FPCA VAR(1,7) + holiday/load regime"

# ------------------------------------------------------------
# Load curves
# ------------------------------------------------------------

curves <- readRDS("data/processed/load_curves.rds")

dates <- as.Date(curves$dates)
hours <- curves$hours

Y <- curves$actual
colnames(Y) <- as.character(hours)

has_official <- !is.null(curves$forecast)

if (has_official) {
  Y_official <- curves$forecast
  colnames(Y_official) <- as.character(hours)
}

ok <- complete.cases(Y)

dates <- dates[ok]
Y <- Y[ok, , drop = FALSE]

if (has_official) {
  Y_official <- Y_official[ok, , drop = FALSE]
}

n_days <- nrow(Y)

hcol <- function(h) {
  out <- match(h, hours)
  if (is.na(out)) stop("Hour ", h, " not found.")
  out
}

H5 <- hcol(5)
H6 <- hcol(6)
H7 <- hcol(7)
H9 <- hcol(9)

# ------------------------------------------------------------
# Calendar and holiday features
# ------------------------------------------------------------

date_seq <- function(start, end) seq(as.Date(start), as.Date(end), by = "day")

make_holiday_lookup <- function(years) {
  fixed_holidays <- bind_rows(lapply(years, function(y) {
    tibble(
      date = as.Date(c(
        paste0(y, "-01-01"), paste0(y, "-04-23"), paste0(y, "-05-01"),
        paste0(y, "-05-19"), paste0(y, "-07-15"), paste0(y, "-08-30"),
        paste0(y, "-10-29")
      )),
      holiday_type = "fixed"
    )
  }))
  
  religious_ranges <- tribble(
    ~start, ~end, ~holiday_type,
    "2019-06-04", "2019-06-06", "ramadan_bayram",
    "2020-05-24", "2020-05-26", "ramadan_bayram",
    "2021-05-13", "2021-05-15", "ramadan_bayram",
    "2022-05-02", "2022-05-04", "ramadan_bayram",
    "2023-04-21", "2023-04-23", "ramadan_bayram",
    "2024-04-10", "2024-04-12", "ramadan_bayram",
    "2025-03-30", "2025-04-01", "ramadan_bayram",
    "2019-08-11", "2019-08-14", "kurban_bayram",
    "2020-07-31", "2020-08-03", "kurban_bayram",
    "2021-07-20", "2021-07-23", "kurban_bayram",
    "2022-07-09", "2022-07-12", "kurban_bayram",
    "2023-06-28", "2023-07-01", "kurban_bayram",
    "2024-06-16", "2024-06-19", "kurban_bayram",
    "2025-06-06", "2025-06-09", "kurban_bayram"
  ) |>
    mutate(
      start = as.Date(start),
      end = as.Date(end)
    )
  
  religious_holidays <- bind_rows(lapply(seq_len(nrow(religious_ranges)), function(i) {
    tibble(
      date = date_seq(religious_ranges$start[i], religious_ranges$end[i]),
      holiday_type = religious_ranges$holiday_type[i]
    )
  }))
  
  bind_rows(fixed_holidays, religious_holidays) |>
    filter(year(date) %in% years) |>
    distinct(date, .keep_all = TRUE)
}

make_ramadan_lookup <- function(years) {
  ramadan_ranges <- tribble(
    ~start, ~end,
    "2019-05-06", "2019-06-03",
    "2020-04-24", "2020-05-23",
    "2021-04-13", "2021-05-12",
    "2022-04-02", "2022-05-01",
    "2023-03-23", "2023-04-20",
    "2024-03-11", "2024-04-09",
    "2025-03-01", "2025-03-29"
  ) |>
    mutate(
      start = as.Date(start),
      end = as.Date(end)
    )
  
  bind_rows(lapply(seq_len(nrow(ramadan_ranges)), function(i) {
    tibble(date = date_seq(ramadan_ranges$start[i], ramadan_ranges$end[i]))
  })) |>
    filter(year(date) %in% years) |>
    distinct(date)
}

holiday_lookup <- make_holiday_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))
ramadan_lookup <- make_ramadan_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))

make_calendar_features <- function(date_vec) {
  doy <- yday(date_vec)
  dow <- wday(date_vec, week_start = 1)
  month_num <- month(date_vec)
  holiday_dates <- holiday_lookup$date
  
  tibble(date = date_vec) |>
    mutate(
      dow = dow,
      weekend = as.integer(dow >= 6),
      monday = as.integer(dow == 1),
      friday = as.integer(dow == 5),
      month = month_num,
      sin_week = sin(2 * pi * dow / 7),
      cos_week = cos(2 * pi * dow / 7),
      sin_year = sin(2 * pi * doy / 365.25),
      cos_year = cos(2 * pi * doy / 365.25),
      summer_school_break = as.integer(
        (month_num == 6 & mday(date) >= 15) |
          month_num %in% c(7, 8) |
          (month_num == 9 & mday(date) <= 15)
      ),
      is_holiday = as.integer(date %in% holiday_dates),
      is_religious_holiday = as.integer(
        date %in% holiday_lookup$date[
          holiday_lookup$holiday_type %in% c("ramadan_bayram", "kurban_bayram")
        ]
      ),
      is_ramadan = as.integer(date %in% ramadan_lookup$date),
      holiday_eve = as.integer((date + days(1)) %in% holiday_dates),
      post_holiday = as.integer((date - days(1)) %in% holiday_dates),
      near_holiday = as.integer(
        date %in% c(
          holiday_dates - days(2), holiday_dates - days(1),
          holiday_dates, holiday_dates + days(1), holiday_dates + days(2)
        )
      ),
      bridge_day = as.integer(
        dow %in% c(1, 5) &
          !(date %in% holiday_dates) &
          ((date - days(1)) %in% holiday_dates | (date + days(1)) %in% holiday_dates)
      )
    ) |>
    dplyr::select(-date)
}

safe_window_mean <- function(Y_mat, start_idx, end_idx) {
  rowMeans(Y_mat[start_idx:end_idx, , drop = FALSE], na.rm = TRUE) |>
    mean(na.rm = TRUE)
}

safe_window_ramp <- function(Y_mat, start_idx, end_idx) {
  mean(Y_mat[start_idx:end_idx, H9] - Y_mat[start_idx:end_idx, H6], na.rm = TRUE)
}

make_load_regime_features <- function(Y_mat, d_index) {
  bind_rows(lapply(d_index, function(d) {
    lag1 <- Y_mat[d - 1, ]
    lag7 <- Y_mat[d - 7, ]
    
    tibble(
      ramp_lag1 = as.numeric(Y_mat[d - 1, H9] - Y_mat[d - 1, H6]),
      ramp_lag7 = as.numeric(Y_mat[d - 7, H9] - Y_mat[d - 7, H6]),
      early_lag1 = mean(Y_mat[d - 1, c(H5, H6, H7)], na.rm = TRUE),
      mean_lag1 = mean(lag1, na.rm = TRUE),
      mean_lag7 = mean(lag7, na.rm = TRUE),
      peak_lag1 = max(lag1, na.rm = TRUE),
      peak_lag7 = max(lag7, na.rm = TRUE),
      min_lag1 = min(lag1, na.rm = TRUE),
      min_lag7 = min(lag7, na.rm = TRUE),
      daily_range_lag1 = max(lag1, na.rm = TRUE) - min(lag1, na.rm = TRUE),
      daily_range_lag7 = max(lag7, na.rm = TRUE) - min(lag7, na.rm = TRUE),
      roll_mean_7 = safe_window_mean(Y_mat, d - 7, d - 1),
      roll_mean_28 = safe_window_mean(Y_mat, d - 28, d - 1),
      roll_ramp_7 = safe_window_ramp(Y_mat, d - 7, d - 1),
      roll_ramp_28 = safe_window_ramp(Y_mat, d - 28, d - 1)
    )
  }))
}

make_features <- function(Y_mat, date_vec, d_index) {
  bind_cols(
    make_load_regime_features(Y_mat, d_index),
    make_calendar_features(date_vec[d_index])
  )
}

standardize_train_new <- function(Z_train, Z_new) {
  mu <- sapply(Z_train, mean, na.rm = TRUE)
  sdv <- sapply(Z_train, sd, na.rm = TRUE)
  
  sdv[is.na(sdv) | sdv == 0] <- 1
  
  Z_train_s <- sweep(as.matrix(Z_train), 2, mu, "-")
  Z_train_s <- sweep(Z_train_s, 2, sdv, "/")
  
  Z_new_s <- sweep(as.matrix(Z_new), 2, mu, "-")
  Z_new_s <- sweep(Z_new_s, 2, sdv, "/")
  
  list(
    Z_train = Z_train_s,
    Z_new = Z_new_s,
    mu = mu,
    sd = sdv
  )
}

clean_coefficients <- function(coef_mat) {
  coef_mat[is.na(coef_mat)] <- 0
  coef_mat
}

feature_indices <- (MIN_LAG_DAYS + 1):n_days

cat("\nPrecomputing load/calendar features...\n")

all_feature_rows <- make_features(Y, dates, feature_indices)

get_feature_rows <- function(index_vec) {
  out <- all_feature_rows[match(index_vec, feature_indices), , drop = FALSE]
  if (any(!complete.cases(out))) {
    stop("Missing precomputed features for requested index.")
  }
  out
}

# ------------------------------------------------------------
# Test sample
# ------------------------------------------------------------

test_idx <- which(dates >= TEST_START)
test_idx <- test_idx[test_idx > max(MIN_LAG_DAYS, MIN_TRAIN_DAYS)]

if (length(test_idx) == 0) {
  first_test <- floor(0.8 * n_days) + 1
  test_idx <- seq(first_test, n_days)
  test_idx <- test_idx[test_idx > max(MIN_LAG_DAYS, MIN_TRAIN_DAYS)]
}

if (is.finite(DEBUG_N)) {
  test_idx <- head(test_idx, DEBUG_N)
}

cat("\nForecast days:", length(test_idx), "\n")
cat("First:", as.character(min(dates[test_idx])), "\n")
cat("Last: ", as.character(max(dates[test_idx])), "\n")

# ------------------------------------------------------------
# One-day forecast
# ------------------------------------------------------------

forecast_one_day <- function(i, Y, dates, hours, K, nbasis, norder) {
  Y_train <- Y[1:(i - 1), , drop = FALSE]
  dates_train <- dates[1:(i - 1)]
  n_train <- nrow(Y_train)
  
  basis <- create.bspline.basis(
    rangeval = range(hours),
    nbasis = nbasis,
    norder = norder
  )
  
  fd_train <- Data2fd(
    argvals = hours,
    y = t(Y_train),
    basisobj = basis
  )
  
  pca <- pca.fd(fd_train, nharm = K)
  scores <- pca$scores[, 1:K, drop = FALSE]
  
  d_index <- (MIN_LAG_DAYS + 1):n_train
  
  Y_score <- scores[d_index, , drop = FALSE]
  X_lag1 <- scores[d_index - 1, , drop = FALSE]
  X_lag7 <- scores[d_index - 7, , drop = FALSE]
  
  Z_train_raw <- get_feature_rows(d_index)
  Z_new_raw <- get_feature_rows(i)
  
  Z_scaled <- standardize_train_new(Z_train_raw, Z_new_raw)
  
  X_score <- cbind(
    intercept = 1,
    X_lag1,
    X_lag7,
    Z_scaled$Z_train
  )
  
  colnames(X_score) <- c(
    "intercept",
    paste0("lag1_score", 1:K),
    paste0("lag7_score", 1:K),
    colnames(Z_train_raw)
  )
  
  mu_hat <- as.vector(eval.fd(hours, pca$meanfd))
  phi_hat <- eval.fd(hours, pca$harmonics)
  
  fit_score <- lm.fit(X_score, Y_score)
  B_hat <- clean_coefficients(fit_score$coefficients)
  
  x_new <- c(
    1,
    scores[n_train, ],
    scores[n_train - 6, ],
    as.numeric(Z_scaled$Z_new)
  )
  
  score_hat <- as.vector(x_new %*% B_hat)
  
  y_base <- as.vector(mu_hat + phi_hat %*% score_hat)
  
  list(
    forecast_base = y_base,
    score_hat = score_hat,
    varprop = pca$varprop[1:K],
    z_new_raw = Z_new_raw
  )
}

# ------------------------------------------------------------
# Rolling forecasts
# ------------------------------------------------------------

pred_base <- matrix(
  NA_real_,
  nrow = length(test_idx),
  ncol = length(hours),
  dimnames = list(as.character(dates[test_idx]), as.character(hours))
)

score_forecasts <- matrix(
  NA_real_,
  nrow = length(test_idx),
  ncol = K,
  dimnames = list(as.character(dates[test_idx]), paste0("score", 1:K))
)

varprops <- matrix(
  NA_real_,
  nrow = length(test_idx),
  ncol = K,
  dimnames = list(as.character(dates[test_idx]), paste0("FPC", 1:K))
)

feature_track <- list()

for (j in seq_along(test_idx)) {
  i <- test_idx[j]
  
  out <- forecast_one_day(
    i = i,
    Y = Y,
    dates = dates,
    hours = hours,
    K = K,
    nbasis = NBASIS,
    norder = NORDER
  )
  
  pred_base[j, ] <- out$forecast_base
  score_forecasts[j, ] <- out$score_hat
  varprops[j, ] <- out$varprop
  feature_track[[j]] <- out$z_new_raw
  
  if (j %% 25 == 0) {
    message("Finished ", j, " / ", length(test_idx),
            " | date: ", dates[i])
  }
}

feature_track_df <- bind_rows(feature_track) |>
  mutate(date = dates[test_idx], .before = 1)

# ------------------------------------------------------------
# Benchmarks
# ------------------------------------------------------------

actual_test <- Y[test_idx, , drop = FALSE]
pred_naive <- Y[test_idx - 7, , drop = FALSE]

if (has_official) {
  pred_official <- Y_official[test_idx, , drop = FALSE]
}

make_eval_df <- function(actual, pred, model_name) {
  tibble(
    date = rep(dates[test_idx], each = length(hours)),
    hour = rep(hours, times = length(test_idx)),
    model = model_name,
    actual = as.vector(t(actual)),
    forecast = as.vector(t(pred))
  ) |>
    mutate(
      year = year(date),
      month = month(date),
      ym = format(date, "%Y-%m"),
      error = actual - forecast,
      abs_error = abs(error),
      sq_error = error^2
    )
}

eval_df <- bind_rows(
  make_eval_df(actual_test, pred_naive, "Seasonal naive: Y[d-7]"),
  make_eval_df(actual_test, pred_base, MODEL_NAME_BASE)
)

if (has_official) {
  eval_df <- bind_rows(
    eval_df,
    make_eval_df(actual_test, pred_official, "Official forecast")
  )
}

# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------

metric_table <- eval_df |>
  group_by(model) |>
  summarise(
    n = sum(!is.na(error)),
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

base_mae <- metric_table$mae[metric_table$model == "Seasonal naive: Y[d-7]"]
base_rmse <- metric_table$rmse[metric_table$model == "Seasonal naive: Y[d-7]"]

official_mae <- metric_table$mae[metric_table$model == "Official forecast"]
official_rmse <- metric_table$rmse[metric_table$model == "Official forecast"]

metric_table <- metric_table |>
  mutate(
    mae_improvement_vs_naive_pct = 100 * (base_mae - mae) / base_mae,
    rmse_improvement_vs_naive_pct = 100 * (base_rmse - rmse) / base_rmse,
    mae_improvement_vs_official_pct = if (length(official_mae) == 1) {
      100 * (official_mae - mae) / official_mae
    } else {
      NA_real_
    },
    rmse_improvement_vs_official_pct = if (length(official_rmse) == 1) {
      100 * (official_rmse - rmse) / official_rmse
    } else {
      NA_real_
    }
  )

print(metric_table)

write_csv(
  metric_table,
  "output/tables/fpca_VAR_1_7_holiday_regime_metrics_overall.csv"
)

metrics_by_hour <- eval_df |>
  group_by(hour, model) |>
  summarise(
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

metrics_by_month <- eval_df |>
  group_by(ym, model) |>
  summarise(
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

metrics_by_year <- eval_df |>
  group_by(year, model) |>
  summarise(
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(metrics_by_hour, "output/tables/fpca_VAR_1_7_holiday_regime_by_hour.csv")
write_csv(metrics_by_month, "output/tables/fpca_VAR_1_7_holiday_regime_by_month.csv")
write_csv(metrics_by_year, "output/tables/fpca_VAR_1_7_holiday_regime_by_year.csv")

write_csv(
  feature_track_df,
  "output/tables/fpca_VAR_1_7_holiday_regime_features_used.csv"
)

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

p_hour <- metrics_by_hour |>
  ggplot(aes(hour, mae, linetype = model)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "MAE",
    title = "Forecast MAE by hour",
    linetype = "Model"
  )

ggsave(
  "output/figs/fpca_VAR_1_7_holiday_regime_mae_by_hour.png",
  p_hour,
  width = 9,
  height = 5
)

p_month <- metrics_by_month |>
  ggplot(aes(ym, mae, group = model, linetype = model)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Month",
    y = "MAE",
    title = "Forecast MAE by month",
    linetype = "Model"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
  "output/figs/fpca_VAR_1_7_holiday_regime_mae_by_month.png",
  p_month,
  width = 10,
  height = 5
)

forecast_objects <- list(
  dates = dates[test_idx],
  hours = hours,
  actual = actual_test,
  holiday_regime = pred_base,
  seasonal_naive = pred_naive,
  official_forecast = if (has_official) pred_official else NULL,
  score_forecasts = score_forecasts,
  varprops = varprops,
  features_used = feature_track_df,
  settings = list(
    K = K,
    NBASIS = NBASIS,
    NORDER = NORDER,
    TEST_START = TEST_START,
    MIN_TRAIN_DAYS = MIN_TRAIN_DAYS,
    MIN_LAG_DAYS = MIN_LAG_DAYS
  ),
  metrics_overall = metric_table,
  metrics_by_hour = metrics_by_hour,
  metrics_by_month = metrics_by_month,
  metrics_by_year = metrics_by_year
)

saveRDS(
  forecast_objects,
  "output/results/fpca_VAR_1_7_holiday_regime_forecast_results.rds"
)

write_csv(
  eval_df,
  "output/tables/fpca_VAR_1_7_holiday_regime_eval_long.csv"
)
