# 16_forecast_fpca_RF.R
#
# Random Forest benchmark for the FPCA score regression.
# Identical rolling-origin design and feature matrix as script 11, but
# replaces OLS (lm.fit) with ranger on each score component separately.
#
# Purpose: test whether the linear score equation misses nonlinear structure
# in z_d. If FPCA-OLS matches or beats FPCA-RF, the linear spec is adequate.

library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)
library(lubridate)
library(readr)
library(ranger)

dir.create("output/figs",    recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Settings  (must match script 11 exactly)
# ------------------------------------------------------------

K            <- 2
K_DERIV      <- 3
NBASIS       <- 12
NORDER       <- 4
TEST_START   <- as.Date("2023-01-01")
MIN_TRAIN_DAYS <- 365
MIN_LAG_DAYS   <- 28
DEBUG_N        <- Inf

RF_TREES     <- 500
RF_SEED      <- 42

MODEL_NAME   <- "FPCA RF score regression"
OLS_NAME     <- "FPCA VAR(1,7) + holiday/load/derivative regime"

# ------------------------------------------------------------
# Load curves  (identical to script 11)
# ------------------------------------------------------------

curves <- readRDS("data/processed/load_curves.rds")

dates <- as.Date(curves$dates)
hours <- curves$hours
Y     <- curves$actual
colnames(Y) <- as.character(hours)

has_official <- !is.null(curves$forecast)
if (has_official) {
  Y_official <- curves$forecast
  colnames(Y_official) <- as.character(hours)
}

ok    <- complete.cases(Y)
dates <- dates[ok]
Y     <- Y[ok, , drop = FALSE]
if (has_official) Y_official <- Y_official[ok, , drop = FALSE]

n_days <- nrow(Y)

hcol <- function(h) { out <- match(h, hours); if (is.na(out)) stop("Hour ", h); out }
H5 <- hcol(5); H6 <- hcol(6); H7 <- hcol(7); H9 <- hcol(9)
H17 <- hcol(17); H20 <- hcol(20); H23 <- hcol(23)

# ------------------------------------------------------------
# Calendar / holiday helpers  (identical to script 11)
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
    "2019-06-04","2019-06-06","ramadan_bayram","2020-05-24","2020-05-26","ramadan_bayram",
    "2021-05-13","2021-05-15","ramadan_bayram","2022-05-02","2022-05-04","ramadan_bayram",
    "2023-04-21","2023-04-23","ramadan_bayram","2024-04-10","2024-04-12","ramadan_bayram",
    "2025-03-30","2025-04-01","ramadan_bayram",
    "2019-08-11","2019-08-14","kurban_bayram","2020-07-31","2020-08-03","kurban_bayram",
    "2021-07-20","2021-07-23","kurban_bayram","2022-07-09","2022-07-12","kurban_bayram",
    "2023-06-28","2023-07-01","kurban_bayram","2024-06-16","2024-06-19","kurban_bayram",
    "2025-06-06","2025-06-09","kurban_bayram"
  ) |> mutate(start = as.Date(start), end = as.Date(end))
  religious_holidays <- bind_rows(lapply(seq_len(nrow(religious_ranges)), function(i) {
    tibble(date = date_seq(religious_ranges$start[i], religious_ranges$end[i]),
           holiday_type = religious_ranges$holiday_type[i])
  }))
  bind_rows(fixed_holidays, religious_holidays) |>
    filter(year(date) %in% years) |> distinct(date, .keep_all = TRUE)
}

make_ramadan_lookup <- function(years) {
  ramadan_ranges <- tribble(
    ~start, ~end,
    "2019-05-06","2019-06-03","2020-04-24","2020-05-23","2021-04-13","2021-05-12",
    "2022-04-02","2022-05-01","2023-03-23","2023-04-20","2024-03-11","2024-04-09",
    "2025-03-01","2025-03-29"
  ) |> mutate(start = as.Date(start), end = as.Date(end))
  bind_rows(lapply(seq_len(nrow(ramadan_ranges)), function(i) {
    tibble(date = date_seq(ramadan_ranges$start[i], ramadan_ranges$end[i]))
  })) |> filter(year(date) %in% years) |> distinct(date)
}

holiday_lookup <- make_holiday_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))
ramadan_lookup <- make_ramadan_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))

make_calendar_features <- function(date_vec) {
  doy <- yday(date_vec); dow <- wday(date_vec, week_start = 1)
  month_num <- month(date_vec); holiday_dates <- holiday_lookup$date
  tibble(date = date_vec) |>
    mutate(
      dow = dow, weekend = as.integer(dow >= 6),
      monday = as.integer(dow == 1), friday = as.integer(dow == 5),
      month = month_num,
      sin_week = sin(2 * pi * dow / 7), cos_week = cos(2 * pi * dow / 7),
      sin_year = sin(2 * pi * doy / 365.25), cos_year = cos(2 * pi * doy / 365.25),
      summer_school_break = as.integer(
        (month_num == 6 & mday(date) >= 15) | month_num %in% c(7, 8) |
          (month_num == 9 & mday(date) <= 15)),
      is_holiday = as.integer(date %in% holiday_dates),
      is_religious_holiday = as.integer(date %in% holiday_lookup$date[
        holiday_lookup$holiday_type %in% c("ramadan_bayram", "kurban_bayram")]),
      is_ramadan = as.integer(date %in% ramadan_lookup$date),
      holiday_eve  = as.integer((date + days(1)) %in% holiday_dates),
      post_holiday = as.integer((date - days(1)) %in% holiday_dates),
      near_holiday = as.integer(date %in% c(
        holiday_dates - days(2), holiday_dates - days(1), holiday_dates,
        holiday_dates + days(1), holiday_dates + days(2))),
      bridge_day = as.integer(
        dow %in% c(1, 5) & !(date %in% holiday_dates) &
          ((date - days(1)) %in% holiday_dates | (date + days(1)) %in% holiday_dates))
    ) |> dplyr::select(-date)
}

safe_window_mean <- function(Y_mat, s, e) mean(rowMeans(Y_mat[s:e, , drop=FALSE], na.rm=TRUE), na.rm=TRUE)
safe_window_ramp <- function(Y_mat, s, e) mean(Y_mat[s:e, H9] - Y_mat[s:e, H6], na.rm=TRUE)

make_fd_object <- function(Y_mat, argvals, nbasis, norder) {
  basis <- create.bspline.basis(rangeval = range(argvals), nbasis = nbasis, norder = norder)
  Data2fd(argvals = argvals, y = t(Y_mat), basisobj = basis)
}

make_derivative_matrix <- function(Y_mat) {
  fd_obj   <- make_fd_object(Y_mat, hours, NBASIS, NORDER)
  deriv_fd <- deriv.fd(fd_obj, Lfdobj = 1)
  out      <- t(eval.fd(hours, deriv_fd))
  colnames(out) <- as.character(hours); out
}

make_load_regime_features <- function(Y_mat, d_index) {
  D_mat <- make_derivative_matrix(Y_mat)
  bind_rows(lapply(d_index, function(d) {
    lag1 <- Y_mat[d-1,]; lag7 <- Y_mat[d-7,]
    dlag1 <- D_mat[d-1,]; dlag7 <- D_mat[d-7,]
    tibble(
      ramp_lag1=as.numeric(Y_mat[d-1,H9]-Y_mat[d-1,H6]),
      ramp_lag7=as.numeric(Y_mat[d-7,H9]-Y_mat[d-7,H6]),
      early_lag1=mean(Y_mat[d-1,c(H5,H6,H7)],na.rm=TRUE),
      early_ramp_lag1=as.numeric(Y_mat[d-1,H7]-Y_mat[d-1,H5]),
      early_ramp_lag7=as.numeric(Y_mat[d-7,H7]-Y_mat[d-7,H5]),
      evening_ramp_lag1=as.numeric(Y_mat[d-1,H20]-Y_mat[d-1,H17]),
      evening_ramp_lag7=as.numeric(Y_mat[d-7,H20]-Y_mat[d-7,H17]),
      day_end_drop_lag1=as.numeric(Y_mat[d-1,H23]-Y_mat[d-1,H20]),
      day_end_drop_lag7=as.numeric(Y_mat[d-7,H23]-Y_mat[d-7,H20]),
      mean_lag1=mean(lag1,na.rm=TRUE), mean_lag7=mean(lag7,na.rm=TRUE),
      peak_lag1=max(lag1,na.rm=TRUE),  peak_lag7=max(lag7,na.rm=TRUE),
      min_lag1=min(lag1,na.rm=TRUE),   min_lag7=min(lag7,na.rm=TRUE),
      daily_range_lag1=max(lag1,na.rm=TRUE)-min(lag1,na.rm=TRUE),
      daily_range_lag7=max(lag7,na.rm=TRUE)-min(lag7,na.rm=TRUE),
      roll_mean_7=safe_window_mean(Y_mat,d-7,d-1),
      roll_mean_28=safe_window_mean(Y_mat,d-28,d-1),
      roll_ramp_7=safe_window_ramp(Y_mat,d-7,d-1),
      roll_ramp_28=safe_window_ramp(Y_mat,d-28,d-1),
      deriv_sd_lag1=sd(dlag1,na.rm=TRUE),   deriv_sd_lag7=sd(dlag7,na.rm=TRUE),
      deriv_max_lag1=max(dlag1,na.rm=TRUE),  deriv_max_lag7=max(dlag7,na.rm=TRUE),
      deriv_min_lag1=min(dlag1,na.rm=TRUE),  deriv_min_lag7=min(dlag7,na.rm=TRUE),
      deriv_range_lag1=max(dlag1,na.rm=TRUE)-min(dlag1,na.rm=TRUE),
      deriv_range_lag7=max(dlag7,na.rm=TRUE)-min(dlag7,na.rm=TRUE)
    )
  }))
}

make_features <- function(Y_mat, date_vec, d_index) {
  bind_cols(make_load_regime_features(Y_mat, d_index),
            make_calendar_features(date_vec[d_index]))
}

standardize_train_new <- function(Z_train, Z_new) {
  mu  <- sapply(Z_train, mean, na.rm = TRUE)
  sdv <- sapply(Z_train, sd,   na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(
    Z_train = sweep(sweep(as.matrix(Z_train), 2, mu, "-"), 2, sdv, "/"),
    Z_new   = sweep(sweep(as.matrix(Z_new),   2, mu, "-"), 2, sdv, "/")
  )
}

# ------------------------------------------------------------
# Precompute features  (identical to script 11)
# ------------------------------------------------------------

feature_indices <- (MIN_LAG_DAYS + 1):n_days

cat("\nPrecomputing load/calendar features...\n")
all_feature_rows <- make_features(Y, dates, feature_indices)

get_feature_rows <- function(index_vec) {
  out <- all_feature_rows[match(index_vec, feature_indices), , drop = FALSE]
  if (any(!complete.cases(out))) stop("Missing features.")
  out
}

# ------------------------------------------------------------
# Test sample
# ------------------------------------------------------------

test_idx <- which(dates >= TEST_START)
test_idx <- test_idx[test_idx > max(MIN_LAG_DAYS, MIN_TRAIN_DAYS)]
if (is.finite(DEBUG_N)) test_idx <- head(test_idx, DEBUG_N)

cat("Forecast days:", length(test_idx), "\n")
cat("First:", as.character(min(dates[test_idx])), "\n")
cat("Last: ", as.character(max(dates[test_idx])), "\n\n")

# ------------------------------------------------------------
# One-day forecast (RF replaces lm.fit)
# ------------------------------------------------------------

forecast_one_day_rf <- function(i) {
  Y_train     <- Y[1:(i-1), , drop = FALSE]
  n_train     <- nrow(Y_train)

  fd_train    <- make_fd_object(Y_train, hours, NBASIS, NORDER)
  pca         <- pca.fd(fd_train, nharm = K)
  scores      <- pca$scores[, 1:K, drop = FALSE]

  fd_deriv    <- deriv.fd(fd_train, Lfdobj = 1)
  deriv_pca   <- pca.fd(fd_deriv, nharm = K_DERIV)
  deriv_scores <- deriv_pca$scores[, 1:K_DERIV, drop = FALSE]

  d_index     <- (MIN_LAG_DAYS + 1):n_train
  Y_score     <- scores[d_index, , drop = FALSE]
  X_lag1      <- scores[d_index - 1, , drop = FALSE]
  X_lag7      <- scores[d_index - 7, , drop = FALSE]
  X_dlag1     <- deriv_scores[d_index - 1, , drop = FALSE]
  X_dlag7     <- deriv_scores[d_index - 7, , drop = FALSE]

  Z_train_raw <- get_feature_rows(d_index)
  Z_new_raw   <- get_feature_rows(i)
  Z_scaled    <- standardize_train_new(Z_train_raw, Z_new_raw)

  # RF feature matrix: lags + external features (no intercept needed)
  X_rf <- cbind(X_lag1, X_lag7, X_dlag1, X_dlag7, Z_scaled$Z_train)
  colnames(X_rf) <- c(
    paste0("lag1_s",  1:K),
    paste0("lag7_s",  1:K),
    paste0("lag1_ds", 1:K_DERIV),
    paste0("lag7_ds", 1:K_DERIV),
    colnames(Z_train_raw)
  )

  x_new <- c(
    scores[n_train, ],
    scores[n_train - 6, ],
    deriv_scores[n_train, ],
    deriv_scores[n_train - 6, ],
    as.numeric(Z_scaled$Z_new)
  )
  x_new_df <- as.data.frame(t(x_new))
  colnames(x_new_df) <- colnames(X_rf)

  # Fit RF per score component
  score_hat <- numeric(K)
  for (k in seq_len(K)) {
    df_train <- as.data.frame(cbind(y = Y_score[, k], X_rf))
    rf_fit   <- ranger(y ~ ., data = df_train, num.trees = RF_TREES,
                       seed = RF_SEED, num.threads = 1, verbose = FALSE)
    score_hat[k] <- predict(rf_fit, data = x_new_df)$predictions
  }

  mu_hat  <- as.vector(eval.fd(hours, pca$meanfd))
  phi_hat <- eval.fd(hours, pca$harmonics)
  y_hat   <- as.vector(mu_hat + phi_hat %*% score_hat)

  list(forecast = y_hat, score_hat = score_hat, varprop = pca$varprop[1:K])
}

# ------------------------------------------------------------
# Rolling loop
# ------------------------------------------------------------

pred_rf     <- matrix(NA_real_, nrow = length(test_idx), ncol = length(hours))
score_preds <- matrix(NA_real_, nrow = length(test_idx), ncol = K)
varprops    <- matrix(NA_real_, nrow = length(test_idx), ncol = K)

cat("Starting rolling-origin RF forecast...\n")
t0 <- proc.time()

for (j in seq_along(test_idx)) {
  i             <- test_idx[j]
  out           <- forecast_one_day_rf(i)
  pred_rf[j, ]     <- out$forecast
  score_preds[j, ] <- out$score_hat
  varprops[j, ]    <- out$varprop

  if (j %% 25 == 0) {
    elapsed <- round((proc.time() - t0)["elapsed"])
    message("  ", j, " / ", length(test_idx),
            "  date: ", dates[i],
            "  elapsed: ", elapsed, "s",
            "  est remaining: ",
            round(elapsed / j * (length(test_idx) - j)), "s")
  }
}

cat("Done. Total time:", round((proc.time() - t0)["elapsed"]), "s\n\n")

rownames(pred_rf) <- as.character(dates[test_idx])
colnames(pred_rf) <- as.character(hours)

# ------------------------------------------------------------
# Evaluation
# ------------------------------------------------------------

actual_test  <- Y[test_idx, , drop = FALSE]
pred_naive   <- Y[test_idx - 7, , drop = FALSE]
if (has_official) pred_official <- Y_official[test_idx, , drop = FALSE]

# Also load OLS results for direct comparison
ols_eval <- read_csv(
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  show_col_types = FALSE,
  col_types = cols(date = col_character())
) |>
  mutate(date = as.Date(date), hour = as.integer(hour)) |>
  filter(model == OLS_NAME)

make_eval_df <- function(actual, pred, model_name) {
  tibble(
    date     = rep(dates[test_idx], each = length(hours)),
    hour     = rep(hours, times = length(test_idx)),
    model    = model_name,
    actual   = as.vector(t(actual)),
    forecast = as.vector(t(pred))
  ) |>
    mutate(
      year = year(date), month = month(date), ym = format(date, "%Y-%m"),
      error = actual - forecast, abs_error = abs(error), sq_error = error^2
    )
}

eval_df <- bind_rows(
  make_eval_df(actual_test, pred_rf,     MODEL_NAME),
  make_eval_df(actual_test, pred_naive,  "Seasonal naive: Y[d-7]")
)
if (has_official)
  eval_df <- bind_rows(eval_df, make_eval_df(actual_test, pred_official, "Official forecast"))

# Append OLS for a unified comparison table
eval_combined <- bind_rows(eval_df, ols_eval |> dplyr::select(all_of(names(eval_df))))

# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------

summarise_metrics <- function(df) {
  df |>
    group_by(model) |>
    summarise(n = n(), mae = mean(abs_error, na.rm=TRUE),
              rmse = sqrt(mean(sq_error, na.rm=TRUE)),
              bias = mean(error, na.rm=TRUE), .groups = "drop")
}

metrics <- summarise_metrics(eval_combined)

official_mae  <- metrics$mae[metrics$model == "Official forecast"]
official_rmse <- metrics$rmse[metrics$model == "Official forecast"]

metrics <- metrics |>
  mutate(
    mae_gain  = if (length(official_mae)  == 1) 100 * (official_mae  - mae)  / official_mae  else NA,
    rmse_gain = if (length(official_rmse) == 1) 100 * (official_rmse - rmse) / official_rmse else NA
  )

cat("=== Overall metrics ===\n")
print(metrics |> select(model, mae, rmse, bias, mae_gain, rmse_gain), width=120)

write_csv(metrics, "output/tables/fpca_RF_metrics_overall.csv")

# By hour
metrics_by_hour <- eval_combined |>
  group_by(hour, model) |>
  summarise(mae = mean(abs_error, na.rm=TRUE),
            rmse = sqrt(mean(sq_error, na.rm=TRUE)), .groups = "drop")

write_csv(metrics_by_hour, "output/tables/fpca_RF_by_hour.csv")

# RF vs OLS gain by hour
rf_hour  <- metrics_by_hour |> filter(model == MODEL_NAME) |> select(hour, rf_mae = mae)
ols_hour <- metrics_by_hour |> filter(model == OLS_NAME)   |> select(hour, ols_mae = mae)

gain_by_hour <- inner_join(rf_hour, ols_hour, by = "hour") |>
  mutate(gain_rf_over_ols = 100 * (ols_mae - rf_mae) / ols_mae)

write_csv(gain_by_hour, "output/tables/fpca_RF_gain_vs_OLS_by_hour.csv")

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

# MAE by hour: RF vs OLS vs official
p_hour <- metrics_by_hour |>
  filter(model %in% c(MODEL_NAME, OLS_NAME, "Official forecast")) |>
  ggplot(aes(hour, mae, colour = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  labs(title = "MAE by hour: RF vs OLS vs Official",
       x = "Hour", y = "MAE", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8))

ggsave("output/figs/fpca_RF_mae_by_hour.png", p_hour, width = 9, height = 5, dpi = 150)
cat("Saved: output/figs/fpca_RF_mae_by_hour.png\n")

# RF gain over OLS by hour
p_gain <- ggplot(gain_by_hour, aes(x = hour, y = gain_rf_over_ols)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_col(fill = "steelblue", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  labs(title = "RF gain over OLS by hour (positive = RF better)",
       x = "Hour", y = "MAE gain (%)") +
  theme_bw(base_size = 11)

ggsave("output/figs/fpca_RF_gain_vs_OLS_by_hour.png", p_gain, width = 8, height = 4, dpi = 150)
cat("Saved: output/figs/fpca_RF_gain_vs_OLS_by_hour.png\n")

# Save eval long for potential DM test use
write_csv(eval_df |> filter(model == MODEL_NAME),
          "output/tables/fpca_RF_eval_long.csv")

cat("\nDone.\n")
