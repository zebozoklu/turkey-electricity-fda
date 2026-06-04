# 16d_forecast_fpca_RF_rolling_subset.R
#
# Rolling-basis RF benchmark — FPCA basis is RE-FIT inside the rolling loop
# on every forecast day, exactly as in script 16, but restricted to a
# manageable sub-period so the run completes in a few hours.
#
# Purpose: check whether the fixed-basis shortcut in script 16b introduces
# any material error. If rolling-RF ≈ fixed-RF, the fixed basis is validated.
# We also compare rolling-RF vs OLS on the same sub-period.
#
# Key differences vs script 16:
#   - TEST_END caps the test window (default: 2023-12-31, ~365 forecast days)
#   - SAMPLE_EVERY_N subsamples test days for an additional speed-up (default 1 = all)
#   - RF_TREES = 100  (vs 500 in script 16)
#   - Output files use the suffix "_rolling" to avoid overwriting script 16b results
#
# Run time with defaults:
#   ~365 test days × rolling FPCA + 100 trees ≈ 2-4 hours (single-threaded)
#   Set SAMPLE_EVERY_N = 3 to cut that to ~45-80 minutes at the cost of ~122 days.

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

# ============================================================
# Settings
# ============================================================

K              <- 2
K_DERIV        <- 3
NBASIS         <- 12
NORDER         <- 4
TEST_START     <- as.Date("2023-01-01")
TEST_END       <- as.Date("2023-12-31")   # <-- narrow window; change to extend
SAMPLE_EVERY_N <- 1                        # 1 = all days; 3 = every 3rd day, etc.
MIN_TRAIN_DAYS <- 365
MIN_LAG_DAYS   <- 28

RF_TREES       <- 100
RF_SEED        <- 42

MODEL_NAME     <- "FPCA RF (rolling basis)"
FIXED_RF_NAME  <- "FPCA RF (fixed basis)"
OLS_NAME       <- "FPCA VAR(1,7) + holiday/load/derivative regime"

# ============================================================
# Load data
# ============================================================

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
H5  <- hcol(5);  H6  <- hcol(6);  H7  <- hcol(7);  H9  <- hcol(9)
H17 <- hcol(17); H20 <- hcol(20); H23 <- hcol(23)

# ============================================================
# Calendar / holiday helpers  (identical to scripts 11 / 16 / 16b)
# ============================================================

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
  rr <- tribble(~start, ~end,
    "2019-05-06","2019-06-03","2020-04-24","2020-05-23","2021-04-13","2021-05-12",
    "2022-04-02","2022-05-01","2023-03-23","2023-04-20","2024-03-11","2024-04-09",
    "2025-03-01","2025-03-29") |> mutate(start = as.Date(start), end = as.Date(end))
  bind_rows(lapply(seq_len(nrow(rr)), function(i)
    tibble(date = date_seq(rr$start[i], rr$end[i])))) |>
    filter(year(date) %in% years) |> distinct(date)
}

holiday_lookup <- make_holiday_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))
ramadan_lookup <- make_ramadan_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))

make_calendar_features <- function(date_vec) {
  doy <- yday(date_vec); dow <- wday(date_vec, week_start = 1)
  mn  <- month(date_vec); hd  <- holiday_lookup$date
  tibble(date = date_vec) |>
    mutate(
      dow = dow, weekend = as.integer(dow >= 6),
      monday = as.integer(dow == 1), friday = as.integer(dow == 5), month = mn,
      sin_week = sin(2 * pi * dow / 7), cos_week = cos(2 * pi * dow / 7),
      sin_year = sin(2 * pi * doy / 365.25), cos_year = cos(2 * pi * doy / 365.25),
      summer_school_break = as.integer(
        (mn == 6 & mday(date) >= 15) | mn %in% c(7, 8) | (mn == 9 & mday(date) <= 15)),
      is_holiday           = as.integer(date %in% hd),
      is_religious_holiday = as.integer(date %in% holiday_lookup$date[
        holiday_lookup$holiday_type %in% c("ramadan_bayram", "kurban_bayram")]),
      is_ramadan   = as.integer(date %in% ramadan_lookup$date),
      holiday_eve  = as.integer((date + days(1)) %in% hd),
      post_holiday = as.integer((date - days(1)) %in% hd),
      near_holiday = as.integer(date %in% c(
        hd - days(2), hd - days(1), hd, hd + days(1), hd + days(2))),
      bridge_day   = as.integer(dow %in% c(1, 5) & !(date %in% hd) &
        ((date - days(1)) %in% hd | (date + days(1)) %in% hd))
    ) |> dplyr::select(-date)
}

safe_window_mean <- function(Ym, s, e) mean(rowMeans(Ym[s:e, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
safe_window_ramp <- function(Ym, s, e) mean(Ym[s:e, H9] - Ym[s:e, H6], na.rm = TRUE)

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
    lag1  <- Y_mat[d - 1, ]; lag7  <- Y_mat[d - 7, ]
    dlag1 <- D_mat[d - 1, ]; dlag7 <- D_mat[d - 7, ]
    tibble(
      ramp_lag1         = as.numeric(Y_mat[d-1, H9] - Y_mat[d-1, H6]),
      ramp_lag7         = as.numeric(Y_mat[d-7, H9] - Y_mat[d-7, H6]),
      early_lag1        = mean(Y_mat[d-1, c(H5, H6, H7)], na.rm = TRUE),
      early_ramp_lag1   = as.numeric(Y_mat[d-1, H7] - Y_mat[d-1, H5]),
      early_ramp_lag7   = as.numeric(Y_mat[d-7, H7] - Y_mat[d-7, H5]),
      evening_ramp_lag1 = as.numeric(Y_mat[d-1, H20] - Y_mat[d-1, H17]),
      evening_ramp_lag7 = as.numeric(Y_mat[d-7, H20] - Y_mat[d-7, H17]),
      day_end_drop_lag1 = as.numeric(Y_mat[d-1, H23] - Y_mat[d-1, H20]),
      day_end_drop_lag7 = as.numeric(Y_mat[d-7, H23] - Y_mat[d-7, H20]),
      mean_lag1         = mean(lag1, na.rm = TRUE),  mean_lag7  = mean(lag7, na.rm = TRUE),
      peak_lag1         = max(lag1,  na.rm = TRUE),  peak_lag7  = max(lag7,  na.rm = TRUE),
      min_lag1          = min(lag1,  na.rm = TRUE),  min_lag7   = min(lag7,  na.rm = TRUE),
      daily_range_lag1  = max(lag1,  na.rm = TRUE) - min(lag1,  na.rm = TRUE),
      daily_range_lag7  = max(lag7,  na.rm = TRUE) - min(lag7,  na.rm = TRUE),
      roll_mean_7       = safe_window_mean(Y_mat, d-7,  d-1),
      roll_mean_28      = safe_window_mean(Y_mat, d-28, d-1),
      roll_ramp_7       = safe_window_ramp(Y_mat, d-7,  d-1),
      roll_ramp_28      = safe_window_ramp(Y_mat, d-28, d-1),
      deriv_sd_lag1     = sd(dlag1,  na.rm = TRUE),  deriv_sd_lag7   = sd(dlag7,  na.rm = TRUE),
      deriv_max_lag1    = max(dlag1, na.rm = TRUE),  deriv_max_lag7  = max(dlag7, na.rm = TRUE),
      deriv_min_lag1    = min(dlag1, na.rm = TRUE),  deriv_min_lag7  = min(dlag7, na.rm = TRUE),
      deriv_range_lag1  = max(dlag1, na.rm = TRUE) - min(dlag1, na.rm = TRUE),
      deriv_range_lag7  = max(dlag7, na.rm = TRUE) - min(dlag7, na.rm = TRUE)
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

# ============================================================
# Precompute external features on full dataset
# (load-regime features only need past lags, no leakage)
# ============================================================

feature_indices  <- (MIN_LAG_DAYS + 1):n_days

cat("Precomputing load/calendar features...\n")
all_feature_rows <- make_features(Y, dates, feature_indices)

get_feature_rows <- function(index_vec) {
  out <- all_feature_rows[match(index_vec, feature_indices), , drop = FALSE]
  if (any(!complete.cases(out))) stop("Missing features for requested indices.")
  out
}

# ============================================================
# Test indices — restricted window + optional subsampling
# ============================================================

test_idx <- which(dates >= TEST_START & dates <= TEST_END)
test_idx <- test_idx[test_idx > max(MIN_LAG_DAYS, MIN_TRAIN_DAYS)]

if (SAMPLE_EVERY_N > 1) {
  test_idx <- test_idx[seq(1, length(test_idx), by = SAMPLE_EVERY_N)]
  cat("Subsampling: keeping every", SAMPLE_EVERY_N, "th day.\n")
}

cat("Forecast days:", length(test_idx), "\n")
cat("First:", as.character(min(dates[test_idx])), "\n")
cat("Last: ", as.character(max(dates[test_idx])), "\n\n")

# ============================================================
# One-day forecast with rolling FPCA refit
# ============================================================

forecast_one_day_rf_rolling <- function(i) {
  Y_train  <- Y[1:(i - 1), , drop = FALSE]
  n_train  <- nrow(Y_train)

  # Refit FPCA on all training data up to i-1
  fd_train    <- make_fd_object(Y_train, hours, NBASIS, NORDER)
  pca         <- pca.fd(fd_train, nharm = K)
  scores      <- pca$scores[, 1:K, drop = FALSE]

  fd_deriv    <- deriv.fd(fd_train, Lfdobj = 1)
  deriv_pca   <- pca.fd(fd_deriv, nharm = K_DERIV)
  deriv_scores <- deriv_pca$scores[, 1:K_DERIV, drop = FALSE]

  d_index <- (MIN_LAG_DAYS + 1):n_train
  Y_score  <- scores[d_index, , drop = FALSE]
  X_lag1   <- scores[d_index - 1, , drop = FALSE]
  X_lag7   <- scores[d_index - 7, , drop = FALSE]
  X_dlag1  <- deriv_scores[d_index - 1, , drop = FALSE]
  X_dlag7  <- deriv_scores[d_index - 7, , drop = FALSE]

  Z_train_raw <- get_feature_rows(d_index)
  Z_new_raw   <- get_feature_rows(i)
  Z_sc        <- standardize_train_new(Z_train_raw, Z_new_raw)

  X_rf <- cbind(X_lag1, X_lag7, X_dlag1, X_dlag7, Z_sc$Z_train)
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
    as.numeric(Z_sc$Z_new)
  )
  x_new_df <- setNames(as.data.frame(t(x_new)), colnames(X_rf))

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

# ============================================================
# Rolling loop
# ============================================================

pred_rf     <- matrix(NA_real_, nrow = length(test_idx), ncol = length(hours))
score_preds <- matrix(NA_real_, nrow = length(test_idx), ncol = K)
varprops    <- matrix(NA_real_, nrow = length(test_idx), ncol = K)

rownames(pred_rf) <- as.character(dates[test_idx])
colnames(pred_rf) <- as.character(hours)

cat("Starting rolling-origin RF forecast (rolling basis)...\n")
t0 <- proc.time()

for (j in seq_along(test_idx)) {
  i              <- test_idx[j]
  out            <- forecast_one_day_rf_rolling(i)
  pred_rf[j, ]      <- out$forecast
  score_preds[j, ]  <- out$score_hat
  varprops[j, ]     <- out$varprop

  if (j %% 25 == 0) {
    el  <- round((proc.time() - t0)["elapsed"])
    rem <- round(el / j * (length(test_idx) - j))
    message("  ", j, " / ", length(test_idx),
            "  date: ", dates[i],
            "  elapsed: ", el, "s",
            "  est remaining: ", rem, "s")
  }
}

cat("Done. Total time:", round((proc.time() - t0)["elapsed"]), "s\n\n")

# ============================================================
# Evaluation
# ============================================================

actual_test <- Y[test_idx, , drop = FALSE]
pred_naive  <- Y[test_idx - 7, , drop = FALSE]
if (has_official) pred_official <- Y_official[test_idx, , drop = FALSE]

# Load OLS errors for the same sub-period
ols_eval <- read_csv(
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  show_col_types = FALSE, col_types = cols(date = col_character())
) |>
  mutate(date = as.Date(date), hour = as.integer(hour)) |>
  filter(model == OLS_NAME, date %in% dates[test_idx])

# Load fixed-basis RF errors for the same sub-period (if available)
fixed_rf_file <- "output/tables/fpca_RF_fixed_eval_long.csv"
if (file.exists(fixed_rf_file)) {
  fixed_rf_eval <- read_csv(fixed_rf_file, show_col_types = FALSE,
                             col_types = cols(date = col_character())) |>
    mutate(date = as.Date(date), hour = as.integer(hour)) |>
    filter(date %in% dates[test_idx])
} else {
  fixed_rf_eval <- NULL
}

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
  make_eval_df(actual_test, pred_rf,    MODEL_NAME),
  make_eval_df(actual_test, pred_naive, "Seasonal naive: Y[d-7]")
)
if (has_official)
  eval_df <- bind_rows(eval_df, make_eval_df(actual_test, pred_official, "Official forecast"))

eval_combined <- bind_rows(
  eval_df,
  ols_eval   |> dplyr::select(all_of(names(eval_df))),
  if (!is.null(fixed_rf_eval)) fixed_rf_eval |> dplyr::select(all_of(names(eval_df))) else NULL
)

# ============================================================
# Metrics
# ============================================================

metrics <- eval_combined |>
  filter(actual > 0) |>
  group_by(model) |>
  summarise(
    n    = n(),
    mae  = mean(abs_error,          na.rm = TRUE),
    rmse = sqrt(mean(sq_error,      na.rm = TRUE)),
    mape = mean(abs_error / actual, na.rm = TRUE) * 100,
    bias = mean(error,              na.rm = TRUE),
    .groups = "drop"
  )

off <- metrics |> filter(model == "Official forecast")
if (nrow(off) == 1) {
  metrics <- metrics |>
    mutate(
      mae_gain  = 100 * (off$mae  - mae)  / off$mae,
      rmse_gain = 100 * (off$rmse - rmse) / off$rmse,
      mape_gain = 100 * (off$mape - mape) / off$mape
    )
}

cat("=== Overall metrics (sub-period:", as.character(TEST_START), "to",
    as.character(TEST_END), ") ===\n")
print(metrics |>
        dplyr::select(model, mae, rmse, mape, bias, mae_gain, mape_gain) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))),
      width = 140)

write_csv(metrics, "output/tables/fpca_RF_rolling_metrics_overall.csv")

# By-hour
metrics_by_hour <- eval_combined |>
  group_by(hour, model) |>
  summarise(mae  = mean(abs_error, na.rm = TRUE),
            mape = mean(abs_error / actual, na.rm = TRUE) * 100,
            .groups = "drop")

write_csv(metrics_by_hour, "output/tables/fpca_RF_rolling_by_hour.csv")

# Rolling vs Fixed RF gain by hour
rolling_h <- metrics_by_hour |> filter(model == MODEL_NAME)     |> dplyr::select(hour, rolling_mae = mae)
fixed_h   <- metrics_by_hour |> filter(model == FIXED_RF_NAME)  |> dplyr::select(hour, fixed_mae   = mae)
ols_h     <- metrics_by_hour |> filter(model == OLS_NAME)       |> dplyr::select(hour, ols_mae     = mae)

if (nrow(fixed_h) > 0) {
  gain_rf <- inner_join(rolling_h, fixed_h, by = "hour") |>
    mutate(gain_rolling_over_fixed = 100 * (fixed_mae - rolling_mae) / fixed_mae)
  write_csv(gain_rf, "output/tables/fpca_RF_rolling_vs_fixed_by_hour.csv")
}

# ============================================================
# Plots
# ============================================================

models_to_plot <- intersect(
  c(MODEL_NAME, FIXED_RF_NAME, OLS_NAME, "Official forecast"),
  unique(metrics_by_hour$model)
)

p_hour <- metrics_by_hour |>
  filter(model %in% models_to_plot) |>
  ggplot(aes(hour, mae, colour = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  labs(title = paste0("MAE by hour (", TEST_START, " – ", TEST_END, ")"),
       subtitle = "Rolling-basis RF vs fixed-basis RF vs OLS vs Official",
       x = "Hour", y = "MAE", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

ggsave("output/figs/fpca_RF_rolling_mae_by_hour.png", p_hour, width = 9, height = 5, dpi = 150)
cat("Saved: output/figs/fpca_RF_rolling_mae_by_hour.png\n")

if (nrow(fixed_h) > 0) {
  p_vs_fixed <- ggplot(gain_rf, aes(x = hour, y = gain_rolling_over_fixed)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_col(aes(fill = gain_rolling_over_fixed > 0), show.legend = FALSE, alpha = 0.8) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
    scale_x_continuous(breaks = seq(0, 23, 3)) +
    labs(title = "Rolling-basis RF vs Fixed-basis RF by hour",
         subtitle = "Blue = rolling better, Red = fixed better",
         x = "Hour", y = "MAE gain (%)") +
    theme_bw(base_size = 11)

  ggsave("output/figs/fpca_RF_rolling_vs_fixed_by_hour.png", p_vs_fixed, width = 8, height = 4, dpi = 150)
  cat("Saved: output/figs/fpca_RF_rolling_vs_fixed_by_hour.png\n")
}

# Save eval for DM test
write_csv(eval_df |> filter(model == MODEL_NAME),
          "output/tables/fpca_RF_rolling_eval_long.csv")

cat("\nDone.\n")
