# R/10_forecast_fpca_VAR_1_7_ramp_calendar_local_correction.R

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

K <- 4
NBASIS <- 12
NORDER <- 4
TEST_START <- as.Date("2023-01-01")
MIN_TRAIN_DAYS <- 365
DEBUG_N <- Inf

MODEL_NAME <- "FPCA VAR(1,7) + ramp/calendar + local correction"
BASE_MODEL_NAME <- "FPCA VAR(1,7) + ramp/calendar"

CORRECT_HOURS <- c(7, 8, 9)

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
CORRECT_COLS <- sapply(CORRECT_HOURS, hcol)

# ------------------------------------------------------------
# Feature functions
# ------------------------------------------------------------

make_calendar_features <- function(date_vec) {
  doy <- yday(date_vec)
  
  tibble(
    weekend = as.integer(wday(date_vec, week_start = 1) >= 6),
    sin_year = sin(2 * pi * doy / 365.25),
    cos_year = cos(2 * pi * doy / 365.25)
  )
}

make_ramp_calendar_features <- function(Y_mat, date_vec, d_index) {
  
  tibble(
    ramp_lag1 = as.numeric(Y_mat[d_index - 1, H9] - Y_mat[d_index - 1, H6]),
    ramp_lag7 = as.numeric(Y_mat[d_index - 7, H9] - Y_mat[d_index - 7, H6]),
    early_lag1 = rowMeans(
      Y_mat[d_index - 1, c(H5, H6, H7), drop = FALSE],
      na.rm = TRUE
    )
  ) |>
    bind_cols(make_calendar_features(date_vec[d_index]))
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

# ------------------------------------------------------------
# Test sample
# ------------------------------------------------------------

test_idx <- which(dates >= TEST_START)
test_idx <- test_idx[test_idx > max(7, MIN_TRAIN_DAYS)]

if (length(test_idx) == 0) {
  first_test <- floor(0.8 * n_days) + 1
  test_idx <- seq(first_test, n_days)
  test_idx <- test_idx[test_idx > max(7, MIN_TRAIN_DAYS)]
}

if (is.finite(DEBUG_N)) {
  test_idx <- head(test_idx, DEBUG_N)
}

cat("\nForecast days:", length(test_idx), "\n")
cat("First:", as.character(min(dates[test_idx])), "\n")
cat("Last: ", as.character(max(dates[test_idx])), "\n")

# ------------------------------------------------------------
# One-day forecast with local ramp-hour residual correction
# ------------------------------------------------------------

forecast_one_day_local_correction <- function(i, Y, dates, hours, K, nbasis, norder) {
  
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
  
  d_index <- 8:n_train
  
  Y_score <- scores[d_index, , drop = FALSE]
  
  X_lag1 <- scores[d_index - 1, , drop = FALSE]
  X_lag7 <- scores[d_index - 7, , drop = FALSE]
  
  Z_train_raw <- make_ramp_calendar_features(
    Y_mat = Y_train,
    date_vec = dates_train,
    d_index = d_index
  )
  
  Z_new_raw <- tibble(
    ramp_lag1 = as.numeric(Y_train[n_train, H9] - Y_train[n_train, H6]),
    ramp_lag7 = as.numeric(Y_train[n_train - 6, H9] - Y_train[n_train - 6, H6]),
    early_lag1 = mean(Y_train[n_train, c(H5, H6, H7)], na.rm = TRUE)
  ) |>
    bind_cols(make_calendar_features(dates[i]))
  
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
  
  fit_score <- lm.fit(X_score, Y_score)
  B_hat <- fit_score$coefficients
  
  x_new <- c(
    1,
    scores[n_train, ],
    scores[n_train - 6, ],
    as.numeric(Z_scaled$Z_new)
  )
  
  score_hat <- as.vector(x_new %*% B_hat)
  
  mu_hat <- as.vector(eval.fd(hours, pca$meanfd))
  phi_hat <- eval.fd(hours, pca$harmonics)
  
  y_base <- as.vector(mu_hat + phi_hat %*% score_hat)
  
  # ----------------------------------------------------------
  # In-sample fitted base curves for residual correction
  # ----------------------------------------------------------
  
  score_fitted <- X_score %*% B_hat
  
  fitted_curves <- score_fitted %*% t(phi_hat)
  fitted_curves <- sweep(fitted_curves, 2, mu_hat, "+")
  
  actual_train_curves <- Y_train[d_index, , drop = FALSE]
  residual_train_curves <- actual_train_curves - fitted_curves
  
  # residual correction model:
  # residual_{d,h} = a_h + b_h' z_d + noise
  X_resid <- cbind(
    intercept = 1,
    Z_scaled$Z_train
  )
  
  x_resid_new <- c(
    1,
    as.numeric(Z_scaled$Z_new)
  )
  
  y_corrected <- y_base
  correction_hat <- rep(0, length(hours))
  names(correction_hat) <- as.character(hours)
  
  for (h in CORRECT_HOURS) {
    
    hc <- hcol(h)
    
    resid_h <- residual_train_curves[, hc]
    
    fit_resid_h <- lm.fit(X_resid, resid_h)
    
    corr_h <- as.numeric(x_resid_new %*% fit_resid_h$coefficients)
    
    y_corrected[hc] <- y_corrected[hc] + corr_h
    correction_hat[as.character(h)] <- corr_h
  }
  
  list(
    forecast_base = y_base,
    forecast_corrected = y_corrected,
    score_hat = score_hat,
    correction_hat = correction_hat,
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

pred_corrected <- matrix(
  NA_real_,
  nrow = length(test_idx),
  ncol = length(hours),
  dimnames = list(as.character(dates[test_idx]), as.character(hours))
)

correction_mat <- matrix(
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
  
  out <- forecast_one_day_local_correction(
    i = i,
    Y = Y,
    dates = dates,
    hours = hours,
    K = K,
    nbasis = NBASIS,
    norder = NORDER
  )
  
  pred_base[j, ] <- out$forecast_base
  pred_corrected[j, ] <- out$forecast_corrected
  correction_mat[j, ] <- out$correction_hat
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
  make_eval_df(actual_test, pred_base, BASE_MODEL_NAME),
  make_eval_df(actual_test, pred_corrected, MODEL_NAME)
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
    mae_improvement_vs_official_pct = 100 * (official_mae - mae) / official_mae,
    rmse_improvement_vs_official_pct = 100 * (official_rmse - rmse) / official_rmse
  )

print(metric_table)

write_csv(
  metric_table,
  "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_metrics_overall.csv"
)

# ------------------------------------------------------------
# By hour/month/year
# ------------------------------------------------------------

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

write_csv(metrics_by_hour, "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_by_hour.csv")
write_csv(metrics_by_month, "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_by_month.csv")
write_csv(metrics_by_year, "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_by_year.csv")

cat("\nRamp hours h = 7,8,9:\n")
print(metrics_by_hour |> filter(hour %in% CORRECT_HOURS))

# ------------------------------------------------------------
# Gain relative to base ramp/calendar
# ------------------------------------------------------------

comparison <- eval_df |>
  filter(model %in% c(BASE_MODEL_NAME, MODEL_NAME)) |>
  dplyr::select(date, hour, model, abs_error, sq_error) |>
  pivot_wider(
    names_from = model,
    values_from = c(abs_error, sq_error)
  )

names(comparison) <- make.names(names(comparison))

abs_base_col <- make.names(paste0("abs_error_", BASE_MODEL_NAME))
abs_new_col  <- make.names(paste0("abs_error_", MODEL_NAME))
sq_base_col  <- make.names(paste0("sq_error_", BASE_MODEL_NAME))
sq_new_col   <- make.names(paste0("sq_error_", MODEL_NAME))

comparison <- comparison |>
  mutate(
    mae_gain_from_local_correction = .data[[abs_base_col]] - .data[[abs_new_col]],
    mse_gain_from_local_correction = .data[[sq_base_col]] - .data[[sq_new_col]]
  )

gain_overall <- comparison |>
  summarise(
    mae_gain_from_local_correction = mean(mae_gain_from_local_correction, na.rm = TRUE),
    mse_gain_from_local_correction = mean(mse_gain_from_local_correction, na.rm = TRUE)
  )

gain_by_hour <- comparison |>
  group_by(hour) |>
  summarise(
    mae_gain_from_local_correction = mean(mae_gain_from_local_correction, na.rm = TRUE),
    mse_gain_from_local_correction = mean(mse_gain_from_local_correction, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nGain from local correction relative to ramp/calendar:\n")
print(gain_overall)

cat("\nGain at corrected hours:\n")
print(gain_by_hour |> filter(hour %in% CORRECT_HOURS))

write_csv(gain_overall, "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_gain_overall.csv")
write_csv(gain_by_hour, "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_gain_by_hour.csv")

# ------------------------------------------------------------
# Save corrections/features
# ------------------------------------------------------------

correction_df <- as_tibble(correction_mat) |>
  mutate(date = dates[test_idx], .before = 1)

write_csv(
  correction_df,
  "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_values.csv"
)

write_csv(
  feature_track_df,
  "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_features_used.csv"
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
  "output/figs/fpca_VAR_1_7_ramp_calendar_local_correction_mae_by_hour.png",
  p_hour,
  width = 9,
  height = 5
)

p_gain <- gain_by_hour |>
  ggplot(aes(hour, mae_gain_from_local_correction)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "MAE gain from local correction",
    title = "Positive values mean local correction improves ramp/calendar model"
  )

ggsave(
  "output/figs/fpca_VAR_1_7_ramp_calendar_local_correction_gain_by_hour.png",
  p_gain,
  width = 8,
  height = 4.8
)

# ------------------------------------------------------------
# Save result object
# ------------------------------------------------------------

forecast_objects <- list(
  dates = dates[test_idx],
  hours = hours,
  actual = actual_test,
  base_ramp_calendar = pred_base,
  corrected_forecast = pred_corrected,
  seasonal_naive = pred_naive,
  official_forecast = if (has_official) pred_official else NULL,
  correction = correction_mat,
  score_forecasts = score_forecasts,
  varprops = varprops,
  features_used = feature_track_df,
  settings = list(
    K = K,
    NBASIS = NBASIS,
    NORDER = NORDER,
    TEST_START = TEST_START,
    MIN_TRAIN_DAYS = MIN_TRAIN_DAYS,
    CORRECT_HOURS = CORRECT_HOURS
  ),
  metrics_overall = metric_table,
  metrics_by_hour = metrics_by_hour,
  metrics_by_month = metrics_by_month,
  metrics_by_year = metrics_by_year,
  gain_overall = gain_overall,
  gain_by_hour = gain_by_hour
)

saveRDS(
  forecast_objects,
  "output/results/fpca_VAR_1_7_ramp_calendar_local_correction_forecast_results.rds"
)

write_csv(
  eval_df,
  "output/tables/fpca_VAR_1_7_ramp_calendar_local_correction_eval_long.csv"
)


