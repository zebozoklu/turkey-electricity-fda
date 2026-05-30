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
TLF_Q <- 0.75
DEBUG_N <- Inf

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

curves <- readRDS("data/processed/load_curves.rds")
tlf_daily <- readRDS("data/processed/tlf_daily.rds")

dates <- as.Date(curves$dates)
hours <- curves$hours

Y <- curves$actual
colnames(Y) <- as.character(hours)

has_official <- !is.null(curves$forecast)

if (has_official) {
  Y_official <- curves$forecast
  colnames(Y_official) <- as.character(hours)
}

# align TLF to load dates
tlf_aligned <- tibble(date = dates) |>
  left_join(tlf_daily |> dplyr::select(date, tlf_mean), by = "date")

ell <- tlf_aligned$tlf_mean

ok <- complete.cases(Y)

dates <- dates[ok]
Y <- Y[ok, , drop = FALSE]
ell <- ell[ok]

if (has_official) {
  Y_official <- Y_official[ok, , drop = FALSE]
}

n_days <- nrow(Y)

cat("\nTLF missing share:", mean(is.na(ell)), "\n")

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
# Forecast one day
# ------------------------------------------------------------

forecast_one_day_tlf <- function(i, Y, ell, hours, K, nbasis, norder, tlf_q) {
  
  Y_train <- Y[1:(i - 1), , drop = FALSE]
  ell_train <- ell[1:(i - 1)]
  
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
  
  # training quantile only
  q_hat <- quantile(ell_train, probs = tlf_q, na.rm = TRUE)
  
  # response days: d = 8,...,n_train
  # use TLF_{d-1}
  S_lag1 <- as.integer(ell_train[7:(n_train - 1)] > q_hat)
  S_lag1[is.na(S_lag1)] <- 0
  
  Y_score <- scores[8:n_train, , drop = FALSE]
  
  X_score <- cbind(
    intercept = 1,
    lag1 = scores[7:(n_train - 1), , drop = FALSE],
    lag7 = scores[1:(n_train - 7), , drop = FALSE],
    high_tlf_lag1 = S_lag1
  )
  
  colnames(X_score) <- c(
    "intercept",
    paste0("lag1_score", 1:K),
    paste0("lag7_score", 1:K),
    "high_tlf_lag1"
  )
  
  fit <- lm.fit(X_score, Y_score)
  B_hat <- fit$coefficients
  
  # forecast day i
  S_new <- as.integer(ell_train[n_train] > q_hat)
  if (is.na(S_new)) S_new <- 0
  
  x_new <- c(
    1,
    scores[n_train, ],
    scores[n_train - 6, ],
    S_new
  )
  
  score_hat <- as.vector(x_new %*% B_hat)
  
  mu_hat <- as.vector(eval.fd(hours, pca$meanfd))
  phi_hat <- eval.fd(hours, pca$harmonics)
  
  y_hat <- as.vector(mu_hat + phi_hat %*% score_hat)
  
  list(
    forecast = y_hat,
    score_hat = score_hat,
    high_tlf_lag1 = S_new,
    q_hat = q_hat
  )
}

# ------------------------------------------------------------
# Rolling forecasts
# ------------------------------------------------------------

pred_tlf <- matrix(
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

regime_track <- numeric(length(test_idx))
q_track <- numeric(length(test_idx))

for (j in seq_along(test_idx)) {
  
  i <- test_idx[j]
  
  out <- forecast_one_day_tlf(
    i = i,
    Y = Y,
    ell = ell,
    hours = hours,
    K = K,
    nbasis = NBASIS,
    norder = NORDER,
    tlf_q = TLF_Q
  )
  
  pred_tlf[j, ] <- out$forecast
  score_forecasts[j, ] <- out$score_hat
  regime_track[j] <- out$high_tlf_lag1
  q_track[j] <- out$q_hat
  
  if (j %% 25 == 0) {
    message("Finished ", j, " / ", length(test_idx),
            " | date: ", dates[i])
  }
}

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
    year = year(date),
    month = month(date),
    ym = format(date, "%Y-%m"),
    hour = rep(hours, times = length(test_idx)),
    model = model_name,
    actual = as.vector(t(actual)),
    forecast = as.vector(t(pred))
  ) |>
    mutate(
      error = actual - forecast,
      abs_error = abs(error),
      sq_error = error^2
    )
}

eval_df <- bind_rows(
  make_eval_df(actual_test, pred_naive, "Seasonal naive: Y[d-7]"),
  make_eval_df(actual_test, pred_tlf, "FPCA VAR(1,7) + TLF regime")
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
  "output/tables/fpca_VAR_1_7_TLF_regime_metrics_overall.csv"
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

write_csv(metrics_by_hour, "output/tables/fpca_VAR_1_7_TLF_regime_by_hour.csv")
write_csv(metrics_by_month, "output/tables/fpca_VAR_1_7_TLF_regime_by_month.csv")
write_csv(metrics_by_year, "output/tables/fpca_VAR_1_7_TLF_regime_by_year.csv")

# ------------------------------------------------------------
# Compare to previous FPCA VAR(1,7), if available
# ------------------------------------------------------------

prev_file <- "output/tables/fpca_VAR_1_7_eval_long.csv"

if (file.exists(prev_file)) {
  
  prev_eval <- read_csv(prev_file, show_col_types = FALSE) |>
    filter(model == "FPCA VAR(1,7) scores") |>
    mutate(
      date = as.Date(date),
      hour = as.integer(hour)
    ) |>
    dplyr::select(date, hour, abs_error_prev = abs_error, sq_error_prev = sq_error)
  
  new_eval <- eval_df |>
    filter(model == "FPCA VAR(1,7) + TLF regime") |>
    dplyr::select(date, hour, abs_error_tlf = abs_error, sq_error_tlf = sq_error)
  
  comparison <- prev_eval |>
    inner_join(new_eval, by = c("date", "hour")) |>
    mutate(
      mae_gain_from_tlf = abs_error_prev - abs_error_tlf,
      mse_gain_from_tlf = sq_error_prev - sq_error_tlf
    )
  
  comp_overall <- comparison |>
    summarise(
      mae_gain_from_tlf = mean(mae_gain_from_tlf, na.rm = TRUE),
      mse_gain_from_tlf = mean(mse_gain_from_tlf, na.rm = TRUE)
    )
  
  comp_hour <- comparison |>
    group_by(hour) |>
    summarise(
      mae_gain_from_tlf = mean(mae_gain_from_tlf, na.rm = TRUE),
      mse_gain_from_tlf = mean(mse_gain_from_tlf, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\nGain from adding TLF relative to FPCA VAR(1,7):\n")
  print(comp_overall)
  
  write_csv(comp_overall, "output/tables/fpca_VAR_1_7_TLF_gain_overall.csv")
  write_csv(comp_hour, "output/tables/fpca_VAR_1_7_TLF_gain_by_hour.csv")
}

# ------------------------------------------------------------
# Regime diagnostics
# ------------------------------------------------------------

regime_df <- tibble(
  date = dates[test_idx],
  high_tlf_lag1 = regime_track,
  q75_train = q_track,
  tlf_lag1 = ell[test_idx - 1]
)

write_csv(regime_df, "output/tables/fpca_VAR_1_7_TLF_regime_indicator.csv")

cat("\nHigh-TLF regime counts:\n")
print(table(regime_df$high_tlf_lag1, useNA = "ifany"))

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
  "output/figs/fpca_VAR_1_7_TLF_regime_mae_by_hour.png",
  p_hour,
  width = 8,
  height = 4.8
)

p_month <- metrics_by_month |>
  ggplot(aes(ym, mae, group = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  labs(
    x = "Month",
    y = "MAE",
    title = "Forecast MAE by month",
    linetype = "Model"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(
  "output/figs/fpca_VAR_1_7_TLF_regime_mae_by_month.png",
  p_month,
  width = 10,
  height = 4.8
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

forecast_objects <- list(
  dates = dates[test_idx],
  hours = hours,
  actual = actual_test,
  fpca_VAR_1_7_TLF_regime = pred_tlf,
  seasonal_naive = pred_naive,
  official_forecast = if (has_official) pred_official else NULL,
  score_forecasts = score_forecasts,
  regime = regime_df,
  settings = list(
    K = K,
    NBASIS = NBASIS,
    NORDER = NORDER,
    TEST_START = TEST_START,
    MIN_TRAIN_DAYS = MIN_TRAIN_DAYS,
    TLF_Q = TLF_Q
  ),
  metrics_overall = metric_table,
  metrics_by_hour = metrics_by_hour,
  metrics_by_month = metrics_by_month,
  metrics_by_year = metrics_by_year
)

saveRDS(
  forecast_objects,
  "output/results/fpca_VAR_1_7_TLF_regime_forecast_results.rds"
)

write_csv(
  eval_df,
  "output/tables/fpca_VAR_1_7_TLF_regime_eval_long.csv"
)


