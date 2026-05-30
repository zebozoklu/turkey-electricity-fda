# 04_forecast_simple_fpca_weekly.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)
library(lubridate)

dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

K <- 2              # number of FPCs
NBASIS <- 12        # B-spline basis size
NORDER <- 4         # cubic B-splines
TEST_START <- as.Date("2023-01-01")
MIN_TRAIN_DAYS <- 365
DEBUG_N <- Inf      # set e.g. 30 for fast trial run

# ------------------------------------------------------------
# Load daily curves
# ------------------------------------------------------------

curves <- readRDS("data/processed/load_curves.rds")

dates <- as.Date(curves$dates)
hours <- curves$hours
Y <- curves$actual

colnames(Y) <- as.character(hours)

# optional official forecast benchmark
has_official <- !is.null(curves$forecast)
if (has_official) {
  Y_official <- curves$forecast
  colnames(Y_official) <- as.character(hours)
}

# keep complete actual curves
ok <- complete.cases(Y)
dates <- dates[ok]
Y <- Y[ok, , drop = FALSE]

if (has_official) {
  Y_official <- Y_official[ok, , drop = FALSE]
}

n_days <- nrow(Y)

# test period
test_idx <- which(dates >= TEST_START)
test_idx <- test_idx[test_idx > max(7, MIN_TRAIN_DAYS)]

# fallback if TEST_START not available
if (length(test_idx) == 0) {
  first_test <- floor(0.8 * n_days) + 1
  test_idx <- seq(first_test, n_days)
  test_idx <- test_idx[test_idx > max(7, MIN_TRAIN_DAYS)]
}

if (is.finite(DEBUG_N)) {
  test_idx <- head(test_idx, DEBUG_N)
}

message("Number of forecast days: ", length(test_idx))
message("First forecast date: ", min(dates[test_idx]))
message("Last forecast date: ", max(dates[test_idx]))

# ------------------------------------------------------------
# One-step-ahead simple FPCA weekly forecast
# ------------------------------------------------------------

forecast_one_day <- function(i, Y, hours, K, nbasis, norder) {
  
  # training sample uses only information before target day i
  Y_train <- Y[1:(i - 1), , drop = FALSE]
  n_train <- nrow(Y_train)
  
  # functional object
  basis <- create.bspline.basis(
    rangeval = range(hours),
    nbasis = nbasis,
    norder = norder
  )
  
  fd_train <- Data2fd(
    argvals = hours,
    y = t(Y_train),          # rows = hours, cols = days
    basisobj = basis
  )
  
  # FPCA estimated only on training curves
  pca <- pca.fd(fd_train, nharm = K)
  
  scores <- pca$scores[, 1:K, drop = FALSE]
  
  # target day i uses same weekday last week: i - 7
  # in training coordinates, that is row n_train - 6
  lag_score <- scores[n_train - 6, ]
  
  score_hat <- numeric(K)
  
  for (k in 1:K) {
    
    y_k <- scores[8:n_train, k]
    x_k <- scores[1:(n_train - 7), k]
    
    X_k <- cbind(1, x_k)
    
    fit_k <- lm.fit(X_k, y_k)
    
    score_hat[k] <- c(1, lag_score[k]) %*% fit_k$coefficients
  }
  
  # reconstruct curve on hourly grid
  mu_hat <- as.vector(eval.fd(hours, pca$meanfd))
  phi_hat <- eval.fd(hours, pca$harmonics)
  
  y_hat <- as.vector(mu_hat + phi_hat %*% score_hat)
  
  list(
    forecast = y_hat,
    score_hat = score_hat,
    varprop = pca$varprop
  )
}

# ------------------------------------------------------------
# Rolling / expanding forecasts
# ------------------------------------------------------------

pred_fpca <- matrix(
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

for (j in seq_along(test_idx)) {
  
  i <- test_idx[j]
  
  out <- forecast_one_day(
    i = i,
    Y = Y,
    hours = hours,
    K = K,
    nbasis = NBASIS,
    norder = NORDER
  )
  
  pred_fpca[j, ] <- out$forecast
  score_forecasts[j, ] <- out$score_hat
  varprops[j, ] <- out$varprop[1:K]
  
  if (j %% 25 == 0) {
    message("Finished ", j, " / ", length(test_idx),
            " forecasts. Current date: ", dates[i])
  }
}

# ------------------------------------------------------------
# Benchmarks
# ------------------------------------------------------------

actual_test <- Y[test_idx, , drop = FALSE]

# seasonal naive: copy same weekday last week
pred_naive <- Y[test_idx - 7, , drop = FALSE]

if (has_official) {
  pred_official <- Y_official[test_idx, , drop = FALSE]
}

# ------------------------------------------------------------
# Evaluation data
# ------------------------------------------------------------

make_eval_df <- function(actual, pred, model_name) {
  tibble(
    date = rep(dates[test_idx], each = length(hours)),
    year = year(date),
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
  make_eval_df(actual_test, pred_fpca, "Simple FPCA weekly score AR")
)

if (has_official) {
  eval_df <- bind_rows(
    eval_df,
    make_eval_df(actual_test, pred_official, "Official forecast")
  )
}

metric_table <- eval_df |>
  group_by(model) |>
  summarise(
    n = sum(!is.na(error)),
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

# improvement relative to seasonal naive
base_mae <- metric_table$mae[metric_table$model == "Seasonal naive: Y[d-7]"]
base_rmse <- metric_table$rmse[metric_table$model == "Seasonal naive: Y[d-7]"]

metric_table <- metric_table |>
  mutate(
    mae_improvement_vs_naive_pct = 100 * (base_mae - mae) / base_mae,
    rmse_improvement_vs_naive_pct = 100 * (base_rmse - rmse) / base_rmse
  )

print(metric_table)

write.csv(
  metric_table,
  "output/tables/simple_fpca_weekly_metrics_overall.csv",
  row.names = FALSE
)

metrics_by_year <- eval_df |>
  group_by(year, model) |>
  summarise(
    n = sum(!is.na(error)),
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  metrics_by_year,
  "output/tables/simple_fpca_weekly_metrics_by_year.csv",
  row.names = FALSE
)

metrics_by_hour <- eval_df |>
  group_by(hour, model) |>
  summarise(
    n = sum(!is.na(error)),
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  metrics_by_hour,
  "output/tables/simple_fpca_weekly_metrics_by_hour.csv",
  row.names = FALSE
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
    title = "Forecast accuracy by hour",
    linetype = "Model"
  )

ggsave(
  "output/figs/simple_fpca_weekly_mae_by_hour.png",
  p_hour,
  width = 8,
  height = 4.5
)

# example curve: first test day
example_day <- 1

example_df <- tibble(
  hour = hours,
  Actual = actual_test[example_day, ],
  `Seasonal naive` = pred_naive[example_day, ],
  `Simple FPCA` = pred_fpca[example_day, ]
)

if (has_official) {
  example_df$`Official forecast` <- pred_official[example_day, ]
}

example_long <- example_df |>
  pivot_longer(-hour, names_to = "series", values_to = "load_mwh")

p_example <- example_long |>
  ggplot(aes(hour, load_mwh, linetype = series)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "Load",
    title = paste("Example forecast curve:", dates[test_idx[example_day]]),
    linetype = "Series"
  )

ggsave(
  "output/figs/simple_fpca_weekly_example_curve.png",
  p_example,
  width = 8,
  height = 4.5
)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

forecast_objects <- list(
  dates = dates[test_idx],
  hours = hours,
  actual = actual_test,
  fpca_forecast = pred_fpca,
  seasonal_naive = pred_naive,
  official_forecast = if (has_official) pred_official else NULL,
  score_forecasts = score_forecasts,
  varprops = varprops,
  settings = list(
    K = K,
    NBASIS = NBASIS,
    NORDER = NORDER,
    TEST_START = TEST_START,
    MIN_TRAIN_DAYS = MIN_TRAIN_DAYS
  ),
  metrics_overall = metric_table,
  metrics_by_year = metrics_by_year,
  metrics_by_hour = metrics_by_hour
)

saveRDS(
  forecast_objects,
  "output/results/simple_fpca_weekly_forecast_results.rds"
)

write.csv(
  eval_df,
  "output/tables/simple_fpca_weekly_eval_long.csv",
  row.names = FALSE
)
