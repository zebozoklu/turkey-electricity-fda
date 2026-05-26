# 04_forecasting.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)
library(lubridate)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)

curves <- readRDS("data/processed/load_curves.rds")

dates <- curves$dates
hours <- curves$hours

Y_actual <- curves$actual
Y_forecast <- curves$forecast
Y_error <- curves$error

eval_forecast <- function(actual, pred, model_name, test_year) {
  err <- actual - pred
  
  tibble(
    test_year = test_year,
    model = model_name,
    mae = mean(abs(err), na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    bias = mean(err, na.rm = TRUE)
  )
}

center_fd <- function(fdobj, meanfd) {
  fdobj$coefs <- sweep(fdobj$coefs, 1, as.vector(meanfd$coefs), "-")
  fdobj
}

run_one_test_year <- function(test_year, K = 4, nbasis = 12) {
  
  train_idx <- year(dates) < test_year
  test_idx  <- year(dates) == test_year
  
  E_train <- Y_error[train_idx, ]
  E_test  <- Y_error[test_idx, ]
  
  Y_test_actual <- Y_actual[test_idx, ]
  Y_test_forecast <- Y_forecast[test_idx, ]
  test_dates <- dates[test_idx]
  
  basis <- create.bspline.basis(
    rangeval = c(0, 23),
    nbasis = nbasis,
    norder = 4
  )
  
  E_train_fd <- Data2fd(
    argvals = hours,
    y = t(E_train),
    basisobj = basis
  )
  
  E_test_fd <- Data2fd(
    argvals = hours,
    y = t(E_test),
    basisobj = basis
  )
  
  # FPCA fitted only on training data: no test leakage
  error_pca_train <- pca.fd(E_train_fd, nharm = K)
  
  mean_error_fd <- mean.fd(E_train_fd)
  phi_fd <- error_pca_train$harmonics
  
  scores_train <- error_pca_train$scores[, 1:K, drop = FALSE]
  
  # project test errors onto training FPCs
  E_test_centered_fd <- center_fd(E_test_fd, mean_error_fd)
  scores_test <- inprod(E_test_centered_fd, phi_fd)
  
  if (nrow(scores_test) != nrow(E_test)) {
    scores_test <- t(scores_test)
  }
  
  # VAR(1) score dynamics: xi_d = a + A xi_{d-1} + eta_d
  X_lag <- scores_train[-nrow(scores_train), , drop = FALSE]
  Y_now <- scores_train[-1, , drop = FALSE]
  X_reg <- cbind(1, X_lag)
  
  B_hat <- qr.solve(X_reg, Y_now)
  
  pred_scores <- matrix(NA_real_, nrow = nrow(scores_test), ncol = K)
  
  for (i in seq_len(nrow(scores_test))) {
    if (i == 1) {
      lag_score <- scores_train[nrow(scores_train), ]
    } else {
      # rolling one-step forecast: yesterday's realized error is known
      lag_score <- scores_test[i - 1, ]
    }
    
    pred_scores[i, ] <- c(1, lag_score) %*% B_hat
  }
  
  mean_error_vec <- as.vector(eval.fd(hours, mean_error_fd))
  phi_mat <- eval.fd(hours, phi_fd)   # 24 x K
  
  Ehat_fpca_ar <- pred_scores %*% t(phi_mat)
  Ehat_fpca_ar <- sweep(Ehat_fpca_ar, 2, mean_error_vec, "+")
  
  # benchmarks
  Ehat_mean <- matrix(
    rep(mean_error_vec, each = nrow(Y_test_forecast)),
    nrow = nrow(Y_test_forecast),
    ncol = length(hours)
  )
  
  Ehat_yesterday <- matrix(NA_real_, nrow = nrow(Y_test_forecast), ncol = length(hours))
  
  for (i in seq_len(nrow(Ehat_yesterday))) {
    if (i == 1) {
      Ehat_yesterday[i, ] <- E_train[nrow(E_train), ]
    } else {
      Ehat_yesterday[i, ] <- E_test[i - 1, ]
    }
  }
  
  # corrected forecasts
  Yhat_official <- Y_test_forecast
  Yhat_mean_corr <- Y_test_forecast + Ehat_mean
  Yhat_yday_corr <- Y_test_forecast + Ehat_yesterday
  Yhat_fpca_ar_corr <- Y_test_forecast + Ehat_fpca_ar
  
  eval_table <- bind_rows(
    eval_forecast(Y_test_actual, Yhat_official, "official", test_year),
    eval_forecast(Y_test_actual, Yhat_mean_corr, "official_plus_train_mean_error", test_year),
    eval_forecast(Y_test_actual, Yhat_yday_corr, "official_plus_yesterday_error", test_year),
    eval_forecast(Y_test_actual, Yhat_fpca_ar_corr, "official_plus_dynamic_fpca_error", test_year)
  )
  
  hourly_eval <- tibble(
    test_year = test_year,
    hour = rep(hours, 4),
    model = rep(eval_table$model, each = length(hours)),
    mae = c(
      colMeans(abs(Y_test_actual - Yhat_official), na.rm = TRUE),
      colMeans(abs(Y_test_actual - Yhat_mean_corr), na.rm = TRUE),
      colMeans(abs(Y_test_actual - Yhat_yday_corr), na.rm = TRUE),
      colMeans(abs(Y_test_actual - Yhat_fpca_ar_corr), na.rm = TRUE)
    )
  )
  
  list(
    eval_table = eval_table,
    hourly_eval = hourly_eval,
    test_dates = test_dates,
    actual = Y_test_actual,
    official = Yhat_official,
    corrected_mean = Yhat_mean_corr,
    corrected_yesterday = Yhat_yday_corr,
    corrected_fpca_ar = Yhat_fpca_ar_corr
  )
}

test_years <- c(2023, 2024, 2025)

results <- lapply(test_years, run_one_test_year)

eval_all <- bind_rows(lapply(results, `[[`, "eval_table"))
hourly_eval_all <- bind_rows(lapply(results, `[[`, "hourly_eval"))

# improvement relative to official within year
official_mae <- eval_all |>
  filter(model == "official") |>
  dplyr::select(test_year, official_mae = mae, official_rmse = rmse)

eval_all <- eval_all |>
  left_join(official_mae, by = "test_year") |>
  mutate(
    mae_improvement_pct = 100 * (official_mae - mae) / official_mae,
    rmse_improvement_pct = 100 * (official_rmse - rmse) / official_rmse
  )

print(eval_all)

write.csv(eval_all, "output/tables/forecast_eval_all_test_years.csv", row.names = FALSE)
write.csv(hourly_eval_all, "output/tables/forecast_hourly_eval_all_test_years.csv", row.names = FALSE)

p_eval <- hourly_eval_all |>
  ggplot(aes(hour, mae, linetype = model)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ test_year, scales = "free_y") +
  labs(
    x = "Hour",
    y = "MAE",
    title = "Hourly MAE: official vs corrected forecasts"
  )

ggsave("output/figs/forecast_eval_hourly_all_test_years.png", p_eval, width = 10, height = 5)

saveRDS(
  list(
    results = results,
    eval_all = eval_all,
    hourly_eval_all = hourly_eval_all
  ),
  "data/processed/forecast_results_all_test_years.rds"
)

