# R/12_diagnostics_holiday_regime_controls.R

library(dplyr)
library(tidyr)
library(readr)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

RESULT_PATH <- "output/results/fpca_VAR_1_7_holiday_regime_local_correction_forecast_results.rds"
BY_HOUR_PATH <- "output/tables/fpca_VAR_1_7_holiday_regime_local_correction_by_hour.csv"

obj <- readRDS(RESULT_PATH)

actual <- obj$actual
base <- obj$base_holiday_regime
correction <- obj$correction
dates <- as.Date(obj$dates)
hours <- obj$hours

stopifnot(all(dim(actual) == dim(base)))
stopifnot(all(dim(actual) == dim(correction)))
stopifnot(length(dates) == nrow(actual))
stopifnot(length(hours) == ncol(actual))

evaluate_matrix <- function(pred) {
  err <- as.vector(actual - pred)
  
  tibble(
    n = sum(!is.na(err)),
    mae = mean(abs(err), na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    bias = mean(err, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# Shrinkage sensitivity
# ------------------------------------------------------------

lambda_grid <- seq(0, 1.2, by = 0.1)

shrinkage_sensitivity <- bind_rows(lapply(lambda_grid, function(lambda) {
  evaluate_matrix(base + lambda * correction) |>
    mutate(lambda = lambda, .before = 1)
}))

write_csv(
  shrinkage_sensitivity,
  "output/tables/fpca_VAR_1_7_holiday_regime_local_correction_shrinkage_sensitivity.csv"
)

print(shrinkage_sensitivity)

# ------------------------------------------------------------
# Placebo corrections
# ------------------------------------------------------------

set.seed(1)

shift_rows <- function(mat, k) {
  n <- nrow(mat)
  mat[c((n - k + 1):n, 1:(n - k)), , drop = FALSE]
}

hour_mean_correction <- matrix(
  rep(colMeans(correction, na.rm = TRUE), each = nrow(correction)),
  nrow = nrow(correction)
)

placebo_metrics <- bind_rows(
  evaluate_matrix(base) |> mutate(control = "base_no_correction"),
  evaluate_matrix(base + 0.6 * correction) |> mutate(control = "true_lambda_0.6"),
  evaluate_matrix(base + 1.0 * correction) |> mutate(control = "true_lambda_1.0"),
  evaluate_matrix(base + 0.6 * shift_rows(correction, 1)) |>
    mutate(control = "shifted_1_day"),
  evaluate_matrix(base + 0.6 * shift_rows(correction, 7)) |>
    mutate(control = "shifted_7_days"),
  evaluate_matrix(base + 0.6 * correction[sample(seq_len(nrow(correction))), , drop = FALSE]) |>
    mutate(control = "permuted_dates"),
  evaluate_matrix(base + 0.6 * hour_mean_correction) |>
    mutate(control = "hour_mean_correction")
) |>
  relocate(control)

write_csv(
  placebo_metrics,
  "output/tables/fpca_VAR_1_7_holiday_regime_local_correction_placebo_controls.csv"
)

print(placebo_metrics)

# ------------------------------------------------------------
# Official-winning hours check
# ------------------------------------------------------------

by_hour <- read_csv(BY_HOUR_PATH, show_col_types = FALSE)

model_name <- "FPCA VAR(1,7) + holiday/load regime + nested adaptive correction"

official_wins <- by_hour |>
  filter(model %in% c(model_name, "Official forecast")) |>
  dplyr::select(hour, model, mae, rmse) |>
  pivot_wider(
    names_from = model,
    values_from = c(mae, rmse),
    names_repair = "unique"
  )

names(official_wins) <- make.names(names(official_wins))

model_mae_col <- make.names(paste0("mae_", model_name))
official_mae_col <- make.names("mae_Official forecast")
model_rmse_col <- make.names(paste0("rmse_", model_name))
official_rmse_col <- make.names("rmse_Official forecast")

official_wins <- official_wins |>
  mutate(
    official_mae_gain = .data[[model_mae_col]] - .data[[official_mae_col]],
    official_rmse_gain = .data[[model_rmse_col]] - .data[[official_rmse_col]]
  ) |>
  filter(official_mae_gain > 0 | official_rmse_gain > 0)

write_csv(
  official_wins,
  "output/tables/fpca_VAR_1_7_holiday_regime_local_correction_official_winning_hours.csv"
)

print(official_wins)
