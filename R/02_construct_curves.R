# 02_construct_curves.R

library(dplyr)
library(tidyr)

load_hourly <- readRDS("data/processed/load_hourly.rds")

make_curve_matrix <- function(data, value_col) {
  data |>
    dplyr::select(date, hour_id, {{ value_col }}) |>
    pivot_wider(names_from = hour_id, values_from = {{ value_col }}) |>
    arrange(date)
}

actual_mat <- make_curve_matrix(load_hourly, load_actual_mwh)
forecast_mat <- make_curve_matrix(load_hourly, load_forecast_mwh)
error_mat <- make_curve_matrix(load_hourly, forecast_error_mwh)

dates <- actual_mat$date

Y_actual <- as.matrix(actual_mat |> dplyr::select(-date))
Y_forecast <- as.matrix(forecast_mat |> dplyr::select(-date))
Y_error <- as.matrix(error_mat |> dplyr::select(-date))

colnames(Y_actual) <- 0:23
colnames(Y_forecast) <- 0:23
colnames(Y_error) <- 0:23

stopifnot(ncol(Y_actual) == 24)
stopifnot(all(dim(Y_actual) == dim(Y_forecast)))
stopifnot(all(dim(Y_actual) == dim(Y_error)))
stopifnot(nrow(Y_actual) == length(unique(load_hourly$date)))

curves <- list(
  dates = dates,
  hours = 0:23,
  actual = Y_actual,
  forecast = Y_forecast,
  error = Y_error
)

saveRDS(curves, "data/processed/load_curves.rds")


