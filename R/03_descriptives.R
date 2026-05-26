# 03_descriptives.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)
library(lubridate)

dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

load_hourly <- readRDS("data/processed/load_hourly.rds")
curves <- readRDS("data/processed/load_curves.rds")

hours <- curves$hours
dates <- curves$dates

Y <- curves$actual
Yhat <- curves$forecast
E <- curves$error

# fda wants: rows = time points, columns = curves
Y_t <- t(Y)
Yhat_t <- t(Yhat)
E_t <- t(E)

# ------------------------------------------------------------
# Simple multi-year descriptive plots
# ------------------------------------------------------------

mean_by_year <- load_hourly |>
  group_by(year, hour_id) |>
  summarise(
    actual = mean(load_actual_mwh, na.rm = TRUE),
    forecast = mean(load_forecast_mwh, na.rm = TRUE),
    error = mean(forecast_error_mwh, na.rm = TRUE),
    mae = mean(abs_error_mwh, na.rm = TRUE),
    .groups = "drop"
  )

p_year_error <- mean_by_year |>
  ggplot(aes(hour_id, error, linetype = factor(year))) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "Actual - forecast",
    linetype = "Year",
    title = "Mean load forecast error curve by year"
  )

ggsave("output/figs/mean_error_by_year.png", p_year_error, width = 8, height = 4.5)

p_year_mae <- mean_by_year |>
  ggplot(aes(hour_id, mae, linetype = factor(year))) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "MAE",
    linetype = "Year",
    title = "Hourly MAE of official load forecast by year"
  )

ggsave("output/figs/mae_by_hour_year.png", p_year_mae, width = 8, height = 4.5)

# ------------------------------------------------------------
# Functional objects
# ------------------------------------------------------------

basis <- create.bspline.basis(
  rangeval = c(0, 23),
  nbasis = 12,
  norder = 4
)

actual_fd <- Data2fd(argvals = hours, y = Y_t, basisobj = basis)
forecast_fd <- Data2fd(argvals = hours, y = Yhat_t, basisobj = basis)
error_fd <- Data2fd(argvals = hours, y = E_t, basisobj = basis)

# ------------------------------------------------------------
# Mean functional curves
# ------------------------------------------------------------

actual_mean_fd <- mean.fd(actual_fd)
forecast_mean_fd <- mean.fd(forecast_fd)
error_mean_fd <- mean.fd(error_fd)

grid <- seq(0, 23, length.out = 100)

mean_df <- tibble(
  hour = grid,
  actual = as.vector(eval.fd(grid, actual_mean_fd)),
  forecast = as.vector(eval.fd(grid, forecast_mean_fd)),
  error = as.vector(eval.fd(grid, error_mean_fd))
)

p1 <- mean_df |>
  dplyr::select(hour, actual, forecast) |>
  pivot_longer(-hour, names_to = "series", values_to = "mwh") |>
  ggplot(aes(hour, mwh, linetype = series)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "MWh",
    title = "Mean functional load curve, 2019–2025"
  )

ggsave("output/figs/fda_mean_actual_forecast_all.png", p1, width = 7, height = 4)

p2 <- mean_df |>
  ggplot(aes(hour, error)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "Actual - forecast",
    title = "Mean functional forecast error curve, 2019–2025"
  )

ggsave("output/figs/fda_mean_error_all.png", p2, width = 7, height = 4)

# ------------------------------------------------------------
# Functional PCA
# ------------------------------------------------------------

actual_pca <- pca.fd(actual_fd, nharm = 4)
error_pca <- pca.fd(error_fd, nharm = 4)

actual_var <- actual_pca$varprop
error_var <- error_pca$varprop

fpca_table <- tibble(
  component = 1:4,
  actual_variance_explained = actual_var,
  actual_cumulative = cumsum(actual_var),
  error_variance_explained = error_var,
  error_cumulative = cumsum(error_var)
)

write.csv(fpca_table, "output/tables/fda_fpca_variance_all.csv", row.names = FALSE)

# actual FPCs
harmonics_eval <- eval.fd(grid, actual_pca$harmonics)

fpc_df <- tibble(
  hour = grid,
  FPC1 = harmonics_eval[, 1],
  FPC2 = harmonics_eval[, 2],
  FPC3 = harmonics_eval[, 3],
  FPC4 = harmonics_eval[, 4]
) |>
  pivot_longer(-hour, names_to = "component", values_to = "loading")

p3 <- fpc_df |>
  ggplot(aes(hour, loading, linetype = component)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "Loading",
    title = "Functional principal components of load curves"
  )

ggsave("output/figs/fda_fpca_loadings_all.png", p3, width = 7, height = 4)

# error FPCs
error_harmonics_eval <- eval.fd(grid, error_pca$harmonics)

error_fpc_df <- tibble(
  hour = grid,
  FPC1 = error_harmonics_eval[, 1],
  FPC2 = error_harmonics_eval[, 2],
  FPC3 = error_harmonics_eval[, 3],
  FPC4 = error_harmonics_eval[, 4]
) |>
  pivot_longer(-hour, names_to = "component", values_to = "loading")

p4 <- error_fpc_df |>
  ggplot(aes(hour, loading, linetype = component)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "Loading",
    title = "Functional principal components of forecast error curves"
  )

ggsave("output/figs/fda_fpca_error_loadings_all.png", p4, width = 7, height = 4)

# ------------------------------------------------------------
# Error FPCA scores over time: useful for dynamic analysis
# ------------------------------------------------------------

error_scores <- as_tibble(error_pca$scores[, 1:4]) |>
  setNames(paste0("score", 1:4)) |>
  mutate(date = dates, year = year(date)) |>
  dplyr::select(date, year, everything())

write.csv(error_scores, "output/tables/error_fpca_scores_all.csv", row.names = FALSE)

score_long <- error_scores |>
  pivot_longer(starts_with("score"), names_to = "component", values_to = "score")

p5 <- score_long |>
  ggplot(aes(date, score)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ component, scales = "free_y", ncol = 1) +
  labs(
    x = "Date",
    y = "Score",
    title = "Forecast error FPCA scores over time"
  )

ggsave("output/figs/error_fpca_scores_time_all.png", p5, width = 8, height = 7)

print(fpca_table)

saveRDS(
  list(
    basis = basis,
    actual_fd = actual_fd,
    forecast_fd = forecast_fd,
    error_fd = error_fd,
    actual_pca = actual_pca,
    error_pca = error_pca,
    error_scores = error_scores
  ),
  "data/processed/fda_objects.rds"
)

