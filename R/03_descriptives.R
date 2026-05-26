# 03_descriptives.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)

curves <- readRDS("data/processed/load_curves_2023.rds")

hours <- curves$hours
Y <- curves$actual      # 365 x 24
Yhat <- curves$forecast # 365 x 24
E <- curves$error       # 365 x 24

# fda wants: rows = time points, columns = curves
Y_t <- t(Y)
Yhat_t <- t(Yhat)
E_t <- t(E)


# B-spline basis on hour domain
basis <- create.bspline.basis(
  rangeval = c(0, 23),
  nbasis = 12,
  norder = 4
)

actual_fd <- Data2fd(argvals = hours, y = Y_t, basisobj = basis)
forecast_fd <- Data2fd(argvals = hours, y = Yhat_t, basisobj = basis)
error_fd <- Data2fd(argvals = hours, y = E_t, basisobj = basis)

# mean functions
actual_mean_fd <- mean.fd(actual_fd)
forecast_mean_fd <- mean.fd(forecast_fd)
error_mean_fd <- mean.fd(error_fd)

# evaluate on grid
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
    title = "Mean functional load curve, 2023"
  )

ggsave("output/figs/fda_mean_actual_forecast_2023.png", p1, width = 7, height = 4)

p2 <- mean_df |>
  ggplot(aes(hour, error)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "Actual - forecast",
    title = "Mean functional forecast error curve, 2023"
  )

ggsave("output/figs/fda_mean_error_2023.png", p2, width = 7, height = 4)

# Functional PCA
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

write.csv(fpca_table, "output/tables/fda_fpca_variance_2023.csv", row.names = FALSE)

# plot actual FPCs
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

ggsave("output/figs/fda_fpca_loadings_2023.png", p3, width = 7, height = 4)

print(fpca_table)

saveRDS(
  list(
    basis = basis,
    actual_fd = actual_fd,
    forecast_fd = forecast_fd,
    error_fd = error_fd,
    actual_pca = actual_pca,
    error_pca = error_pca
  ),
  "data/processed/fda_objects_2023.rds"
)

# plot error FPCs
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

ggsave("output/figs/fda_fpca_error_loadings_2023.png", p4, width = 7, height = 4)





