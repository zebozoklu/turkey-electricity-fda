# 15_eigenfunction_stability.R
#
# Checks whether FPCA eigenfunctions are stable over time.
# Two tests:
#   (1) Variance proportions (FPC1, FPC2) from rolling-origin FPCA over the
#       test period — already saved in the forecast RDS.
#   (2) Eigenfunction shape comparison: fit FPCA separately on two sub-periods
#       (2019-2021 vs 2022-2024) and compare phi_1, phi_2 visually and by
#       inner product. Inner product near 1 = same direction; near 0 = rotated.
#
# Conclusion informs whether dynamic FPCA is worth pursuing.

library(fda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

dir.create("output/figs",   recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

NBASIS <- 12
NORDER <- 4
K      <- 2

curves    <- readRDS("data/processed/load_curves.rds")
dates_all <- as.Date(curves$dates)
hours     <- curves$hours
Y_all     <- curves$actual
ok        <- complete.cases(Y_all)
dates_all <- dates_all[ok]
Y_all     <- Y_all[ok, ]

basis_obj <- create.bspline.basis(rangeval = range(hours),
                                   nbasis   = NBASIS,
                                   norder   = NORDER)

# ============================================================
# 1. Variance proportions over the rolling test period
# ============================================================

res    <- readRDS("output/results/fpca_VAR_1_7_holiday_regime_derivative_forecast_results.rds")
vp     <- as.data.frame(res$varprops)
vp$date <- res$dates

vp_long <- vp |>
  pivot_longer(c(FPC1, FPC2), names_to = "component", values_to = "var_prop")

p_vp <- ggplot(vp_long, aes(x = date, y = var_prop, colour = component)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~component, ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(
    title    = "Rolling-origin FPCA variance proportions over test period",
    subtitle = "Each point = variance explained by FPC on training window up to that date",
    x = NULL, y = "Variance proportion"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figs/eigenfunction_varprops_rolling.png", p_vp,
       width = 8, height = 5, dpi = 150)
cat("Saved: output/figs/eigenfunction_varprops_rolling.png\n")

cat("\nVariance proportion summary:\n")
vp |>
  summarise(
    FPC1_min = min(FPC1), FPC1_max = max(FPC1), FPC1_sd = sd(FPC1),
    FPC2_min = min(FPC2), FPC2_max = max(FPC2), FPC2_sd = sd(FPC2)
  ) |>
  print()

# ============================================================
# 2. Eigenfunction shape: three sub-periods
# ============================================================
# Split the full 2019-2024 data into three windows and compare eigenfunctions.

windows <- list(
  "2019-2021" = c(as.Date("2019-01-01"), as.Date("2021-12-31")),
  "2020-2022" = c(as.Date("2020-01-01"), as.Date("2022-12-31")),
  "2022-2024" = c(as.Date("2022-01-01"), as.Date("2024-12-31"))
)

fit_fpca_window <- function(start, end) {
  mask  <- dates_all >= start & dates_all <= end
  Y_sub <- Y_all[mask, ]
  fd    <- Data2fd(argvals = hours, y = t(Y_sub), basisobj = basis_obj)
  pca   <- pca.fd(fd, nharm = K)
  list(
    harmonics  = eval.fd(hours, pca$harmonics),   # 24 x K matrix
    varprop    = pca$varprop[1:K],
    n_days     = sum(mask)
  )
}

cat("\nFitting FPCA on sub-periods...\n")
fits <- lapply(windows, function(w) fit_fpca_window(w[1], w[2]))

# Sign-align each window's harmonics to the first window
ref_harmonics <- fits[[1]]$harmonics
for (nm in names(fits)[-1]) {
  for (k in 1:K) {
    ip <- sum(fits[[nm]]$harmonics[, k] * ref_harmonics[, k])
    if (ip < 0) fits[[nm]]$harmonics[, k] <- -fits[[nm]]$harmonics[, k]
  }
}

# Build tidy data frame for plotting
eigen_df <- bind_rows(lapply(names(fits), function(nm) {
  h <- fits[[nm]]$harmonics
  bind_rows(lapply(1:K, function(k) {
    tibble(
      window    = nm,
      component = paste0("FPC ", k),
      hour      = hours,
      value     = h[, k]
    )
  }))
}))

p_eigen <- ggplot(eigen_df, aes(x = hour, y = value,
                                 colour = window, linetype = window)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~component, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  labs(
    title    = "FPCA eigenfunctions across sub-periods",
    subtitle = "Overlapping lines = stable eigenfunctions; divergence = structural change",
    x = "Hour", y = expression(phi[k](h)),
    colour = "Window", linetype = "Window"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/figs/eigenfunction_shape_comparison.png", p_eigen,
       width = 8, height = 6, dpi = 150)
cat("Saved: output/figs/eigenfunction_shape_comparison.png\n")

# ============================================================
# 3. Inner products between sub-period eigenfunctions
# ============================================================
# Inner product of unit-norm eigenfunctions = cosine similarity.
# Value near 1 = same direction; near 0 = orthogonal (fully rotated).

window_names <- names(fits)
pairs <- combn(window_names, 2, simplify = FALSE)

ip_table <- bind_rows(lapply(pairs, function(p) {
  a <- fits[[p[1]]]$harmonics
  b <- fits[[p[2]]]$harmonics
  bind_rows(lapply(1:K, function(k) {
    # Normalize (pca.fd harmonics are already unit-norm in L2, but
    # eval.fd gives discrete values; normalize for safety)
    ak <- a[, k] / sqrt(sum(a[, k]^2))
    bk <- b[, k] / sqrt(sum(b[, k]^2))
    tibble(
      pair      = paste(p[1], "vs", p[2]),
      component = paste0("FPC ", k),
      inner_product = abs(sum(ak * bk))  # abs: sign flip is arbitrary
    )
  }))
}))

cat("\n=== Inner products between sub-period eigenfunctions ===\n")
cat("(1 = identical direction, 0 = orthogonal)\n\n")
print(ip_table, n = 20)
write_csv(ip_table, "output/tables/eigenfunction_inner_products.csv")

# ============================================================
# 3b. Full K×K cross-component inner product matrix
# ============================================================
# Verifies FPC1 and FPC2 have not swapped or rotated between windows.
# Diagonal near 1, off-diagonal near 0 = components are consistent.

cross_ip <- bind_rows(lapply(pairs, function(p) {
  a <- fits[[p[1]]]$harmonics
  b <- fits[[p[2]]]$harmonics
  bind_rows(lapply(1:K, function(j) {
    aj <- a[, j] / sqrt(sum(a[, j]^2))
    bind_rows(lapply(1:K, function(l) {
      bl <- b[, l] / sqrt(sum(b[, l]^2))
      tibble(
        pair          = paste(p[1], "vs", p[2]),
        from          = paste0("FPC ", j),
        to            = paste0("FPC ", l),
        inner_product = abs(sum(aj * bl))
      )
    }))
  }))
}))

cat("\n=== Full K×K cross-component inner product matrices ===\n")
cat("(diagonal ≈ 1 = same component; off-diagonal ≈ 0 = no rotation/swap)\n\n")
cross_ip |>
  pivot_wider(names_from = to, values_from = inner_product) |>
  print(n = 30)

write_csv(cross_ip, "output/tables/eigenfunction_cross_inner_products.csv")
cat("Saved: output/tables/eigenfunction_cross_inner_products.csv\n")

# ============================================================
# 4. Variance proportions per window
# ============================================================

varprop_table <- bind_rows(lapply(names(fits), function(nm) {
  tibble(
    window  = nm,
    n_days  = fits[[nm]]$n_days,
    FPC1_varprop = fits[[nm]]$varprop[1],
    FPC2_varprop = fits[[nm]]$varprop[2]
  )
}))

cat("\n=== Variance proportions by sub-period ===\n")
print(varprop_table)
write_csv(varprop_table, "output/tables/eigenfunction_varprop_by_window.csv")

cat("\nDone.\n")
