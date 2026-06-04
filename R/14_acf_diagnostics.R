# 14_acf_diagnostics.R
#
# ACF diagnostics on forecast errors of the main model.
#
# Two views:
#   (1) ACF of daily MAE series — checks for day-level serial dependence
#       in the overall error magnitude.
#   (2) ACF of hourly error series for selected hours — checks whether
#       the model leaves autocorrelated residuals at specific hours.
#
# Both plots are saved to output/figs/.

library(dplyr)
library(ggplot2)
library(readr)
library(lubridate)

dir.create("output/figs",   recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

MAIN_MODEL <- "FPCA VAR(1,7) + holiday/load/derivative regime"
BENCHMARK  <- "Official forecast"

eval_df <- read_csv(
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  show_col_types = FALSE,
  col_types = cols(date = col_character())
) |>
  mutate(date = as.Date(as.character(.data$date)),
         hour = as.integer(.data$hour))

# ============================================================
# 1. ACF of daily MAE series
# ============================================================

daily_mae <- eval_df |>
  filter(model %in% c(MAIN_MODEL, BENCHMARK)) |>
  group_by(date, model) |>
  summarise(daily_mae = mean(abs_error, na.rm = TRUE), .groups = "drop")

acf_from_series <- function(x, lag.max = 30, label) {
  ac  <- acf(x, lag.max = lag.max, plot = FALSE)
  ci  <- qnorm(0.975) / sqrt(length(x))
  tibble(
    lag      = as.integer(ac$lag[-1]),
    acf      = as.numeric(ac$acf[-1]),
    upper_ci =  ci,
    lower_ci = -ci,
    series   = label
  )
}

main_daily <- daily_mae |> filter(model == MAIN_MODEL) |> arrange(date) |> pull(daily_mae)
off_daily  <- daily_mae |> filter(model == BENCHMARK)  |> arrange(date) |> pull(daily_mae)

acf_daily <- bind_rows(
  acf_from_series(main_daily, label = MAIN_MODEL),
  acf_from_series(off_daily,  label = BENCHMARK)
)

p_daily_acf <- ggplot(acf_daily, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = upper_ci), linetype = "dashed", colour = "steelblue", linewidth = 0.4) +
  geom_hline(aes(yintercept = lower_ci), linetype = "dashed", colour = "steelblue", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~series, ncol = 1) +
  scale_x_continuous(breaks = seq(0, 30, 7)) +
  labs(
    title = "ACF of daily MAE series",
    x     = "Lag (days)",
    y     = "ACF"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 8))

ggsave("output/figs/acf_daily_mae.png", p_daily_acf,
       width = 7, height = 5, dpi = 150)
cat("Saved: output/figs/acf_daily_mae.png\n")

# ============================================================
# 2. ACF of hourly error series for all 24 hours
# ============================================================

main_hourly <- eval_df |>
  filter(model == MAIN_MODEL) |>
  arrange(hour, date)

acf_by_hour <- bind_rows(lapply(0:23, function(h) {
  series <- main_hourly |> filter(hour == h) |> arrange(date) |> pull(error)
  if (length(series) < 10) return(NULL)
  ac  <- acf(series, lag.max = 14, plot = FALSE)
  ci  <- qnorm(0.975) / sqrt(length(series))
  tibble(
    hour     = h,
    lag      = as.integer(ac$lag[-1]),
    acf      = as.numeric(ac$acf[-1]),
    upper_ci =  ci,
    lower_ci = -ci
  )
}))

p_hourly_acf <- ggplot(acf_by_hour, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = upper_ci), linetype = "dashed", colour = "steelblue", linewidth = 0.3) +
  geom_hline(aes(yintercept = lower_ci), linetype = "dashed", colour = "steelblue", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.5) +
  geom_point(size = 0.8) +
  facet_wrap(~hour, ncol = 6, labeller = label_both) +
  scale_x_continuous(breaks = c(1, 7, 14)) +
  labs(
    title = paste("ACF of hourly forecast errors —", MAIN_MODEL),
    x     = "Lag (days)",
    y     = "ACF"
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7))

ggsave("output/figs/acf_hourly_errors.png", p_hourly_acf,
       width = 10, height = 8, dpi = 150)
cat("Saved: output/figs/acf_hourly_errors.png\n")

# ============================================================
# 3. Ljung-Box test for each hour
# ============================================================
# Tests H0: no autocorrelation up to lag 14.
# Rejection means remaining serial structure the VAR hasn't captured.

lb_by_hour <- bind_rows(lapply(0:23, function(h) {
  series <- main_hourly |> filter(hour == h) |> arrange(date) |> pull(error)
  if (length(series) < 20) return(NULL)
  lb <- Box.test(series, lag = 14, type = "Ljung-Box")
  tibble(
    hour      = h,
    LB_stat   = round(lb$statistic, 3),
    p_value   = round(lb$p.value, 4),
    reject_5pct = lb$p.value < 0.05
  )
}))

cat("\n=== Ljung-Box test (lag=14) by hour — main model ===\n")
print(lb_by_hour, n = 24)
write_csv(lb_by_hour, "output/tables/acf_ljung_box_by_hour.csv")

cat("\nHours with significant autocorrelation (5%):",
    sum(lb_by_hour$reject_5pct), "of 24\n")

# ============================================================
# 3b. ACF of TLF (official) hourly errors — 24-hour grid
# ============================================================
# If TLF errors show significant ACF at lags 1 and 7 where our model does not,
# that is direct evidence of curve-shape structure that FDA exploits but
# point-by-point methods leave on the table.

off_hourly <- eval_df |>
  filter(model == BENCHMARK) |>
  arrange(hour, date)

acf_by_hour_off <- bind_rows(lapply(0:23, function(h) {
  series <- off_hourly |> filter(hour == h) |> arrange(date) |> pull(error)
  if (length(series) < 10) return(NULL)
  ac  <- acf(series, lag.max = 14, plot = FALSE)
  ci  <- qnorm(0.975) / sqrt(length(series))
  tibble(
    hour     = h,
    lag      = as.integer(ac$lag[-1]),
    acf      = as.numeric(ac$acf[-1]),
    upper_ci =  ci,
    lower_ci = -ci
  )
}))

p_hourly_acf_off <- ggplot(acf_by_hour_off, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = upper_ci), linetype = "dashed", colour = "steelblue", linewidth = 0.3) +
  geom_hline(aes(yintercept = lower_ci), linetype = "dashed", colour = "steelblue", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.5) +
  geom_point(size = 0.8) +
  facet_wrap(~hour, ncol = 6, labeller = label_both) +
  scale_x_continuous(breaks = c(1, 7, 14)) +
  labs(
    title = paste("ACF of hourly forecast errors —", BENCHMARK),
    x     = "Lag (days)",
    y     = "ACF"
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7))

ggsave("output/figs/acf_hourly_errors_official.png", p_hourly_acf_off,
       width = 10, height = 8, dpi = 150)
cat("Saved: output/figs/acf_hourly_errors_official.png\n")

# ============================================================
# 3c. Ljung-Box for TLF by hour + side-by-side comparison
# ============================================================

lb_by_hour_off <- bind_rows(lapply(0:23, function(h) {
  series <- off_hourly |> filter(hour == h) |> arrange(date) |> pull(error)
  if (length(series) < 20) return(NULL)
  lb <- Box.test(series, lag = 14, type = "Ljung-Box")
  tibble(
    hour        = h,
    LB_stat     = round(lb$statistic, 3),
    p_value     = round(lb$p.value, 4),
    reject_5pct = lb$p.value < 0.05
  )
}))

cat("\n=== Ljung-Box test (lag=14) by hour —", BENCHMARK, "===\n")
print(lb_by_hour_off, n = 24)
cat("\nHours with significant autocorrelation (5%):",
    sum(lb_by_hour_off$reject_5pct), "of 24\n")

lb_compare <- bind_rows(
  lb_by_hour     |> mutate(model = "FDA model"),
  lb_by_hour_off |> mutate(model = "Official forecast")
)
write_csv(lb_compare, "output/tables/acf_ljung_box_comparison.csv")
cat("Saved: output/tables/acf_ljung_box_comparison.csv\n")

# All p-values are ~0 so plot LB statistic magnitude instead.
# Larger LB stat = more unexploited serial structure remaining in errors.
p_lb_compare <- ggplot(lb_compare, aes(x = hour, y = LB_stat, colour = model)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("FDA model" = "#2166ac", "Official forecast" = "#d6604d"),
                      name   = NULL) +
  scale_x_continuous(breaks = 0:23) +
  labs(
    title    = "Ljung-Box statistic (lag = 14) by hour: FDA model vs official forecast",
    subtitle = "Higher = more unexploited serial structure remaining in forecast errors",
    x        = "Hour of day",
    y        = "Ljung-Box Q statistic"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave("output/figs/acf_lb_pvalue_comparison.png", p_lb_compare,
       width = 9, height = 4.5, dpi = 150)
cat("Saved: output/figs/acf_lb_pvalue_comparison.png\n")

# ============================================================
# 3d. ACF heatmap: |ACF| by lag x hour for both models
# ============================================================
# Visual summary: warmer colour = stronger remaining autocorrelation.
# More structure (darker) in the TLF panel supports the FDA advantage claim.

acf_heat <- bind_rows(
  acf_by_hour     |> mutate(model = "FDA model"),
  acf_by_hour_off |> mutate(model = "Official forecast")
)

p_acf_heat <- ggplot(acf_heat, aes(x = factor(lag), y = factor(hour), fill = abs(acf))) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~model, ncol = 2) +
  scale_fill_distiller(palette = "OrRd", direction = 1,
                       name   = "|ACF|",
                       limits = c(0, NA)) +
  scale_y_discrete(breaks = as.character(seq(0, 23, 4))) +
  labs(
    title    = "ACF heatmap: absolute autocorrelation by lag and hour of day",
    subtitle = "Darker = stronger remaining serial structure in forecast errors",
    x        = "Lag (days)",
    y        = "Hour of day"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text   = element_text(size = 9),
        panel.grid   = element_blank(),
        axis.text.x  = element_text(size = 7))

ggsave("output/figs/acf_heatmap_comparison.png", p_acf_heat,
       width = 10, height = 5, dpi = 150)
cat("Saved: output/figs/acf_heatmap_comparison.png\n")

cat("\nDone.\n")

# ============================================================
# 4. ACF of in-sample VAR(1,7) score residuals
# ============================================================
# Fit FPCA + VAR(1,7) on the full pre-test training window.
# Residuals from the score regression show whether lags 1 and 7
# adequately capture serial dependence in the FPCA scores.
# External features (calendar, holiday, etc.) are not included here:
# we are testing the lag structure of the VAR, not the full spec.

library(fda)

NBASIS      <- 12
NORDER      <- 4
K_SCORE     <- 2
TEST_START  <- as.Date("2023-01-01")
LAG_MAX_ACF <- 30

curves     <- readRDS("data/processed/load_curves.rds")
dates_all  <- as.Date(curves$dates)
hours_all  <- curves$hours
Y_all      <- curves$actual
ok         <- complete.cases(Y_all)
dates_all  <- dates_all[ok]
Y_all      <- Y_all[ok, ]

train_mask  <- dates_all < TEST_START
Y_train     <- Y_all[train_mask, ]
dates_train <- dates_all[train_mask]
n_train     <- nrow(Y_train)

# Fit FPCA on full training window
basis_obj  <- fda::create.bspline.basis(rangeval = range(hours_all),
                                         nbasis   = NBASIS,
                                         norder   = NORDER)
fd_train   <- fda::Data2fd(argvals = hours_all, y = t(Y_train), basisobj = basis_obj)
pca_train  <- fda::pca.fd(fd_train, nharm = K_SCORE)
scores_tr  <- pca_train$scores[, 1:K_SCORE, drop = FALSE]

# Build VAR(1,7) design matrix — rows 8:n_train have both lags
lag_start  <- 8
n_reg      <- n_train - lag_start + 1

Y_scores_reg <- scores_tr[lag_start:n_train,         , drop = FALSE]
X_lag1       <- scores_tr[(lag_start - 1):(n_train - 1), , drop = FALSE]
X_lag7       <- scores_tr[(lag_start - 7):(n_train - 7), , drop = FALSE]

X_var <- cbind(intercept = 1, X_lag1, X_lag7)
colnames(X_var) <- c("intercept",
                     paste0("lag1_s", 1:K_SCORE),
                     paste0("lag7_s", 1:K_SCORE))

# Fit OLS equation by equation, collect residuals
score_resid <- matrix(NA_real_, nrow = n_reg, ncol = K_SCORE)
for (k in seq_len(K_SCORE)) {
  score_resid[, k] <- lm.fit(X_var, Y_scores_reg[, k])$residuals
}

# ACF of score residuals
acf_score_df <- bind_rows(lapply(seq_len(K_SCORE), function(k) {
  series <- score_resid[, k]
  ac     <- acf(series, lag.max = LAG_MAX_ACF, plot = FALSE)
  ci     <- qnorm(0.975) / sqrt(length(series))
  tibble(
    component = paste0("Score ", k),
    lag       = as.integer(ac$lag[-1]),
    acf_val   = as.numeric(ac$acf[-1]),
    upper_ci  =  ci,
    lower_ci  = -ci
  )
}))

p_score_acf <- ggplot(acf_score_df, aes(x = lag, y = acf_val)) +
  geom_hline(aes(yintercept = upper_ci), linetype = "dashed",
             colour = "steelblue", linewidth = 0.4) +
  geom_hline(aes(yintercept = lower_ci), linetype = "dashed",
             colour = "steelblue", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~component, ncol = 1) +
  scale_x_continuous(breaks = seq(0, LAG_MAX_ACF, 7)) +
  labs(
    title    = "ACF of VAR(1,7) score residuals (in-sample)",
    subtitle = paste0("Training: ", min(dates_train), " to ", max(dates_train),
                      "  |  K = ", K_SCORE, "  |  ", n_reg, " obs"),
    x = "Lag (days)",
    y = "ACF"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))

ggsave("output/figs/acf_var_score_residuals.png", p_score_acf,
       width = 7, height = 5, dpi = 150)
cat("Saved: output/figs/acf_var_score_residuals.png\n")

# Cross-component ACF: does u_{d1} predict u_{d2}?
ccf_obj <- ccf(score_resid[, 1], score_resid[, 2],
               lag.max = LAG_MAX_ACF, plot = FALSE)
ci_ccf  <- qnorm(0.975) / sqrt(n_reg)

ccf_df <- tibble(
  lag     = as.integer(ccf_obj$lag),
  ccf_val = as.numeric(ccf_obj$acf),
  upper_ci =  ci_ccf,
  lower_ci = -ci_ccf
)

p_ccf <- ggplot(ccf_df, aes(x = lag, y = ccf_val)) +
  geom_hline(aes(yintercept = upper_ci), linetype = "dashed",
             colour = "steelblue", linewidth = 0.4) +
  geom_hline(aes(yintercept = lower_ci), linetype = "dashed",
             colour = "steelblue", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_x_continuous(breaks = seq(-LAG_MAX_ACF, LAG_MAX_ACF, 7)) +
  labs(
    title = "Cross-correlation of VAR(1,7) score residuals: Score 1 vs Score 2",
    x     = "Lag (days, positive = Score 1 leads)",
    y     = "CCF"
  ) +
  theme_bw(base_size = 11)

ggsave("output/figs/acf_var_score_ccf.png", p_ccf,
       width = 7, height = 3.5, dpi = 150)
cat("Saved: output/figs/acf_var_score_ccf.png\n")

# Ljung-Box for score residuals
lb_scores <- bind_rows(lapply(seq_len(K_SCORE), function(k) {
  lb <- Box.test(score_resid[, k], lag = 14, type = "Ljung-Box")
  tibble(
    component   = paste0("Score ", k),
    n_obs       = n_reg,
    LB_stat     = round(lb$statistic, 3),
    p_value     = round(lb$p.value, 4),
    reject_5pct = lb$p.value < 0.05
  )
}))

cat("\n=== Ljung-Box test (lag=14) for VAR(1,7) score residuals ===\n")
print(lb_scores)
write_csv(lb_scores, "output/tables/acf_var_score_ljung_box.csv")

cat("\nAll ACF diagnostics done.\n")
