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

cat("\nDone.\n")
