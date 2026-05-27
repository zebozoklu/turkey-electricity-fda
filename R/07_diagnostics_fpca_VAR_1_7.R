# 07_diagnostics_fpca_VAR_1_7.R

library(dplyr)
library(ggplot2)
library(lubridate)

dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

eval_df <- read.csv("output/tables/fpca_VAR_1_7_eval_long.csv")

eval_df <- eval_df |>
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    ym = format(date, "%Y-%m"),
    hour = as.integer(hour)
  )

# ------------------------------------------------------------
# Overall table with official comparison
# ------------------------------------------------------------

overall <- eval_df |>
  group_by(model) |>
  summarise(
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

official_mae <- overall$mae[overall$model == "Official forecast"]
official_rmse <- overall$rmse[overall$model == "Official forecast"]

overall <- overall |>
  mutate(
    mae_improvement_vs_official_pct =
      100 * (official_mae - mae) / official_mae,
    rmse_improvement_vs_official_pct =
      100 * (official_rmse - rmse) / official_rmse
  )

print(overall)

write.csv(
  overall,
  "output/tables/fpca_VAR_1_7_overall_vs_official.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# By year
# ------------------------------------------------------------

by_year <- eval_df |>
  group_by(year, model) |>
  summarise(
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

print(by_year)

write.csv(
  by_year,
  "output/tables/fpca_VAR_1_7_by_year.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# By month
# ------------------------------------------------------------

by_month <- eval_df |>
  group_by(ym, model) |>
  summarise(
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  by_month,
  "output/tables/fpca_VAR_1_7_by_month.csv",
  row.names = FALSE
)

p_month_mae <- by_month |>
  ggplot(aes(x = ym, y = mae, group = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  labs(
    x = "Month",
    y = "MAE",
    title = "Forecast MAE by month",
    linetype = "Model"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(
  "output/figs/fpca_VAR_1_7_mae_by_month.png",
  p_month_mae,
  width = 10,
  height = 4.8
)

# ------------------------------------------------------------
# By hour
# ------------------------------------------------------------

by_hour <- eval_df |>
  group_by(hour, model) |>
  summarise(
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    bias = mean(error, na.rm = TRUE),
    .groups = "drop"
  )

print(by_hour)

write.csv(
  by_hour,
  "output/tables/fpca_VAR_1_7_by_hour.csv",
  row.names = FALSE
)

p_hour_mae <- by_hour |>
  ggplot(aes(x = hour, y = mae, group = model, linetype = model)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "MAE",
    title = "Forecast MAE by hour",
    linetype = "Model"
  )

ggsave(
  "output/figs/fpca_VAR_1_7_mae_by_hour_diagnostic.png",
  p_hour_mae,
  width = 8,
  height = 4.8
)

# ------------------------------------------------------------
# Direct FPCA vs Official differences: robust version
# ------------------------------------------------------------

print(unique(eval_df$model))

official_err <- eval_df |>
  filter(model == "Official forecast") |>
  dplyr::select(date, hour, abs_error_official = abs_error, sq_error_official = sq_error)

fpca_err <- eval_df |>
  filter(model == "FPCA VAR(1,7) scores") |>
  dplyr::select(date, hour, abs_error_fpca = abs_error, sq_error_fpca = sq_error)

wide_model <- official_err |>
  inner_join(fpca_err, by = c("date", "hour")) |>
  mutate(
    abs_gain_vs_official = abs_error_official - abs_error_fpca,
    sq_gain_vs_official = sq_error_official - sq_error_fpca
  )

# Positive gain means FPCA is better than official
gain_by_hour <- wide_model |>
  group_by(hour) |>
  summarise(
    mae_gain_vs_official = mean(abs_gain_vs_official, na.rm = TRUE),
    mse_gain_vs_official = mean(sq_gain_vs_official, na.rm = TRUE),
    .groups = "drop"
  )

print(gain_by_hour)

write.csv(
  gain_by_hour,
  "output/tables/fpca_VAR_1_7_gain_vs_official_by_hour.csv",
  row.names = FALSE
)

p_gain_hour <- gain_by_hour |>
  ggplot(aes(x = hour, y = mae_gain_vs_official)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  labs(
    x = "Hour",
    y = "MAE gain vs official",
    title = "Positive values mean FPCA beats official forecast"
  )

ggsave(
  "output/figs/fpca_VAR_1_7_gain_vs_official_by_hour.png",
  p_gain_hour,
  width = 8,
  height = 4.8
)

gain_by_month <- wide_model |>
  mutate(ym = format(date, "%Y-%m")) |>
  group_by(ym) |>
  summarise(
    mae_gain_vs_official = mean(abs_gain_vs_official, na.rm = TRUE),
    mse_gain_vs_official = mean(sq_gain_vs_official, na.rm = TRUE),
    .groups = "drop"
  )

print(gain_by_month)

write.csv(
  gain_by_month,
  "output/tables/fpca_VAR_1_7_gain_vs_official_by_month.csv",
  row.names = FALSE
)

p_gain_month <- gain_by_month |>
  ggplot(aes(x = ym, y = mae_gain_vs_official, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  labs(
    x = "Month",
    y = "MAE gain vs official",
    title = "Positive values mean FPCA beats official forecast"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(
  "output/figs/fpca_VAR_1_7_gain_vs_official_by_month.png",
  p_gain_month,
  width = 10,
  height = 4.8
)


