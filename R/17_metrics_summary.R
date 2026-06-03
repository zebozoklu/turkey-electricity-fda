# 17_metrics_summary.R
#
# Unified metrics table across all models.
# Reads eval_long CSVs — no forecasting scripts need to be rerun.
# Reports MAE, RMSE, MAPE, Bias, and gains vs official forecast.
#
# MAPE = mean(|error| / actual) * 100
# Requires actual > 0; electricity consumption satisfies this.

library(dplyr)
library(tidyr)
library(readr)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load all eval_long files
# ============================================================

eval_files <- c(
  "output/tables/simple_fpca_weekly_eval_long.csv",
  "output/tables/fpca_weekly_VAR_eval_long.csv",
  "output/tables/fpca_VAR_1_7_eval_long.csv",
  "output/tables/fpca_VAR_1_7_TLF_regime_eval_long.csv",
  "output/tables/fpca_VAR_1_7_ramp_calendar_eval_long.csv",
  "output/tables/fpca_VAR_1_7_holiday_regime_eval_long.csv",
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_residual_step_eval_long.csv",
  "output/tables/fpca_RF_eval_long.csv"
)

eval_all <- bind_rows(lapply(eval_files, function(f) {
  if (!file.exists(f)) { message("Skipping (not found): ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE,
           col_types = cols(date = col_character())) |>
    mutate(date = as.Date(date), hour = as.integer(hour))
})) |>
  distinct(date, hour, model, .keep_all = TRUE)

# Add official forecast rows from the main eval file if not already present
official_file <- "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv"
if (file.exists(official_file)) {
  official_rows <- read_csv(official_file, show_col_types = FALSE,
                            col_types = cols(date = col_character())) |>
    mutate(date = as.Date(date), hour = as.integer(hour)) |>
    filter(model == "Official forecast")
  eval_all <- bind_rows(eval_all, official_rows) |>
    distinct(date, hour, model, .keep_all = TRUE)
}

cat("Models in eval data:\n")
print(unique(eval_all$model))
cat("\n")

# ============================================================
# 2. Compute metrics
# ============================================================

metrics <- eval_all |>
  filter(actual > 0) |>           # MAPE guard (electricity: always true)
  group_by(model) |>
  summarise(
    n        = n(),
    mae      = mean(abs_error,            na.rm = TRUE),
    rmse     = sqrt(mean(sq_error,        na.rm = TRUE)),
    mape     = mean(abs_error / actual,   na.rm = TRUE) * 100,
    bias     = mean(error,                na.rm = TRUE),
    .groups  = "drop"
  )

# Gains relative to official forecast
off <- metrics |> filter(model == "Official forecast")

if (nrow(off) == 1) {
  metrics <- metrics |>
    mutate(
      mae_gain  = 100 * (off$mae  - mae)  / off$mae,
      rmse_gain = 100 * (off$rmse - rmse) / off$rmse,
      mape_gain = 100 * (off$mape - mape) / off$mape
    )
}

# ============================================================
# 3. Model ordering for display
# ============================================================

model_order <- c(
  "Seasonal naive: Y[d-7]",
  "Official forecast",
  "Simple FPCA weekly AR",
  "FPCA weekly VAR",
  "FPCA VAR(1,7) scores",
  "FPCA VAR(1,7) + TLF regime",
  "FPCA VAR(1,7) + ramp/calendar",
  "FPCA VAR(1,7) + holiday/load-shape",
  "FPCA VAR(1,7) + holiday/load/derivative regime",
  "FPCA RF score regression",
  "Base model + residual second step",
  "Base model (K=2) + residual step"
)

metrics <- metrics |>
  mutate(model = factor(model, levels = c(model_order,
                                           setdiff(model, model_order)))) |>
  arrange(model)

# ============================================================
# 4. Print and save
# ============================================================

cat("=== Full metrics table ===\n\n")
metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = 50, width = 140)

write_csv(metrics, "output/tables/metrics_all_models.csv")
cat("\nSaved: output/tables/metrics_all_models.csv\n")

# ============================================================
# 5. Thesis-facing table: main model path only
# ============================================================

main_path <- c(
  "Seasonal naive: Y[d-7]",
  "Official forecast",
  "FPCA VAR(1,7) scores",
  "FPCA VAR(1,7) + ramp/calendar",
  "FPCA VAR(1,7) + holiday/load-shape",
  "FPCA VAR(1,7) + holiday/load/derivative regime",
  "FPCA RF score regression"
)

thesis_table <- metrics |>
  filter(model %in% main_path) |>
  select(model, mae, rmse, mape, bias, mae_gain, rmse_gain, mape_gain)

cat("\n=== Thesis main-path table ===\n\n")
thesis_table |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = 20, width = 140)

write_csv(thesis_table, "output/tables/metrics_thesis_table.csv")
cat("Saved: output/tables/metrics_thesis_table.csv\n")
