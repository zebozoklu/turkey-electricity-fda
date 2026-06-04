# 13_diebold_mariano.R
#
# Diebold-Mariano tests comparing each FPCA model against the official forecast
# and against each other.
#
# The DM test is applied in two ways:
#   (1) Pooled hourly: stack the day x hour panel chronologically and run one
#       DM test on the full hourly error series.
#   (2) Daily aggregated: collapse the panel into a daily mean absolute error
#       series. This is a robustness check that reduces within-day dependence.
#   (3) Hour-by-hour: run a separate DM test for each of the 24 hours and report
#       a table of p-values. Useful for diagnosing where gains are significant.
#
# Loss function: absolute error (consistent with reporting MAE).
# DM statistic uses Harvey, Leybourne & Newbold (1997) small-sample correction.
#
# Requires: forecast (for dm.test), dplyr, tidyr

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(forecast)   # dm.test

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figs",   recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load evaluation files
# ============================================================
# Each eval_long CSV has columns: date, hour, model, actual, forecast,
# error, abs_error, sq_error.
# We load every available eval file and stack them, keeping only models
# we care about.

eval_files <- list(
  "output/tables/fpca_VAR_1_7_eval_long.csv",
  "output/tables/fpca_weekly_VAR_eval_long.csv",
  "output/tables/simple_fpca_weekly_eval_long.csv",
  "output/tables/fpca_VAR_1_7_TLF_regime_eval_long.csv",
  "output/tables/fpca_VAR_1_7_ramp_calendar_eval_long.csv",
  "output/tables/fpca_VAR_1_7_holiday_regime_eval_long.csv",
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_residual_step_K2_eval_long.csv",
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_residual_step_K4_eval_long.csv",
  "output/tables/fpca_RF_fixed_eval_long.csv",
  "output/tables/fpca_RF_rolling_eval_long.csv",
  "output/tables/direct_RF_no_temp_eval_long.csv",
  "output/tables/direct_RF_temp_eval_long.csv"
)

# Read whichever files exist
eval_all <- bind_rows(lapply(eval_files, function(f) {
  if (file.exists(f)) read_csv(f, show_col_types = FALSE, col_types = cols(date = col_character())) else NULL
})) |>
  mutate(date = as.Date(as.character(.data$date), format = "%Y-%m-%d"), hour = as.integer(.data$hour)) |>
  distinct(date, hour, model, .keep_all = TRUE)

cat("Models found in eval data:\n")
print(unique(eval_all$model))
cat("\n")

# ============================================================
# 2. Define the comparison set
# ============================================================
# Benchmark = "Official forecast"
# Challengers = all other models except seasonal naive

BENCHMARK   <- "Official forecast"
NAIVE       <- "Seasonal naive: Y[d-7]"

challengers <- setdiff(unique(eval_all$model), c(BENCHMARK, NAIVE))
cat("Challengers:\n"); print(challengers); cat("\n")

# ============================================================
# 3. Helper: run DM test (Harvey-Leybourne-Newbold corrected)
# ============================================================
# e1 = errors of model 1 (challenger), e2 = errors of benchmark
# loss = "AE" (absolute error) or "SE" (squared error)
# alternative = "less" means H1: model 1 has SMALLER loss than model 2

run_dm <- function(e1, e2, loss = "AE", h = 1, alternative = "less") {
  if (loss == "AE") {
    d <- abs(e1) - abs(e2)
  } else {
    d <- e1^2 - e2^2
  }
  
  # Remove NAs
  ok <- complete.cases(d)
  d  <- d[ok]
  
  if (length(d) < 10) return(list(statistic = NA, p.value = NA, mean_d = NA))
  
  result <- tryCatch(
    dm.test(e1[ok], e2[ok], alternative = alternative,
            h = h, power = if (loss == "AE") 1 else 2),
    error = function(e) list(statistic = NA_real_, p.value = NA_real_)
  )
  
  list(
    statistic = as.numeric(result$statistic),
    p.value   = result$p.value,
    mean_d    = mean(d)   # negative = challenger is better
  )
}

# ============================================================
# 4. Pooled hourly DM test
# ============================================================
# Stack all date-hour forecast errors in chronological order and run one DM test
# per challenger-official pair. Negative mean_d means the challenger has lower
# mean absolute error than the official forecast.

run_pooled_hourly_dm <- function(challenger_name) {
  bench_series <- eval_all |>
    filter(model == BENCHMARK) |>
    arrange(date, hour) |>
    select(date, hour, bench_err = error)

  chall_series <- eval_all |>
    filter(model == challenger_name) |>
    arrange(date, hour) |>
    select(date, hour, chall_err = error)

  joined <- inner_join(bench_series, chall_series, by = c("date", "hour")) |>
    arrange(date, hour)

  dm_res <- tryCatch(
    dm.test(
      joined$chall_err,
      joined$bench_err,
      alternative = "less",
      h = 1,
      power = 1
    ),
    error = function(e) list(statistic = NA_real_, p.value = NA_real_)
  )

  tibble(
    challenger = challenger_name,
    benchmark = BENCHMARK,
    n_obs = nrow(joined),
    mean_d = mean(abs(joined$chall_err) - abs(joined$bench_err), na.rm = TRUE),
    DM_stat = as.numeric(dm_res$statistic),
    p_value = dm_res$p.value,
    sig_10pct = !is.na(dm_res$p.value) & dm_res$p.value < 0.10,
    sig_5pct = !is.na(dm_res$p.value) & dm_res$p.value < 0.05,
    sig_1pct = !is.na(dm_res$p.value) & dm_res$p.value < 0.01
  )
}

pooled_results <- bind_rows(lapply(challengers, run_pooled_hourly_dm))

cat("=== Pooled hourly DM test (H1: challenger < official) ===\n")
print(pooled_results, width = 120)
cat("\n")

write_csv(pooled_results, "output/tables/dm_test_pooled.csv")

# ============================================================
# 5. Daily aggregated DM test
# ============================================================
# For each day, compute mean absolute error across 24 hours.
# Then run DM on this daily series (length ~ number of evaluation days).
# This reduces the within-day correlation problem.

daily_loss <- eval_all |>
  group_by(date, model) |>
  summarise(
    daily_mae = mean(abs_error, na.rm = TRUE),
    daily_mse = mean(sq_error,  na.rm = TRUE),
    .groups = "drop"
  )

run_daily_dm <- function(challenger_name) {
  bench_series <- daily_loss |>
    filter(model == BENCHMARK) |>
    arrange(date) |>
    select(date, bench_mae = daily_mae)
  
  chall_series <- daily_loss |>
    filter(model == challenger_name) |>
    arrange(date) |>
    select(date, chall_mae = daily_mae)
  
  joined <- inner_join(bench_series, chall_series, by = "date")
  
  # d_t = loss(challenger) - loss(benchmark)
  # H1: challenger < benchmark  =>  alternative = "less"
  n  <- nrow(joined)
  d  <- joined$chall_mae - joined$bench_mae
  
  dm_res <- tryCatch(
    dm.test(
      joined$chall_mae,
      joined$bench_mae,
      alternative = "less",   # challenger is better
      h = 1,
      power = 1               # MAE
    ),
    error = function(e) list(statistic = NA_real_, p.value = NA_real_)
  )
  
  tibble(
    challenger    = challenger_name,
    benchmark     = BENCHMARK,
    n_days        = n,
    mean_d        = mean(d, na.rm = TRUE),   # negative = challenger better
    DM_stat       = as.numeric(dm_res$statistic),
    p_value       = dm_res$p.value,
    sig_10pct     = !is.na(dm_res$p.value) & dm_res$p.value < 0.10,
    sig_5pct      = !is.na(dm_res$p.value) & dm_res$p.value < 0.05,
    sig_1pct      = !is.na(dm_res$p.value) & dm_res$p.value < 0.01
  )
}

daily_results <- bind_rows(lapply(challengers, run_daily_dm))

cat("=== Daily aggregated DM test (daily MAE series, H1: challenger < official) ===\n")
print(daily_results, width = 120)
cat("\n")

write_csv(daily_results, "output/tables/dm_test_daily.csv")

# ============================================================
# 6. Hour-by-hour DM test
# ============================================================
# For each hour h in 0:23, extract the time series of hourly errors
# and run DM. Reports a 24-row table for each challenger.

run_hourly_dm <- function(challenger_name) {
  bench_df <- eval_all |>
    filter(model == BENCHMARK) |>
    select(date, hour, bench_err = error)
  
  chall_df <- eval_all |>
    filter(model == challenger_name) |>
    select(date, hour, chall_err = error)
  
  joined <- inner_join(bench_df, chall_df, by = c("date", "hour"))
  
  bind_rows(lapply(0:23, function(h) {
    sub <- joined |> filter(hour == h) |> arrange(date)
    
    if (nrow(sub) < 10) {
      return(tibble(hour = h, DM_stat = NA, p_value = NA, mean_d = NA))
    }
    
    dm_res <- tryCatch(
      dm.test(sub$chall_err, sub$bench_err,
              alternative = "less", h = 1, power = 1),
      error = function(e) list(statistic = NA_real_, p.value = NA_real_)
    )
    
    tibble(
      hour      = h,
      n         = nrow(sub),
      mean_d    = mean(abs(sub$chall_err) - abs(sub$bench_err), na.rm = TRUE),
      DM_stat   = as.numeric(dm_res$statistic),
      p_value   = dm_res$p.value
    )
  })) |>
    mutate(
      challenger = challenger_name,
      benchmark  = BENCHMARK,
      sig_5pct   = !is.na(p_value) & p_value < 0.05
    )
}

hourly_results <- bind_rows(lapply(challengers, run_hourly_dm))

cat("=== Hour-by-hour DM tests ===\n")
print(hourly_results |> select(challenger, hour, DM_stat, p_value, sig_5pct),
      n = 200, width = 120)

write_csv(hourly_results, "output/tables/dm_test_by_hour.csv")

# ============================================================
# 7. Summary p-value heatmap (text-friendly table)
# ============================================================

pval_wide <- hourly_results |>
  select(challenger, hour, p_value) |>
  pivot_wider(names_from = hour, values_from = p_value,
              names_prefix = "h")

cat("\n=== P-value table (rows=model, cols=hour) ===\n")
print(pval_wide, width = 200)
write_csv(pval_wide, "output/tables/dm_test_pvalue_by_hour_wide.csv")

# ============================================================
# 8. Count of hours where challenger significantly beats official
# ============================================================

sig_summary <- hourly_results |>
  group_by(challenger) |>
  summarise(
    hours_sig_5pct  = sum(sig_5pct, na.rm = TRUE),
    hours_better    = sum(mean_d < 0, na.rm = TRUE),   # lower abs error
    hours_worse     = sum(mean_d > 0, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Hours where challenger significantly (5%) beats official ===\n")
print(sig_summary, width = 120)
write_csv(sig_summary, "output/tables/dm_test_sig_hours_summary.csv")

# ============================================================
# 9. Optional: pairwise DM between challengers
# ============================================================
# Compares your best model vs each other challenger.
# Uses the best model in the pooled hourly results (lowest mean_d) as reference.

if (nrow(pooled_results) > 1) {
  best_model <- pooled_results |>
    filter(!is.na(DM_stat)) |>
    slice_min(mean_d, n = 1) |>
    pull(challenger)
  
  cat("\nBest challenger model:", best_model, "\n")
  
  other_challengers <- setdiff(challengers, best_model)
  
  run_pairwise_dm <- function(rival_name) {
    best_hourly <- eval_all |>
      filter(model == best_model) |>
      arrange(date, hour) |>
      select(date, hour, best_err = error)
    
    rival_hourly <- eval_all |>
      filter(model == rival_name) |>
      arrange(date, hour) |>
      select(date, hour, rival_err = error)
    
    joined <- inner_join(best_hourly, rival_hourly, by = c("date", "hour")) |>
      arrange(date, hour)
    
    if (nrow(joined) < 10) return(NULL)
    
    dm_res <- tryCatch(
      dm.test(joined$best_err, joined$rival_err,
              alternative = "less", h = 1, power = 1),
      error = function(e) list(statistic = NA_real_, p.value = NA_real_)
    )
    
    tibble(
      model_1   = best_model,
      model_2   = rival_name,
      n_obs     = nrow(joined),
      mean_d    = mean(abs(joined$best_err) - abs(joined$rival_err), na.rm = TRUE),
      DM_stat   = as.numeric(dm_res$statistic),
      p_value   = dm_res$p.value
    )
  }
  
  pairwise_results <- bind_rows(lapply(other_challengers, run_pairwise_dm))
  
  cat("\n=== Pairwise DM: best model vs others (H1: best < rival) ===\n")
  print(pairwise_results, width = 120)
  write_csv(pairwise_results, "output/tables/dm_test_pairwise.csv")
}

cat("\nDone. Tables saved to output/tables/dm_test_*.csv\n")
