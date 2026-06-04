# 16e_forecast_direct_RF.R
#
# Direct multi-output Random Forest benchmark.
#
# WHY THIS IS DIFFERENT FROM SCRIPTS 16 / 16b / 16d
# --------------------------------------------------
# Scripts 16/16b/16d predict FPCA scores and reconstruct via
#   Ŷ = μ + φ₁·RF₁(X) + φ₂·RF₂(X)
# This is nonlinear in X but constrained to a K=2 FPCA subspace that
# captures 97% of curve variance. OLS is equally constrained. So those
# scripts correctly test "can RF find nonlinearity in score dynamics?"
# but cannot test nonlinearity in the remaining 3% of variance, and
# both models are handicapped identically by the subspace.
#
# THIS script bypasses FPCA entirely and predicts all 24 hours jointly
# using ranger's multi-output regression. The forecast is completely
# unconstrained — it can lie anywhere in R^24. This is the general
# nonlinearity test: can RF beat FPCA-OLS when it is free to exploit
# the full curve space?
#
# TWO VARIANTS run back-to-back:
#   (A) without temperature  → tests pure structural nonlinearity
#   (B) with temperature     → adds HDD/CDD/rolling temp as features
#
# Comparing (A) vs (B) isolates the marginal value of temperature.
# Comparing both vs FPCA-OLS answers: is the OLS spec genuinely adequate
# or does a free nonlinear model with richer features do better?
#
# FEATURE MATRIX (per training/test day d)
#   lag-1 curve      Y[d-1, 0:23]         24 features
#   lag-7 curve      Y[d-7, 0:23]         24 features
#   calendar                               14 features
#   temperature (B only)  temp_mean, hdd, cdd, roll_temp_7   4 features
#
# SPEED: multi-output ranger (one fit per day, 24-dim response).
# With 100 trees and a 2023-only test window this runs in ~15-25 min.

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(readr)
library(ranger)

dir.create("output/figs",    recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Settings
# ============================================================

TEST_START     <- as.Date("2023-01-01")
TEST_END       <- as.Date("2023-12-31")
MIN_TRAIN_DAYS <- 365
MIN_LAG_DAYS   <- 7          # only need lag-7, not lag-28
RF_TREES       <- 100
RF_SEED        <- 42

MODEL_NO_TEMP  <- "Direct RF (no temperature)"
MODEL_TEMP     <- "Direct RF (with temperature)"
OLS_NAME       <- "FPCA VAR(1,7) + holiday/load/derivative regime"

WEATHER_FILE   <- "data/processed/weather_daily.rds"
HAS_WEATHER    <- file.exists(WEATHER_FILE)

if (!HAS_WEATHER) {
  warning("Weather file not found at ", WEATHER_FILE,
          ". Run 00b_download_weather.R first. Running without temperature.")
}

# ============================================================
# Load curves
# ============================================================

curves <- readRDS("data/processed/load_curves.rds")

dates <- as.Date(curves$dates)
hours <- curves$hours
Y     <- curves$actual
colnames(Y) <- as.character(hours)

has_official <- !is.null(curves$forecast)
if (has_official) {
  Y_official <- curves$forecast
  colnames(Y_official) <- as.character(hours)
}

ok    <- complete.cases(Y)
dates <- dates[ok]
Y     <- Y[ok, , drop = FALSE]
if (has_official) Y_official <- Y_official[ok, , drop = FALSE]

n_days <- nrow(Y)
n_hours <- length(hours)

cat("Data loaded:", n_days, "days,", n_hours, "hours\n")

# ============================================================
# Calendar features  (same as all other scripts)
# ============================================================

date_seq <- function(start, end) seq(as.Date(start), as.Date(end), by = "day")

make_holiday_lookup <- function(years) {
  fixed_holidays <- bind_rows(lapply(years, function(y) {
    tibble(
      date = as.Date(c(
        paste0(y, "-01-01"), paste0(y, "-04-23"), paste0(y, "-05-01"),
        paste0(y, "-05-19"), paste0(y, "-07-15"), paste0(y, "-08-30"),
        paste0(y, "-10-29")
      )),
      holiday_type = "fixed"
    )
  }))
  religious_ranges <- tribble(
    ~start, ~end, ~holiday_type,
    "2019-06-04","2019-06-06","ramadan_bayram","2020-05-24","2020-05-26","ramadan_bayram",
    "2021-05-13","2021-05-15","ramadan_bayram","2022-05-02","2022-05-04","ramadan_bayram",
    "2023-04-21","2023-04-23","ramadan_bayram","2024-04-10","2024-04-12","ramadan_bayram",
    "2025-03-30","2025-04-01","ramadan_bayram",
    "2019-08-11","2019-08-14","kurban_bayram","2020-07-31","2020-08-03","kurban_bayram",
    "2021-07-20","2021-07-23","kurban_bayram","2022-07-09","2022-07-12","kurban_bayram",
    "2023-06-28","2023-07-01","kurban_bayram","2024-06-16","2024-06-19","kurban_bayram",
    "2025-06-06","2025-06-09","kurban_bayram"
  ) |> mutate(start = as.Date(start), end = as.Date(end))
  religious_holidays <- bind_rows(lapply(seq_len(nrow(religious_ranges)), function(i)
    tibble(date = date_seq(religious_ranges$start[i], religious_ranges$end[i]),
           holiday_type = religious_ranges$holiday_type[i])))
  bind_rows(fixed_holidays, religious_holidays) |>
    filter(year(date) %in% years) |> distinct(date, .keep_all = TRUE)
}

make_ramadan_lookup <- function(years) {
  rr <- tribble(~start, ~end,
    "2019-05-06","2019-06-03","2020-04-24","2020-05-23","2021-04-13","2021-05-12",
    "2022-04-02","2022-05-01","2023-03-23","2023-04-20","2024-03-11","2024-04-09",
    "2025-03-01","2025-03-29") |> mutate(start = as.Date(start), end = as.Date(end))
  bind_rows(lapply(seq_len(nrow(rr)), function(i)
    tibble(date = date_seq(rr$start[i], rr$end[i])))) |>
    filter(year(date) %in% years) |> distinct(date)
}

holiday_lookup <- make_holiday_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))
ramadan_lookup <- make_ramadan_lookup(seq(min(year(dates)) - 1, max(year(dates)) + 1))

make_calendar_features <- function(date_vec) {
  doy <- yday(date_vec); dow <- wday(date_vec, week_start = 1)
  mn  <- month(date_vec); hd  <- holiday_lookup$date
  tibble(date = date_vec) |>
    mutate(
      dow = dow, weekend = as.integer(dow >= 6),
      monday = as.integer(dow == 1), friday = as.integer(dow == 5),
      month = mn,
      sin_week = sin(2 * pi * dow / 7), cos_week = cos(2 * pi * dow / 7),
      sin_year = sin(2 * pi * doy / 365.25), cos_year = cos(2 * pi * doy / 365.25),
      summer_school_break = as.integer(
        (mn == 6 & mday(date) >= 15) | mn %in% c(7, 8) | (mn == 9 & mday(date) <= 15)),
      is_holiday           = as.integer(date %in% hd),
      is_religious_holiday = as.integer(date %in% holiday_lookup$date[
        holiday_lookup$holiday_type %in% c("ramadan_bayram", "kurban_bayram")]),
      is_ramadan   = as.integer(date %in% ramadan_lookup$date),
      holiday_eve  = as.integer((date + days(1)) %in% hd),
      post_holiday = as.integer((date - days(1)) %in% hd),
      near_holiday = as.integer(date %in% c(
        hd - days(2), hd - days(1), hd, hd + days(1), hd + days(2))),
      bridge_day   = as.integer(dow %in% c(1, 5) & !(date %in% hd) &
        ((date - days(1)) %in% hd | (date + days(1)) %in% hd))
    ) |> dplyr::select(-date)
}

# ============================================================
# Weather features
# ============================================================

if (HAS_WEATHER) {
  weather <- readRDS(WEATHER_FILE) |>
    dplyr::select(date, temp_mean, hdd, cdd, roll_temp_7) |>
    filter(!is.na(temp_mean)) |>
    # Shift ALL weather features by 1 day: every feature for day d uses
    # information from day d-1 or earlier — strictly no temperature leakage.
    # temp_lag1        = temp_mean[d-1]  (yesterday's temperature)
    # hdd_lag1         = hdd[d-1]
    # cdd_lag1         = cdd[d-1]
    # roll_temp_7_lag1 = mean(temp[d-7:d-1])  (7-day trailing mean ending yesterday)
    mutate(
      temp_lag1        = lag(temp_mean,   1),
      hdd_lag1         = lag(hdd,         1),
      cdd_lag1         = lag(cdd,         1),
      roll_temp_7_lag1 = lag(roll_temp_7, 1)
    ) |>
    dplyr::select(date, temp_lag1, hdd_lag1, cdd_lag1, roll_temp_7_lag1)
  cat("Weather loaded:", nrow(weather), "days,",
      as.character(min(weather$date)), "to", as.character(max(weather$date)), "\n")
} else {
  weather <- NULL
}

make_weather_features <- function(date_vec) {
  if (is.null(weather)) return(NULL)
  tibble(date = date_vec) |>
    left_join(weather, by = "date") |>
    dplyr::select(-date)
}

# ============================================================
# Precompute feature matrix for all days
# ============================================================
# Each row d (d >= 8) has: lag-1 curve + lag-7 curve + calendar + weather
# Lag features use only past observations — no leakage.

cat("Precomputing feature matrix...\n")

d_all <- (MIN_LAG_DAYS + 1):n_days   # days with valid lag-7

# Lag curves (columns named h0_lag1 ... h23_lag1, etc.)
lag1_mat <- Y[d_all - 1, , drop = FALSE]
lag7_mat <- Y[d_all - 7, , drop = FALSE]

colnames(lag1_mat) <- paste0("h", hours, "_lag1")
colnames(lag7_mat) <- paste0("h", hours, "_lag7")

# Calendar
cal_mat <- make_calendar_features(dates[d_all])

# Weather (may be NULL if file missing)
wthr_mat <- make_weather_features(dates[d_all])

# Combine (without temperature first; we'll add it optionally inside the loop)
if (!is.null(wthr_mat)) {
  X_all <- cbind(lag1_mat, lag7_mat, as.matrix(cal_mat), as.matrix(wthr_mat))
  n_temp_features <- ncol(wthr_mat)
  base_features <- ncol(lag1_mat) + ncol(lag7_mat) + ncol(cal_mat)
  cat("Features: lag-1(", ncol(lag1_mat), ") + lag-7(", ncol(lag7_mat),
      ") + calendar(", ncol(cal_mat), ") + weather(", n_temp_features, ")\n")
} else {
  X_all <- cbind(lag1_mat, lag7_mat, as.matrix(cal_mat))
  n_temp_features <- 0
  base_features <- ncol(X_all)
  cat("Features: lag-1(", ncol(lag1_mat), ") + lag-7(", ncol(lag7_mat),
      ") + calendar(", ncol(cal_mat), ") [no weather]\n")
}

# Response matrix (same indexing as d_all)
Y_all <- Y[d_all, , drop = FALSE]
colnames(Y_all) <- paste0("y_h", hours)

# Map from global day index to row in X_all / Y_all
day_to_row <- function(day_idx) match(day_idx, d_all)

# ============================================================
# Test indices
# ============================================================

test_idx <- which(dates >= TEST_START & dates <= TEST_END)
test_idx <- test_idx[test_idx > max(MIN_LAG_DAYS, MIN_TRAIN_DAYS)]

cat("\nForecast days:", length(test_idx),
    "| First:", as.character(min(dates[test_idx])),
    "| Last:", as.character(max(dates[test_idx])), "\n\n")

# ============================================================
# Rolling-origin forecast function
# ============================================================
# For each test day i:
#   - Training rows: all d in d_all with d < i
#   - For variant A (no temp): drop weather columns
#   - For variant B (with temp): use all columns

forecast_one_day <- function(i, use_temp = FALSE) {

  # Training rows: days strictly before i that have valid lag features
  train_rows <- which(d_all < i)
  if (length(train_rows) < 50) return(NULL)

  if (use_temp && n_temp_features > 0) {
    X_tr <- X_all[train_rows, , drop = FALSE]
    x_new <- X_all[day_to_row(i), , drop = FALSE]
  } else {
    # Drop weather columns
    keep_cols <- seq_len(base_features)
    X_tr  <- X_all[train_rows, keep_cols, drop = FALSE]
    x_new <- X_all[day_to_row(i), keep_cols, drop = FALSE]
  }

  Y_tr <- Y_all[train_rows, , drop = FALSE]

  # Handle any NAs in weather features
  complete_rows <- complete.cases(X_tr) & complete.cases(Y_tr)
  X_tr <- X_tr[complete_rows, , drop = FALSE]
  Y_tr <- Y_tr[complete_rows, , drop = FALSE]

  if (nrow(X_tr) < 50) return(NULL)

  # Fit one single-output RF per hour (24 models).
  # More robust than multi-output formula: avoids ranger's column-name
  # mismatch and competing-risks misclassification on matrix y.
  X_tr_df  <- as.data.frame(X_tr)
  x_new_df <- as.data.frame(x_new)
  # Guarantee identical column names between train and predict
  colnames(x_new_df) <- colnames(X_tr_df)

  pred_vec <- numeric(n_hours)
  for (h in seq_len(n_hours)) {
    df_tr <- cbind(y = Y_tr[, h], X_tr_df)
    rf_h  <- ranger(
      y ~ .,
      data        = df_tr,
      num.trees   = RF_TREES,
      seed        = RF_SEED,
      num.threads = 1,
      verbose     = FALSE
    )
    pred_vec[h] <- predict(rf_h, data = x_new_df)$predictions
  }
  pred_vec
}

# ============================================================
# Run variant A: no temperature
# ============================================================

cat("=== Variant A: Direct RF without temperature ===\n")
pred_no_temp <- matrix(NA_real_, nrow = length(test_idx), ncol = n_hours,
                       dimnames = list(as.character(dates[test_idx]),
                                       as.character(hours)))
t0 <- proc.time()

for (j in seq_along(test_idx)) {
  i <- test_idx[j]
  out <- forecast_one_day(i, use_temp = FALSE)
  if (!is.null(out)) pred_no_temp[j, ] <- out

  if (j %% 25 == 0) {
    el  <- round((proc.time() - t0)["elapsed"])
    rem <- round(el / j * (length(test_idx) - j))
    message("  A: ", j, "/", length(test_idx),
            "  date: ", dates[i],
            "  elapsed: ", el, "s  est remaining: ", rem, "s")
  }
}
cat("Variant A done.", round((proc.time() - t0)["elapsed"]), "s\n\n")

# ============================================================
# Run variant B: with temperature  (only if weather available)
# ============================================================

if (HAS_WEATHER && n_temp_features > 0) {
  cat("=== Variant B: Direct RF with temperature ===\n")
  pred_temp <- matrix(NA_real_, nrow = length(test_idx), ncol = n_hours,
                      dimnames = list(as.character(dates[test_idx]),
                                      as.character(hours)))
  t0 <- proc.time()

  for (j in seq_along(test_idx)) {
    i <- test_idx[j]
    out <- forecast_one_day(i, use_temp = TRUE)
    if (!is.null(out)) pred_temp[j, ] <- out

    if (j %% 25 == 0) {
      el  <- round((proc.time() - t0)["elapsed"])
      rem <- round(el / j * (length(test_idx) - j))
      message("  B: ", j, "/", length(test_idx),
              "  date: ", dates[i],
              "  elapsed: ", el, "s  est remaining: ", rem, "s")
    }
  }
  cat("Variant B done.", round((proc.time() - t0)["elapsed"]), "s\n\n")
} else {
  pred_temp <- NULL
  cat("Skipping variant B (no weather file).\n\n")
}

# ============================================================
# Evaluation
# ============================================================

actual_test <- Y[test_idx, , drop = FALSE]
pred_naive  <- Y[test_idx - 7, , drop = FALSE]
if (has_official) pred_official <- Y_official[test_idx, , drop = FALSE]

# OLS for the same sub-period
ols_eval <- read_csv(
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  show_col_types = FALSE, col_types = cols(date = col_character())
) |>
  mutate(date = as.Date(date), hour = as.integer(hour)) |>
  filter(model == OLS_NAME, date %in% dates[test_idx])

make_eval_df <- function(actual, pred, model_name) {
  tibble(
    date     = rep(dates[test_idx], each = n_hours),
    hour     = rep(hours, times = length(test_idx)),
    model    = model_name,
    actual   = as.vector(t(actual)),
    forecast = as.vector(t(pred))
  ) |>
    mutate(
      year = year(date), month = month(date), ym = format(date, "%Y-%m"),
      error = actual - forecast, abs_error = abs(error), sq_error = error^2
    )
}

eval_df <- bind_rows(
  make_eval_df(actual_test, pred_no_temp, MODEL_NO_TEMP),
  make_eval_df(actual_test, pred_naive,   "Seasonal naive: Y[d-7]")
)
if (!is.null(pred_temp))
  eval_df <- bind_rows(eval_df, make_eval_df(actual_test, pred_temp, MODEL_TEMP))
if (has_official)
  eval_df <- bind_rows(eval_df, make_eval_df(actual_test, pred_official, "Official forecast"))

eval_combined <- bind_rows(
  eval_df,
  ols_eval |> dplyr::select(all_of(names(eval_df)))
)

# ============================================================
# Metrics
# ============================================================

metrics <- eval_combined |>
  filter(actual > 0) |>
  group_by(model) |>
  summarise(
    n    = n(),
    mae  = mean(abs_error,          na.rm = TRUE),
    rmse = sqrt(mean(sq_error,      na.rm = TRUE)),
    mape = mean(abs_error / actual, na.rm = TRUE) * 100,
    bias = mean(error,              na.rm = TRUE),
    .groups = "drop"
  )

off <- metrics |> filter(model == "Official forecast")
if (nrow(off) == 1) {
  metrics <- metrics |>
    mutate(
      mae_gain  = 100 * (off$mae  - mae)  / off$mae,
      rmse_gain = 100 * (off$rmse - rmse) / off$rmse,
      mape_gain = 100 * (off$mape - mape) / off$mape
    )
}

cat("=== Overall metrics (", as.character(TEST_START), "to",
    as.character(TEST_END), ") ===\n")
print(metrics |>
        dplyr::select(model, mae, rmse, mape, bias, mae_gain, mape_gain) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))),
      width = 160)

write_csv(metrics, "output/tables/direct_RF_metrics_overall.csv")

# By hour
metrics_by_hour <- eval_combined |>
  group_by(hour, model) |>
  summarise(mae  = mean(abs_error, na.rm = TRUE),
            mape = mean(abs_error / actual, na.rm = TRUE) * 100,
            .groups = "drop")

write_csv(metrics_by_hour, "output/tables/direct_RF_by_hour.csv")

# ============================================================
# Plots
# ============================================================

models_to_plot <- intersect(
  c(MODEL_NO_TEMP, MODEL_TEMP, OLS_NAME, "Official forecast"),
  unique(metrics_by_hour$model)
)

p_hour <- metrics_by_hour |>
  filter(model %in% models_to_plot) |>
  ggplot(aes(hour, mae, colour = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  labs(title = paste0("MAE by hour — Direct RF vs FPCA-OLS (",
                      TEST_START, " – ", TEST_END, ")"),
       x = "Hour", y = "MAE", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

ggsave("output/figs/direct_RF_mae_by_hour.png", p_hour,
       width = 9, height = 5, dpi = 150)
cat("Saved: output/figs/direct_RF_mae_by_hour.png\n")

# Temperature gain by hour (if both variants ran)
if (!is.null(pred_temp)) {
  no_temp_h  <- metrics_by_hour |> filter(model == MODEL_NO_TEMP) |>
    dplyr::select(hour, mae_no_temp = mae)
  with_temp_h <- metrics_by_hour |> filter(model == MODEL_TEMP) |>
    dplyr::select(hour, mae_temp = mae)

  temp_gain_h <- inner_join(no_temp_h, with_temp_h, by = "hour") |>
    mutate(temp_gain = 100 * (mae_no_temp - mae_temp) / mae_no_temp)

  write_csv(temp_gain_h, "output/tables/direct_RF_temperature_gain_by_hour.csv")

  p_temp <- ggplot(temp_gain_h, aes(x = hour, y = temp_gain)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_col(aes(fill = temp_gain > 0), show.legend = FALSE, alpha = 0.8) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
    scale_x_continuous(breaks = seq(0, 23, 3)) +
    labs(title = "Temperature gain in Direct RF by hour",
         subtitle = "Blue = temperature helps, Red = temperature hurts",
         x = "Hour", y = "MAE gain (%)") +
    theme_bw(base_size = 11)

  ggsave("output/figs/direct_RF_temperature_gain_by_hour.png", p_temp,
         width = 8, height = 4, dpi = 150)
  cat("Saved: output/figs/direct_RF_temperature_gain_by_hour.png\n")
}

# ============================================================
# Save eval files for DM test
# ============================================================

write_csv(eval_df |> filter(model == MODEL_NO_TEMP),
          "output/tables/direct_RF_no_temp_eval_long.csv")

if (!is.null(pred_temp)) {
  write_csv(eval_df |> filter(model == MODEL_TEMP),
            "output/tables/direct_RF_temp_eval_long.csv")
}

cat("\nDone.\n")
