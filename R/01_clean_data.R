# 01_clean_data.R

library(readr)
library(dplyr)
library(janitor)
library(lubridate)
library(stringr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

read_epias_hourly <- function(path, value_name) {
  read_delim(
    file = path,
    delim = ";",
    locale = locale(grouping_mark = ",", decimal_mark = "."),
    show_col_types = FALSE
  ) |>
    clean_names() |>
    rename(
      date = 1,
      hour = 2,
      value = 3
    ) |>
    mutate(
      date = dmy(date),
      hour = str_sub(hour, 1, 5),
      hour_id = as.integer(str_sub(hour, 1, 2)),
      value = parse_number(
        as.character(value),
        locale = locale(grouping_mark = ",", decimal_mark = ".")
      ),
      source_file = basename(path)
    ) |>
    rename(!!value_name := value)
}

actual_files <- list.files(
  "data/raw",
  pattern = "Real_Time_Consumption.*\\.csv$",
  full.names = TRUE
)

forecast_files <- list.files(
  "data/raw",
  pattern = "Load_Forecast_Plan.*\\.csv$",
  full.names = TRUE
)

print(actual_files)
print(forecast_files)

load_actual <- bind_rows(
  lapply(actual_files, read_epias_hourly, value_name = "load_actual_mwh")
) |>
  distinct(date, hour, hour_id, .keep_all = TRUE)

load_forecast <- bind_rows(
  lapply(forecast_files, read_epias_hourly, value_name = "load_forecast_mwh")
) |>
  distinct(date, hour, hour_id, .keep_all = TRUE)

load_hourly <- load_actual |>
  dplyr::select(date, hour, hour_id, load_actual_mwh) |>
  left_join(
    load_forecast |> dplyr::select(date, hour, hour_id, load_forecast_mwh),
    by = c("date", "hour", "hour_id")
  ) |>
  mutate(
    forecast_error_mwh = load_actual_mwh - load_forecast_mwh,
    abs_error_mwh = abs(forecast_error_mwh),
    year = year(date),
    month = month(date),
    weekday = wday(date, label = TRUE, week_start = 1)
  ) |>
  arrange(date, hour_id)

# checks
check_days <- load_hourly |>
  count(date) |>
  summarise(
    min_hours = min(n),
    max_hours = max(n),
    bad_days = sum(n != 24)
  )

print(check_days)

summary_table <- load_hourly |>
  group_by(year) |>
  summarise(
    n = n(),
    days = n_distinct(date),
    start = min(date),
    end = max(date),
    missing_actual = sum(is.na(load_actual_mwh)),
    missing_forecast = sum(is.na(load_forecast_mwh)),
    mae = mean(abs_error_mwh, na.rm = TRUE),
    rmse = sqrt(mean(forecast_error_mwh^2, na.rm = TRUE)),
    bias = mean(forecast_error_mwh, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

stopifnot(check_days$bad_days == 0)
stopifnot(all(load_hourly$hour_id %in% 0:23))

saveRDS(load_hourly, "data/processed/load_hourly.rds")
write.csv(summary_table, "output/tables/load_summary_by_year.csv", row.names = FALSE)


