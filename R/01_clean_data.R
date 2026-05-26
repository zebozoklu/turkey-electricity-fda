# 01_clean_data.R

library(readr)
library(dplyr)
library(janitor)
library(lubridate)
library(stringr)

# paths
forecast_path <- "data/raw/Load_Forecast_Plan-01012023-31122023.csv"
actual_path   <- "data/raw/Real_Time_Consumption-01012023-31122023.csv"

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
      )
    ) |>
    rename(!!value_name := value)
}

load_forecast <- read_epias_hourly(forecast_path, "load_forecast_mwh")
load_actual   <- read_epias_hourly(actual_path, "load_actual_mwh")

load_hourly <- load_actual |>
  left_join(load_forecast, by = c("date", "hour", "hour_id")) |>
  mutate(
    forecast_error_mwh = load_actual_mwh - load_forecast_mwh,
    abs_error_mwh = abs(forecast_error_mwh),
    year = year(date),
    month = month(date),
    weekday = wday(date, label = TRUE, week_start = 1)
  ) |>
  arrange(date, hour_id)

# basic checks
print(dim(load_hourly))
print(head(load_hourly))
print(summary(load_hourly))

# should be 8760 for 2023
stopifnot(nrow(load_hourly) == 365 * 24)
stopifnot(all(load_hourly$hour_id %in% 0:23))

saveRDS(load_hourly, "data/processed/load_hourly_2023.rds")



