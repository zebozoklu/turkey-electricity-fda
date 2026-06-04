# 00b_download_weather.R
#
# Downloads daily temperature data for Turkey from the Open-Meteo
# historical weather API (free, no API key required).
#
# Strategy: fetch daily temperature for the four largest load centres
# (Istanbul, Ankara, İzmir, Bursa) and compute a population-weighted
# average. This approximates the demand-weighted national temperature
# that drives electricity consumption.
#
# Population weights (2023 estimates, Turkish Statistical Institute):
#   Istanbul  ~15.8M  → 0.44
#   Ankara    ~ 5.8M  → 0.16
#   İzmir     ~ 4.4M  → 0.12
#   Bursa     ~ 3.2M  → 0.09
#   Remainder ~36.0M weighted via Istanbul proxy → 0.19
#   (Using Istanbul for the remainder is conservative; sensitivity is low
#    because these cities are also on the western/Marmara grid.)
#
# Features produced:
#   temp_mean   weighted daily mean temperature (°C)
#   temp_max    weighted daily maximum
#   temp_min    weighted daily minimum
#   hdd         heating degree days  = max(BASE - temp_mean, 0)
#   cdd         cooling degree days  = max(temp_mean - BASE, 0)
#   roll_temp_7 7-day trailing mean of temp_mean
#   temp_sq     temp_mean² (captures U-shaped demand curve)
#
# Output: data/processed/weather_daily.rds  and  data/processed/weather_daily.csv

library(jsonlite)
library(dplyr)
library(lubridate)
library(readr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Config
# ============================================================

BASE_TEMP  <- 18     # standard HDD/CDD base for Turkey
START_DATE <- "2018-12-01"   # one month before data start for rolling features
END_DATE   <- "2025-12-31"

cities <- tribble(
  ~name,       ~lat,     ~lon,     ~weight,
  "Istanbul",  41.0082,  28.9784,  0.44,
  "Ankara",    39.9334,  32.8597,  0.16,
  "Izmir",     38.4192,  27.1287,  0.12,
  "Bursa",     40.1826,  29.0665,  0.09,
  "Adana",     37.0000,  35.3213,  0.05,
  "Konya",     37.8667,  32.4833,  0.05,
  "Antalya",   36.8841,  30.7056,  0.09
)

# ============================================================
# Download helper
# ============================================================

fetch_city_temp <- function(city_name, lat, lon) {
  cat("  Downloading:", city_name, "\n")
  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    "latitude=",  lat, "&longitude=", lon, "&",
    "start_date=", START_DATE, "&end_date=", END_DATE, "&",
    "daily=temperature_2m_mean,temperature_2m_max,temperature_2m_min&",
    "timezone=Europe%2FIstanbul"
  )
  Sys.sleep(0.5)   # be polite to the API
  raw <- tryCatch(fromJSON(url), error = function(e) {
    warning("Failed to fetch ", city_name, ": ", e$message); NULL
  })
  if (is.null(raw) || is.null(raw$daily)) return(NULL)
  tibble(
    city       = city_name,
    date       = as.Date(raw$daily$time),
    temp_mean  = raw$daily$temperature_2m_mean,
    temp_max   = raw$daily$temperature_2m_max,
    temp_min   = raw$daily$temperature_2m_min
  )
}

# ============================================================
# Fetch all cities
# ============================================================

cat("Fetching temperature data from Open-Meteo...\n")

city_data <- bind_rows(mapply(
  fetch_city_temp,
  cities$name, cities$lat, cities$lon,
  SIMPLIFY = FALSE
))

cat("Rows fetched:", nrow(city_data), "\n")
cat("Cities obtained:", paste(unique(city_data$city), collapse = ", "), "\n")

# ============================================================
# Weighted average
# ============================================================

# Merge weights
city_data <- city_data |>
  left_join(cities |> dplyr::select(name, weight), by = c("city" = "name"))

# Renormalise weights in case some cities failed
city_data <- city_data |>
  group_by(date) |>
  mutate(weight = weight / sum(weight, na.rm = TRUE)) |>
  ungroup()

weather <- city_data |>
  group_by(date) |>
  summarise(
    temp_mean = sum(temp_mean * weight, na.rm = TRUE),
    temp_max  = sum(temp_max  * weight, na.rm = TRUE),
    temp_min  = sum(temp_min  * weight, na.rm = TRUE),
    n_cities  = n(),
    .groups   = "drop"
  ) |>
  arrange(date)

cat("Date range:", as.character(min(weather$date)),
    "to", as.character(max(weather$date)), "\n")
cat("Missing days:", sum(is.na(weather$temp_mean)), "\n")

# ============================================================
# Derived features
# ============================================================

weather <- weather |>
  mutate(
    hdd         = pmax(BASE_TEMP - temp_mean, 0),
    cdd         = pmax(temp_mean - BASE_TEMP, 0),
    temp_sq     = temp_mean^2,
    roll_temp_7 = as.numeric(
      stats::filter(temp_mean, rep(1/7, 7), sides = 1)
    ),
    roll_temp_28 = as.numeric(
      stats::filter(temp_mean, rep(1/28, 28), sides = 1)
    )
  )

# ============================================================
# Save
# ============================================================

saveRDS(weather, "data/processed/weather_daily.rds")
write_csv(weather, "data/processed/weather_daily.csv")

cat("\nSummary of weather features:\n")
print(summary(weather |> dplyr::select(temp_mean, hdd, cdd, roll_temp_7)))
cat("\nSaved to data/processed/weather_daily.rds and .csv\n")
