# Research plan

## Topic

Dynamic functional analysis of Turkey electricity load curves and day-ahead forecast errors.

## Core objects

Each day is treated as one functional observation over hour of day:

L_d(h), h = 0,...,23

where L_d(h) is realized hourly electricity consumption.

The official day-ahead forecast is:

Lhat_d(h)

The forecast error curve is:

e_d(h) = L_d(h) - Lhat_d(h)

## Main question

Are official day-ahead load forecast errors in Turkey systematically structured over the day, and can dynamic functional methods improve or correct these forecasts?

## Current empirical fact

Using 2023 hourly EPİAŞ data:

- Realized load curves are highly low-dimensional.
- FPC1 explains 92.7% of load curve variation.
- FPC1-FPC4 explain 99.4%.
- Forecast error curves are also structured.
- FPC1-FPC4 explain 95.3% of forecast-error variation.

## Next steps

1. Download multiple years.
2. Study seasonality and weekday effects.
3. Model forecast error curves dynamically.
4. Construct corrected forecasts:

Lhat_corrected_d(h) = Lhat_official_d(h) + ehat_d(h)

5. Compare official vs corrected forecasts using MAE/RMSE.