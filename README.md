# turkey-electricity-fda
Functional data analysis of Turkey electricity load curves and day-ahead forecast errors.

## Current thesis framing

The advisor-facing baseline is the parsimonious FPCA VAR(1,7) model with
calendar, holiday, load-shape, and derivative features. It keeps `K = 2` load
FPCA components because the first two components explain more than 95% of load
curve variation.

Use the residual second step as an empirical extension/robustness exercise, not
as the main thesis claim, unless its methodological justification is expanded.
When running that extension, outputs are now written with an explicit `K` suffix
so `K = 2` and `K = 4` results are not mixed.

## Pipeline

Run scripts from the repository root:

```sh
Rscript R/01_clean_data.R
Rscript R/02_construct_curves.R
Rscript R/03_descriptives.R
Rscript R/06_forecast_fpca_VAR_1_7.R
Rscript R/09_forecast_fpca_VAR_1_7_ramp_calendar.R
Rscript R/10_forecast_fpca_VAR_1_7_holiday_regime.R
Rscript R/11_forecast_fpca_VAR_1_7_holiday_regime_derivative.R
FPCA_K=2 Rscript R/12_forecast_fpca_VAR_1_7_holiday_regime_derivative_residual_step.R
FPCA_K=4 Rscript R/12_forecast_fpca_VAR_1_7_holiday_regime_derivative_residual_step.R
Rscript R/13_diebold_mariano.R
Rscript R/14_acf_diagnostics.R
```

Key outputs are saved under `output/tables`, `output/figs`, and
`output/results`. The main report is `reports/advisor_progress_report.tex`.

## Next steps after advisor discussion

- Keep the `K = 2` derivative model as the conservative baseline.
- Report `K = 4` as a robustness check for richer intraday shape variation.
- Present the residual second step separately, with emphasis on rolling-origin
  validation and leakage prevention.
- Refresh the report tables only after a clean rerun of the final scripts.
- Add a short robustness table comparing `K = 2`, `K = 4`, and residual-step
  variants using the explicitly suffixed outputs.
