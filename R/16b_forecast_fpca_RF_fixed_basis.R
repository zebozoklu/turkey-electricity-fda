# 16b_forecast_fpca_RF_fixed_basis.R
#
# RF benchmark with a FIXED FPCA basis (estimated once on 2019-2022 training
# data). The rolling loop re-fits only the RF score regression each day —
# no FPCA refit inside the loop.
#
# Methodological justification: script 15 showed that pairwise inner products
# between eigenfunctions fitted on 2019-2021, 2020-2022, and 2022-2024
# sub-periods all exceed 0.993. The basis is empirically stable, so a fixed
# basis is appropriate and the fixed-basis scores are essentially the same as
# the rolling-basis scores used in script 11.
#
# This is ~10x faster than script 16 (no FPCA refit per day).

library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)
library(lubridate)
library(readr)
library(ranger)

dir.create("output/figs",    recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Settings  (match script 11 / 16)
# ============================================================

K          <- 2
K_DERIV    <- 3
NBASIS     <- 12
NORDER     <- 4
TEST_START <- as.Date("2023-01-01")
MIN_TRAIN_DAYS <- 365
MIN_LAG_DAYS   <- 28
RF_TREES   <- 500
RF_SEED    <- 42

MODEL_NAME <- "FPCA RF (fixed basis)"
OLS_NAME   <- "FPCA VAR(1,7) + holiday/load/derivative regime"

# ============================================================
# Load data
# ============================================================

curves <- readRDS("data/processed/load_curves.rds")
dates  <- as.Date(curves$dates)
hours  <- curves$hours
Y      <- curves$actual
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

hcol <- function(h) { out <- match(h, hours); if (is.na(out)) stop("Hour ", h); out }
H5  <- hcol(5);  H6  <- hcol(6);  H7  <- hcol(7);  H9  <- hcol(9)
H17 <- hcol(17); H20 <- hcol(20); H23 <- hcol(23)

# ============================================================
# Calendar / holiday helpers  (identical to scripts 11 / 16)
# ============================================================

date_seq <- function(start, end) seq(as.Date(start), as.Date(end), by = "day")

make_holiday_lookup <- function(years) {
  fixed_holidays <- bind_rows(lapply(years, function(y) {
    tibble(
      date = as.Date(c(
        paste0(y,"-01-01"), paste0(y,"-04-23"), paste0(y,"-05-01"),
        paste0(y,"-05-19"), paste0(y,"-07-15"), paste0(y,"-08-30"),
        paste0(y,"-10-29")
      )),
      holiday_type = "fixed"
    )
  }))
  religious_ranges <- tribble(
    ~start,~end,~holiday_type,
    "2019-06-04","2019-06-06","ramadan_bayram","2020-05-24","2020-05-26","ramadan_bayram",
    "2021-05-13","2021-05-15","ramadan_bayram","2022-05-02","2022-05-04","ramadan_bayram",
    "2023-04-21","2023-04-23","ramadan_bayram","2024-04-10","2024-04-12","ramadan_bayram",
    "2025-03-30","2025-04-01","ramadan_bayram",
    "2019-08-11","2019-08-14","kurban_bayram","2020-07-31","2020-08-03","kurban_bayram",
    "2021-07-20","2021-07-23","kurban_bayram","2022-07-09","2022-07-12","kurban_bayram",
    "2023-06-28","2023-07-01","kurban_bayram","2024-06-16","2024-06-19","kurban_bayram",
    "2025-06-06","2025-06-09","kurban_bayram"
  ) |> mutate(start=as.Date(start), end=as.Date(end))
  religious_holidays <- bind_rows(lapply(seq_len(nrow(religious_ranges)), function(i)
    tibble(date=date_seq(religious_ranges$start[i], religious_ranges$end[i]),
           holiday_type=religious_ranges$holiday_type[i])))
  bind_rows(fixed_holidays, religious_holidays) |>
    filter(year(date) %in% years) |> distinct(date, .keep_all=TRUE)
}

make_ramadan_lookup <- function(years) {
  rr <- tribble(~start,~end,
    "2019-05-06","2019-06-03","2020-04-24","2020-05-23","2021-04-13","2021-05-12",
    "2022-04-02","2022-05-01","2023-03-23","2023-04-20","2024-03-11","2024-04-09",
    "2025-03-01","2025-03-29") |> mutate(start=as.Date(start), end=as.Date(end))
  bind_rows(lapply(seq_len(nrow(rr)), function(i)
    tibble(date=date_seq(rr$start[i], rr$end[i])))) |>
    filter(year(date) %in% years) |> distinct(date)
}

holiday_lookup <- make_holiday_lookup(seq(min(year(dates))-1, max(year(dates))+1))
ramadan_lookup <- make_ramadan_lookup(seq(min(year(dates))-1, max(year(dates))+1))

make_calendar_features <- function(date_vec) {
  doy <- yday(date_vec); dow <- wday(date_vec, week_start=1)
  mn  <- month(date_vec); hd  <- holiday_lookup$date
  tibble(date=date_vec) |>
    mutate(
      dow=dow, weekend=as.integer(dow>=6),
      monday=as.integer(dow==1), friday=as.integer(dow==5), month=mn,
      sin_week=sin(2*pi*dow/7), cos_week=cos(2*pi*dow/7),
      sin_year=sin(2*pi*doy/365.25), cos_year=cos(2*pi*doy/365.25),
      summer_school_break=as.integer(
        (mn==6 & mday(date)>=15)|mn%in%c(7,8)|(mn==9 & mday(date)<=15)),
      is_holiday=as.integer(date %in% hd),
      is_religious_holiday=as.integer(date %in% holiday_lookup$date[
        holiday_lookup$holiday_type %in% c("ramadan_bayram","kurban_bayram")]),
      is_ramadan=as.integer(date %in% ramadan_lookup$date),
      holiday_eve=as.integer((date+days(1)) %in% hd),
      post_holiday=as.integer((date-days(1)) %in% hd),
      near_holiday=as.integer(date %in% c(hd-days(2),hd-days(1),hd,hd+days(1),hd+days(2))),
      bridge_day=as.integer(dow%in%c(1,5) & !(date%in%hd) &
        ((date-days(1))%in%hd|(date+days(1))%in%hd))
    ) |> dplyr::select(-date)
}

safe_window_mean <- function(Ym,s,e) mean(rowMeans(Ym[s:e,,drop=FALSE],na.rm=TRUE),na.rm=TRUE)
safe_window_ramp <- function(Ym,s,e) mean(Ym[s:e,H9]-Ym[s:e,H6],na.rm=TRUE)

make_fd_object <- function(Y_mat, argvals, nbasis, norder) {
  basis <- create.bspline.basis(rangeval=range(argvals), nbasis=nbasis, norder=norder)
  Data2fd(argvals=argvals, y=t(Y_mat), basisobj=basis)
}

make_load_regime_features <- function(Y_mat, D_mat, d_index) {
  bind_rows(lapply(d_index, function(d) {
    lag1=Y_mat[d-1,]; lag7=Y_mat[d-7,]
    dlag1=D_mat[d-1,]; dlag7=D_mat[d-7,]
    tibble(
      ramp_lag1=as.numeric(Y_mat[d-1,H9]-Y_mat[d-1,H6]),
      ramp_lag7=as.numeric(Y_mat[d-7,H9]-Y_mat[d-7,H6]),
      early_lag1=mean(Y_mat[d-1,c(H5,H6,H7)],na.rm=TRUE),
      early_ramp_lag1=as.numeric(Y_mat[d-1,H7]-Y_mat[d-1,H5]),
      early_ramp_lag7=as.numeric(Y_mat[d-7,H7]-Y_mat[d-7,H5]),
      evening_ramp_lag1=as.numeric(Y_mat[d-1,H20]-Y_mat[d-1,H17]),
      evening_ramp_lag7=as.numeric(Y_mat[d-7,H20]-Y_mat[d-7,H17]),
      day_end_drop_lag1=as.numeric(Y_mat[d-1,H23]-Y_mat[d-1,H20]),
      day_end_drop_lag7=as.numeric(Y_mat[d-7,H23]-Y_mat[d-7,H20]),
      mean_lag1=mean(lag1,na.rm=TRUE), mean_lag7=mean(lag7,na.rm=TRUE),
      peak_lag1=max(lag1,na.rm=TRUE),  peak_lag7=max(lag7,na.rm=TRUE),
      min_lag1=min(lag1,na.rm=TRUE),   min_lag7=min(lag7,na.rm=TRUE),
      daily_range_lag1=max(lag1,na.rm=TRUE)-min(lag1,na.rm=TRUE),
      daily_range_lag7=max(lag7,na.rm=TRUE)-min(lag7,na.rm=TRUE),
      roll_mean_7=safe_window_mean(Y_mat,d-7,d-1),
      roll_mean_28=safe_window_mean(Y_mat,d-28,d-1),
      roll_ramp_7=safe_window_ramp(Y_mat,d-7,d-1),
      roll_ramp_28=safe_window_ramp(Y_mat,d-28,d-1),
      deriv_sd_lag1=sd(dlag1,na.rm=TRUE),   deriv_sd_lag7=sd(dlag7,na.rm=TRUE),
      deriv_max_lag1=max(dlag1,na.rm=TRUE),  deriv_max_lag7=max(dlag7,na.rm=TRUE),
      deriv_min_lag1=min(dlag1,na.rm=TRUE),  deriv_min_lag7=min(dlag7,na.rm=TRUE),
      deriv_range_lag1=max(dlag1,na.rm=TRUE)-min(dlag1,na.rm=TRUE),
      deriv_range_lag7=max(dlag7,na.rm=TRUE)-min(dlag7,na.rm=TRUE)
    )
  }))
}

standardize_train_new <- function(Z_train, Z_new) {
  mu  <- sapply(Z_train, mean, na.rm=TRUE)
  sdv <- sapply(Z_train, sd,   na.rm=TRUE); sdv[is.na(sdv)|sdv==0] <- 1
  list(Z_train=sweep(sweep(as.matrix(Z_train),2,mu,"-"),2,sdv,"/"),
       Z_new  =sweep(sweep(as.matrix(Z_new),  2,mu,"-"),2,sdv,"/"))
}

# ============================================================
# Step 1: Fit FIXED FPCA on training data only
# ============================================================

cat("Fitting fixed FPCA on training data (", as.character(TEST_START), " cutoff)...\n")

train_mask <- dates < TEST_START
Y_train    <- Y[train_mask, ]
fd_train   <- make_fd_object(Y_train, hours, NBASIS, NORDER)

pca_fixed        <- pca.fd(fd_train, nharm = K)
mu_fixed         <- as.vector(eval.fd(hours, pca_fixed$meanfd))
phi_fixed        <- eval.fd(hours, pca_fixed$harmonics)   # 24 x K

fd_deriv_train   <- deriv.fd(fd_train, Lfdobj = 1)
deriv_pca_fixed  <- pca.fd(fd_deriv_train, nharm = K_DERIV)
dmu_fixed        <- as.vector(eval.fd(hours, deriv_pca_fixed$meanfd))
dphi_fixed       <- eval.fd(hours, deriv_pca_fixed$harmonics)  # 24 x K_DERIV

cat("  FPC variance (load):    ", round(pca_fixed$varprop[1:K]*100, 2), "\n")
cat("  FPC variance (deriv):   ", round(deriv_pca_fixed$varprop[1:K_DERIV]*100, 2), "\n")

# ============================================================
# Step 2: Project ALL days onto fixed basis
# ============================================================

cat("Projecting all", n_days, "days onto fixed basis...\n")

# Derivative matrix for all days
fd_all       <- make_fd_object(Y, hours, NBASIS, NORDER)
fd_deriv_all <- deriv.fd(fd_all, Lfdobj = 1)
D_all        <- t(eval.fd(hours, fd_deriv_all))   # n x 24
colnames(D_all) <- as.character(hours)

# Scores: discrete inner product approximation
# xi_dk = integral (Y_d(h) - mu(h)) phi_k(h) dh  ≈  sum_h (Y_d(h) - mu(h)) phi_k(h)
scores_all       <- sweep(Y, 2, mu_fixed,  "-") %*% phi_fixed   # n x K
deriv_scores_all <- sweep(D_all, 2, dmu_fixed, "-") %*% dphi_fixed  # n x K_DERIV

colnames(scores_all)       <- paste0("score",       1:K)
colnames(deriv_scores_all) <- paste0("deriv_score", 1:K_DERIV)

cat("  Score range (S1):", round(range(scores_all[,1])), "\n")
cat("  Score range (S2):", round(range(scores_all[,2])), "\n")

# ============================================================
# Step 3: Precompute external features
# ============================================================

cat("Precomputing external features...\n")

feature_indices  <- (MIN_LAG_DAYS + 1):n_days
all_feature_rows <- bind_cols(
  make_load_regime_features(Y, D_all, feature_indices),
  make_calendar_features(dates[feature_indices])
)

get_feature_rows <- function(idx_vec) {
  all_feature_rows[match(idx_vec, feature_indices), , drop=FALSE]
}

# ============================================================
# Step 4: Test indices
# ============================================================

test_idx <- which(dates >= TEST_START)
test_idx <- test_idx[test_idx > max(MIN_LAG_DAYS, MIN_TRAIN_DAYS)]

cat("\nForecast days:", length(test_idx),
    "| First:", as.character(min(dates[test_idx])),
    "| Last:", as.character(max(dates[test_idx])), "\n\n")

# ============================================================
# Step 5: Rolling RF loop  (no FPCA refit inside)
# ============================================================

pred_rf     <- matrix(NA_real_, nrow=length(test_idx), ncol=length(hours),
                      dimnames=list(as.character(dates[test_idx]), as.character(hours)))
score_preds <- matrix(NA_real_, nrow=length(test_idx), ncol=K)

cat("Starting rolling-origin RF forecast (fixed basis)...\n")
t0 <- proc.time()

for (j in seq_along(test_idx)) {
  i <- test_idx[j]

  # Training indices: days with enough lag history, strictly before i
  d_index <- intersect(feature_indices, (MIN_LAG_DAYS + 1):(i - 1))
  if (length(d_index) < 50) next   # skip if too little training data

  # Response: scores for training days
  Y_score <- scores_all[d_index, , drop=FALSE]

  # Predictors
  X_lag1  <- scores_all[d_index - 1, , drop=FALSE]
  X_lag7  <- scores_all[d_index - 7, , drop=FALSE]
  X_dlag1 <- deriv_scores_all[d_index - 1, , drop=FALSE]
  X_dlag7 <- deriv_scores_all[d_index - 7, , drop=FALSE]

  Z_train_raw <- get_feature_rows(d_index)
  Z_new_raw   <- get_feature_rows(i)
  Z_sc        <- standardize_train_new(Z_train_raw, Z_new_raw)

  X_rf <- cbind(X_lag1, X_lag7, X_dlag1, X_dlag7, Z_sc$Z_train)
  colnames(X_rf) <- c(
    paste0("lag1_s",  1:K),
    paste0("lag7_s",  1:K),
    paste0("lag1_ds", 1:K_DERIV),
    paste0("lag7_ds", 1:K_DERIV),
    colnames(Z_train_raw)
  )

  x_new <- c(scores_all[i-1,], scores_all[i-7,],
              deriv_scores_all[i-1,], deriv_scores_all[i-7,],
              as.numeric(Z_sc$Z_new))
  x_new_df <- setNames(as.data.frame(t(x_new)), colnames(X_rf))

  # Fit RF per score component
  score_hat <- numeric(K)
  for (k in seq_len(K)) {
    df_tr  <- as.data.frame(cbind(y = Y_score[, k], X_rf))
    rf_fit <- ranger(y ~ ., data=df_tr, num.trees=RF_TREES,
                     seed=RF_SEED, num.threads=1, verbose=FALSE)
    score_hat[k] <- predict(rf_fit, data=x_new_df)$predictions
  }

  # Reconstruct curve from fixed basis
  pred_rf[j, ]     <- mu_fixed + phi_fixed %*% score_hat
  score_preds[j, ] <- score_hat

  if (j %% 50 == 0) {
    el  <- round((proc.time()-t0)["elapsed"])
    rem <- round(el/j*(length(test_idx)-j))
    message("  ", j, "/", length(test_idx),
            "  date: ", dates[i],
            "  elapsed: ", el, "s",
            "  est remaining: ", rem, "s")
  }
}

cat("Done. Total time:", round((proc.time()-t0)["elapsed"]), "s\n\n")

# ============================================================
# Step 6: Evaluation
# ============================================================

actual_test  <- Y[test_idx, , drop=FALSE]
pred_naive   <- Y[test_idx - 7, , drop=FALSE]
if (has_official) pred_official <- Y_official[test_idx, , drop=FALSE]

ols_eval <- read_csv(
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  show_col_types=FALSE, col_types=cols(date=col_character())
) |>
  mutate(date=as.Date(date), hour=as.integer(hour)) |>
  filter(model == OLS_NAME)

make_eval_df <- function(actual, pred, model_name) {
  tibble(
    date     = rep(dates[test_idx], each=length(hours)),
    hour     = rep(hours, times=length(test_idx)),
    model    = model_name,
    actual   = as.vector(t(actual)),
    forecast = as.vector(t(pred))
  ) |> mutate(
    year=year(date), month=month(date), ym=format(date,"%Y-%m"),
    error=actual-forecast, abs_error=abs(error), sq_error=error^2
  )
}

eval_df <- bind_rows(
  make_eval_df(actual_test, pred_rf,    MODEL_NAME),
  make_eval_df(actual_test, pred_naive, "Seasonal naive: Y[d-7]")
)
if (has_official)
  eval_df <- bind_rows(eval_df, make_eval_df(actual_test, pred_official, "Official forecast"))

eval_combined <- bind_rows(eval_df, ols_eval |> select(names(eval_df)))

# Metrics
metrics <- eval_combined |>
  filter(actual > 0) |>
  group_by(model) |>
  summarise(n=n(), mae=mean(abs_error,na.rm=TRUE),
            rmse=sqrt(mean(sq_error,na.rm=TRUE)),
            mape=mean(abs_error/actual,na.rm=TRUE)*100,
            bias=mean(error,na.rm=TRUE), .groups="drop")

off <- metrics |> filter(model=="Official forecast")
if (nrow(off)==1) {
  metrics <- metrics |>
    mutate(mae_gain  = 100*(off$mae  - mae)  / off$mae,
           rmse_gain = 100*(off$rmse - rmse) / off$rmse,
           mape_gain = 100*(off$mape - mape) / off$mape)
}

cat("=== Overall metrics ===\n")
print(metrics |> select(model,mae,rmse,mape,bias,mae_gain,mape_gain) |>
      mutate(across(where(is.numeric),\(x) round(x,2))), width=130)

write_csv(metrics, "output/tables/fpca_RF_fixed_metrics_overall.csv")

# By-hour comparison
metrics_by_hour <- eval_combined |>
  group_by(hour, model) |>
  summarise(mae=mean(abs_error,na.rm=TRUE),
            mape=mean(abs_error/actual,na.rm=TRUE)*100, .groups="drop")

write_csv(metrics_by_hour, "output/tables/fpca_RF_fixed_by_hour.csv")

# RF vs OLS gain by hour
rf_h  <- metrics_by_hour |> filter(model==MODEL_NAME) |> select(hour, rf_mae=mae, rf_mape=mape)
ols_h <- metrics_by_hour |> filter(model==OLS_NAME)   |> select(hour, ols_mae=mae, ols_mape=mape)
gain_h <- inner_join(rf_h, ols_h, by="hour") |>
  mutate(mae_gain  = 100*(ols_mae  - rf_mae)  / ols_mae,
         mape_gain = 100*(ols_mape - rf_mape) / ols_mape)

write_csv(gain_h, "output/tables/fpca_RF_fixed_gain_vs_OLS_by_hour.csv")

# ============================================================
# Plots
# ============================================================

p_hour <- metrics_by_hour |>
  filter(model %in% c(MODEL_NAME, OLS_NAME, "Official forecast")) |>
  ggplot(aes(hour, mae, colour=model, linetype=model)) +
  geom_line(linewidth=0.9) +
  scale_x_continuous(breaks=seq(0,23,3)) +
  labs(title="MAE by hour: RF (fixed basis) vs OLS vs Official",
       x="Hour", y="MAE", colour=NULL, linetype=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=8))

ggsave("output/figs/fpca_RF_fixed_mae_by_hour.png", p_hour, width=9, height=5, dpi=150)

p_gain <- ggplot(gain_h, aes(x=hour, y=mae_gain)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_col(aes(fill=mae_gain > 0), show.legend=FALSE, alpha=0.8) +
  scale_fill_manual(values=c("TRUE"="steelblue","FALSE"="tomato")) +
  scale_x_continuous(breaks=seq(0,23,3)) +
  labs(title="RF (fixed basis) gain over OLS by hour",
       subtitle="Blue = RF better, Red = OLS better",
       x="Hour", y="MAE gain (%)") +
  theme_bw(base_size=11)

ggsave("output/figs/fpca_RF_fixed_gain_vs_OLS_by_hour.png", p_gain, width=8, height=4, dpi=150)

cat("Saved figures.\n")

# Save eval for DM test
write_csv(eval_df |> filter(model==MODEL_NAME),
          "output/tables/fpca_RF_fixed_eval_long.csv")

cat("\nDone.\n")
