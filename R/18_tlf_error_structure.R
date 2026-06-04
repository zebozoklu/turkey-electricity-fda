# 18_tlf_error_structure.R
#
# Does TLF leave behind curve-shape structure that FDA exploits?
#
# Two analyses:
#   (1) Project TLF error curves onto the load FPCA eigenfunctions.
#       Systematic non-zero mean scores = TLF misses that shape mode.
#       Compare with FDA model error projections.
#   (2) Cross-hour correlation matrix of TLF errors (and FDA errors).
#       Low effective rank = errors co-move as curves, not independent scalars.

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(fda)

dir.create("output/figs",   recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

NBASIS     <- 12
NORDER     <- 4
K          <- 2
TEST_START <- as.Date("2023-01-01")
MAIN_MODEL <- "FPCA VAR(1,7) + holiday/load/derivative regime"
BENCHMARK  <- "Official forecast"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

curves    <- readRDS("data/processed/load_curves.rds")
dates_all <- as.Date(curves$dates)
hours     <- curves$hours
Y_all     <- curves$actual
ok        <- complete.cases(Y_all)
dates_all <- dates_all[ok]
Y_all     <- Y_all[ok, ]

eval_df <- read_csv(
  "output/tables/fpca_VAR_1_7_holiday_regime_derivative_eval_long.csv",
  show_col_types = FALSE,
  col_types = cols(date = col_character())
) |>
  mutate(date = as.Date(as.character(date)),
         hour = as.integer(hour))

# ------------------------------------------------------------
# Fit FPCA on training data (same basis as main model)
# ------------------------------------------------------------

train_mask  <- dates_all < TEST_START
Y_train     <- Y_all[train_mask, ]
dates_train <- dates_all[train_mask]

basis_obj <- fda::create.bspline.basis(rangeval = range(hours),
                                        nbasis   = NBASIS,
                                        norder   = NORDER)
fd_train  <- fda::Data2fd(argvals = hours, y = t(Y_train), basisobj = basis_obj)
pca_train <- fda::pca.fd(fd_train, nharm = K)

cat(sprintf("Variance explained by first %d FPCs: %.1f%%\n", K,
            sum(pca_train$varprop[1:K]) * 100))

# Eigenfunctions evaluated on 0:23
ef_mat <- fda::eval.fd(hours, pca_train$harmonics)  # 24 x K

# ------------------------------------------------------------
# Build error curve matrices for test period
# ------------------------------------------------------------

to_error_matrix <- function(df, model_label) {
  df |>
    dplyr::filter(model == model_label, date >= TEST_START) |>
    dplyr::select(date, hour, error) |>
    dplyr::arrange(date, hour) |>
    pivot_wider(names_from = hour, values_from = error) |>
    dplyr::arrange(date)
}

tlf_wide <- to_error_matrix(eval_df, BENCHMARK)
fda_wide <- to_error_matrix(eval_df, MAIN_MODEL)

shared_dates <- as.Date(intersect(as.character(tlf_wide$date), as.character(fda_wide$date)))
tlf_mat <- as.matrix(tlf_wide |> dplyr::filter(date %in% shared_dates) |> dplyr::select(-date))
fda_mat <- as.matrix(fda_wide |> dplyr::filter(date %in% shared_dates) |> dplyr::select(-date))
n_test  <- nrow(tlf_mat)
cat(sprintf("Test days used: %d\n", n_test))

# ------------------------------------------------------------
# (1a) Project error curves onto FPCA eigenfunctions
# ------------------------------------------------------------
# Score for day d on component k = integral of error_d(t) * ef_k(t) dt
# Approximated as row dot product (hours equally spaced, width = 1h).

tlf_scores <- tlf_mat %*% ef_mat   # n_test x K
fda_scores <- fda_mat %*% ef_mat   # n_test x K

scores_df <- bind_rows(
  as_tibble(tlf_scores, .name_repair = "minimal") |>
    setNames(paste0("PC", 1:K)) |>
    mutate(model = BENCHMARK,  day = shared_dates) |>
    pivot_longer(starts_with("PC"), names_to = "component", values_to = "score"),
  as_tibble(fda_scores, .name_repair = "minimal") |>
    setNames(paste0("PC", 1:K)) |>
    mutate(model = "FDA model", day = shared_dates) |>
    pivot_longer(starts_with("PC"), names_to = "component", values_to = "score")
)

# t-test: is the mean score significantly different from zero?
score_tests <- scores_df |>
  group_by(model, component) |>
  summarise(
    mean_score = mean(score),
    sd_score   = sd(score),
    t_stat     = mean(score) / (sd(score) / sqrt(n())),
    p_value    = 2 * pt(-abs(mean(score) / (sd(score) / sqrt(n()))),
                         df = n() - 1),
    .groups    = "drop"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))

cat("\n=== Mean FPCA score of error curves (H0: mean = 0) ===\n")
print(score_tests)
write_csv(score_tests, "output/tables/tlf_error_fpca_score_tests.csv")

# Boxplot of scores
p_scores <- ggplot(scores_df, aes(x = model, y = score, fill = model)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "black") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.5) +
  facet_wrap(~component, scales = "free_y") +
  scale_fill_manual(values = c("FDA model" = "#2166ac",
                                "Official forecast" = "#d6604d"),
                    guide = "none") +
  labs(
    title    = "Distribution of error-curve projections onto load FPCA eigenfunctions",
    subtitle = "Non-zero mean = systematic bias along that shape mode",
    x        = NULL,
    y        = "Projection score (MWh·h)"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))

ggsave("output/figs/tlf_error_fpca_scores.png", p_scores,
       width = 7, height = 4, dpi = 150)
cat("Saved: output/figs/tlf_error_fpca_scores.png\n")

# Mean error curves overlaid on eigenfunctions (visual alignment check)
mean_tlf_err <- colMeans(tlf_mat)
mean_fda_err <- colMeans(fda_mat)

ef_df <- as_tibble(ef_mat, .name_repair = "minimal") |>
  setNames(paste0("PC", 1:K)) |>
  mutate(hour = hours) |>
  pivot_longer(starts_with("PC"), names_to = "component", values_to = "ef_value")

mean_err_df <- tibble(
  hour    = hours,
  `Official forecast` = mean_tlf_err,
  `FDA model`         = mean_fda_err
) |>
  pivot_longer(-hour, names_to = "model", values_to = "mean_error")

p_mean_err <- ggplot(mean_err_df, aes(x = hour, y = mean_error, colour = model)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("FDA model" = "#2166ac",
                                  "Official forecast" = "#d6604d"),
                      name = NULL) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  labs(
    title    = "Mean forecast error curve by hour of day (test period)",
    subtitle = "Shape of the mean error indicates systematic bias across the day",
    x        = "Hour of day",
    y        = "Mean error (MWh)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave("output/figs/tlf_error_mean_curve.png", p_mean_err,
       width = 7, height = 4, dpi = 150)
cat("Saved: output/figs/tlf_error_mean_curve.png\n")

# ------------------------------------------------------------
# (1b) How much variance of TLF error curves do load FPCs explain?
# ------------------------------------------------------------

var_explained <- function(err_mat, ef) {
  proj     <- err_mat %*% ef          # scores
  reconstr <- proj %*% t(ef)          # reconstruction
  resid    <- err_mat - reconstr
  1 - sum(resid^2) / sum(err_mat^2)
}

ve_tlf <- var_explained(tlf_mat, ef_mat)
ve_fda <- var_explained(fda_mat, ef_mat)
cat(sprintf("\nVariance of TLF error curves explained by load PCs 1-%d: %.1f%%\n", K, ve_tlf * 100))
cat(sprintf("Variance of FDA error curves explained by load PCs 1-%d: %.1f%%\n", K, ve_fda * 100))

ve_df <- tibble(
  model          = c(BENCHMARK, "FDA model"),
  var_explained  = round(c(ve_tlf, ve_fda) * 100, 2)
)
write_csv(ve_df, "output/tables/tlf_error_fpca_var_explained.csv")

# ------------------------------------------------------------
# (2) Cross-hour correlation matrix and effective rank
# ------------------------------------------------------------

cor_tlf <- cor(tlf_mat)
cor_fda <- cor(fda_mat)

# Effective rank: entropy-based (sum of singular values normalized)
eff_rank <- function(mat) {
  sv  <- svd(cov(mat))$d
  sv  <- sv[sv > 0]
  p   <- sv / sum(sv)
  exp(-sum(p * log(p)))
}

er_tlf <- eff_rank(tlf_mat)
er_fda <- eff_rank(fda_mat)
cat(sprintf("\nEffective rank of TLF error curves: %.2f / 24\n", er_tlf))
cat(sprintf("Effective rank of FDA error curves:  %.2f / 24\n", er_fda))

rank_df <- tibble(
  model          = c(BENCHMARK, "FDA model"),
  effective_rank = round(c(er_tlf, er_fda), 3)
)
write_csv(rank_df, "output/tables/tlf_error_effective_rank.csv")

# Cross-hour correlation heatmap (TLF)
cor_long <- function(cm, label) {
  as.data.frame(cm) |>
    mutate(hour_row = 0:23) |>
    pivot_longer(-hour_row, names_to = "hour_col", values_to = "correlation") |>
    mutate(hour_col = as.integer(hour_col),
           model    = label)
}

cor_df <- bind_rows(
  cor_long(cor_tlf, BENCHMARK),
  cor_long(cor_fda, "FDA model")
)

p_cor_heat <- ggplot(cor_df, aes(x = hour_col, y = hour_row, fill = correlation)) +
  geom_tile() +
  facet_wrap(~model, ncol = 2) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-1, 1), name = "Corr") +
  scale_x_continuous(breaks = seq(0, 23, 4)) +
  scale_y_continuous(breaks = seq(0, 23, 4)) +
  labs(
    title    = "Cross-hour correlation of forecast errors",
    subtitle = "Low effective rank = errors co-move as curves, not independent scalars",
    x        = "Hour",
    y        = "Hour"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text  = element_text(size = 9),
        panel.grid  = element_blank())

ggsave("output/figs/tlf_error_crosshour_correlation.png", p_cor_heat,
       width = 10, height = 5, dpi = 150)
cat("Saved: output/figs/tlf_error_crosshour_correlation.png\n")

# Eigenvalue spectrum of TLF error covariance (how many dimensions dominate?)
sv_tlf <- svd(cov(tlf_mat))$d
sv_fda <- svd(cov(fda_mat))$d

sv_df <- bind_rows(
  tibble(component = 1:24, var_share = sv_tlf / sum(sv_tlf), model = BENCHMARK),
  tibble(component = 1:24, var_share = sv_fda / sum(sv_fda), model = "FDA model")
)

p_spectrum <- ggplot(sv_df |> filter(component <= 10),
                     aes(x = component, y = cumsum(var_share), colour = model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("FDA model" = "#2166ac",
                                  "Official forecast" = "#d6604d"),
                      name = NULL) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = "Cumulative variance explained by eigenvalues of error covariance",
    subtitle = "Steeper rise = errors are more curve-structured (low-dimensional)",
    x        = "Number of components",
    y        = "Cumulative variance share"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave("output/figs/tlf_error_eigenspectrum.png", p_spectrum,
       width = 7, height = 4, dpi = 150)
cat("Saved: output/figs/tlf_error_eigenspectrum.png\n")

cat("\nDone.\n")
