# =============================================================================
# STUDY OF DEPENDENCE IN EXCHANGEABLE GUMBEL–MARGIN MODEL
# Minimal, coherent simulation structure
# =============================================================================

# ------------------------------
# 1. Load required packages
# ------------------------------
source("libs/packages.R")

# ------------------------------
# 2. Source model components
# ------------------------------
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")
source("libs/models/builders/simulator.R")

library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------
# 3. Define parameter mapping
# ------------------------------
param_map <- list(
  margin = c("mu", "sigma"),
  copula = "theta"
)

# ------------------------------
# 4. Build the simulator
# ------------------------------
simulator <- build_simulator(
  copula = copula_gumbel,
  margin = margin_lognormal,
  param_map = param_map
)

# ------------------------------------------------------------
# Study setup
# ------------------------------------------------------------
set.seed(123)

n        <- 1000
n_rep    <- 1000
theta_vec <- c(1, 1.5, 2, 3, 4, 5)

param_base <- c(mu = 0, sigma = 1)
u_prob <- 0.95

# ------------------------------------------------------------
# Unified statistics extractor (ONE dataset → ALL statistics)
# ------------------------------------------------------------
compute_stats <- function(X, u_prob, param_base) {

  n <- length(X)
  X_ord <- sort(X)

  # Fixed marginal threshold
  u <- qlnorm(
    u_prob,
    meanlog = param_base["mu"],
    sdlog   = param_base["sigma"]
  )

  exceed <- which(X > u)

  c(
    max            = X_ord[n],
    spacing        = X_ord[n] - X_ord[n - 1],
    uncond         = mean(X > u),
    cond           = if (length(exceed) >= 2)
                       (length(exceed) - 1) / (n - 1)
                     else NA,
    exceed_count   = length(exceed),
    sample_mean    = mean(X)
  )
}

# ------------------------------------------------------------
# Run simulations for ONE theta
# ------------------------------------------------------------
run_theta_study <- function(theta) {

  stats_mat <- matrix(NA, nrow = n_rep, ncol = 6)
  colnames(stats_mat) <- c(
    "max", "spacing", "uncond", "cond",
    "exceed_count", "sample_mean"
  )

  for (r in 1:n_rep) {
    X <- simulator(c(param_base, theta = theta), n)
    stats_mat[r, ] <- compute_stats(X, u_prob, param_base)
  }

  as.data.frame(stats_mat)
}

# ------------------------------------------------------------
# Run full study
# ------------------------------------------------------------
results <- lapply(theta_vec, run_theta_study)
names(results) <- theta_vec

# ------------------------------------------------------------
# Aggregate results (means only — no bands)
# ------------------------------------------------------------
comparison <- data.frame(
  theta = theta_vec,
  max_mean = sapply(results, function(df) mean(df$max)),
  spacing_mean = sapply(results, function(df) mean(df$spacing)),
  uncond_mean = sapply(results, function(df) mean(df$uncond)),
  cond_mean = sapply(results, function(df) mean(df$cond, na.rm = TRUE)),
  dispersion = sapply(results, function(df)
    var(df$exceed_count) / mean(df$exceed_count)),
  var_mean = sapply(results, function(df) var(df$sample_mean))
)

# Ratio of conditional to unconditional exceedance probability
comparison$ratio_prob <- comparison$cond_mean / comparison$uncond_mean

print(comparison)

# ------------------------------------------------------------
# Prepare data for plotting
# ------------------------------------------------------------
plot_df <- comparison %>%
  select(
    theta,
    max_mean,
    spacing_mean,
    dispersion,
    cond_mean,
    var_mean,
    ratio_prob
  ) %>%
  pivot_longer(
    cols = -theta,
    names_to = "stat",
    values_to = "value"
  )

stat_labels <- c(
  max_mean     = "Mean Maximum",
  spacing_mean = "Mean Upper Spacing",
  dispersion   = "Dispersion of Exceedances",
  cond_mean    = "Conditional Exceedance Probability",
  var_mean     = "Variance of Sample Mean",
  ratio_prob   = "Conditional / Unconditional Exceedance Probability"
)

plot_df$stat_label <- stat_labels[plot_df$stat]

# ------------------------------------------------------------
# Final plot
# ------------------------------------------------------------
ggplot(plot_df, aes(x = theta, y = value,
                    color = stat_label,
                    group = stat_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ stat_label, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression(theta),
    y = "Statistic Value",
    title = "Effect of Dependence (θ) on Key Sample Statistics",
    color = "Statistic"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
