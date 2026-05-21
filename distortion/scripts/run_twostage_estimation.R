# ==============================================================================
# TWO-STAGE ESTIMATION SCRIPT
# ==============================================================================
# setwd("distortion")
source("R/utils.R")
load_all()

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

n <- 100
k <- 500
theta_true <- 2.5
family_choice <- "gumbel"
seed <- 123

# ------------------------------------------------------------------------------
# Simulate
# ------------------------------------------------------------------------------

sim <- simulate_data(k, n, theta_true, family_choice, seed = seed)
pseudo <- compute_pseudo_obs(sim$X)

# ------------------------------------------------------------------------------
# Stage 1: DMLE for theta
# ------------------------------------------------------------------------------

theta_dmle <- estimate_theta_dmle(pseudo$Yhat, n, k, family_choice)
cat("\nStage 1 - DMLE theta:", theta_dmle, "\n")

# ------------------------------------------------------------------------------
# Stage 2: GEV estimation with fixed theta
# ------------------------------------------------------------------------------

fit <- fit_twostage_gev(sim$M, theta_dmle, family_choice)

print_results(fit, theta_true)

# ------------------------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------------------------

# plot_diagnostics(sim$M, fit$mu, fit$sigma, fit$xi, fit$theta, family_choice)

# Bootstrap QQ bands
plot_qq_bootstrap(sim$M, fit, n, k, B = 100)
