# ==============================================================================
# BAYESIAN DISTORTED GEV WITH RSTAN
# ==============================================================================
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("R/utils.R")
load_all()
# ------------------------------------------------------------------------------
# 1. Setup Data (Joe-Fréchet Simulation)
# ------------------------------------------------------------------------------
n <- 1000
k <- 100
theta_true <- 3
family_choice <- "joe" # Joe is non-GEV, so the package fit should be worse
seed <- 123

qfrechet <- function(u, xi = 0.25) {
  # alpha = 1/xi
  return( (-log(u))^(-xi) )
}

# ------------------------------------------------------------------------------
# Simulate
# ------------------------------------------------------------------------------
sim <- simulate_data(k, n, theta_true, family_choice, margin_fn = function(u) qfrechet(u, xi = 0.25), seed = seed)

# ------------------------------------------------------------------------------
# 3. Fit the Model (Two-Stage Bayesian)
# ------------------------------------------------------------------------------
# We use theta_fixed = 3 (In practice, use your theta_dmle result here)
stan_data <- list(
  K = k,
  M = sim$M,
  family = 2,           # Joe
  use_fixed_theta = 0,  # Fixed theta (Two-Stage)
  theta_fixed = 3.0
)

fit <- stan(
  file = "distortion.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)

# ------------------------------------------------------------------------------
# 4. Results and Diagnostics
# ------------------------------------------------------------------------------
# Print parameter summary
print(fit, pars = c("mu", "sigma", "xi", "theta", "return_level_100"))

# Traceplots for convergence
traceplot(fit, pars = c("mu", "sigma", "xi", "theta", "return_level_100"))

# Posterior distribution of the Return Level
plot(fit, pars = "return_level_100", show_density = TRUE)

# Extract samples for manual analysis
post_samples <- as.data.frame(fit)
cat("Posterior Median for xi:", median(post_samples$xi), "\n")
