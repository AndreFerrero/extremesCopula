# =============================================================================
# Test script: small MCMC run with Gumbel copula and Lognormal margin
# =============================================================================

# Load packages
source("libs/packages.R")

# Load models
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")

# Load builders
source("libs/models/builders/simulator.R")
source("libs/models/builders/logposterior.R")

# Load MCMC machinery
source("libs/mcmc/run_chain.R")
source("libs/mcmc/engines/metropolis_hastings.R")
source("libs/mcmc/proposals/gaussian_rw.R")
source("libs/mcmc/adaptation/none.R")
source("libs/mcmc/adaptation/haario.R")
source("libs/mcmc/adaptation/robbins_monro.R")

# =============================================================================
# 1. Simulate data using the simulator
# =============================================================================
param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

set.seed(123)
true_param <- c(mu = 0, sigma = 1, theta = 4)
n_obs <- 150
X <- simulator(true_param, n_obs)

# =============================================================================
# 2. Build log-posterior
# =============================================================================
g <- function(param) {
  c(param["mu"], log(param["sigma"]), log(param["theta"] - 1))
}

g_inv <- function(phi) {
  param <- c(
    mu    = phi[1],
    sigma = exp(phi[2]),
    theta = exp(phi[3]) + 1
  )
  # Ensure names are exactly what we want
  names(param) <- c("mu", "sigma", "theta")
  param
}


log_jacobian <- function(phi) phi[2] + phi[3]

logpost <- build_logposterior(
  copula = copula_gumbel,
  margin = margin_lognormal,
  param_map = param_map,
  data = X,
  inverse_transform = g_inv,
  log_jacobian = log_jacobian
)

# =============================================================================
# 3. Define proposal and run chain
# =============================================================================
phi_init <- g(c(mu = 0, sigma = 1, theta = 2))
Sigma0 <- matrix(c(
  1, 0, 0,
  0, 0.1, 0,
  0, 0, 0.1
), byrow = TRUE, ncol = 3, nrow = 3)

proposal <- proposal_gaussian_rw(Sigma = Sigma0)

res <- run_chain(
  log_target = logpost,
  init = phi_init,
  n_iter = 5000,
  proposal = proposal,
  burn_in = 0,
  adapt = adapt_none()
)

full_lik_dir <- here("sims", "estim", "joint", "full_lik_bayes")
full_lik_res_dir <- here(full_lik_dir, "res")
# save(res, file = here(full_lik_res_dir, "lognorm_gumbel_10kruns_150obs_adapthaario.Rdata"))
# load(here(full_lik_res_dir, "lognorm_gumbel_10kruns_150obs_adapthaario.Rdata"))

# Convert samples back to natural space
samples_natural <- t(apply(res$samples, 1, g_inv))
mcmc_samples <- mcmc(samples_natural)

# burn_in <- nrow(samples_natural)/2
burn_in <- 0
mcmc_clean <- window(mcmc_samples, start = burn_in + 1, thin = 1)

# Traceplot

plot(mcmc_clean)

res$conv_state

summary(mcmc_clean)
# =============================================================================
# 4. Quick summaries
# =============================================================================

# cat("Acceptance rate:", res$accept_rate, "\n")
# cat("Posterior means:\n")
# print(colMeans(as.matrix(mcmc_clean)))

# effectiveSize(mcmc_clean)

# # Arrange 1 row and 3 columns
# par(mfrow = c(1, 3), mar = c(4, 4, 2, 1)) # smaller margins for tighter plots

# plot(mcmc_clean[, "mu"],
#   type = "l",
#   ylab = "mu", xlab = "Iteration", main = "Trace of mu", density = FALSE, auto.layout = FALSE,
# )

# plot(mcmc_clean[, "sigma"],
#   type = "l",
#   ylab = "sigma", xlab = "Iteration", main = "Trace of sigma", density = FALSE, auto.layout = FALSE,
# )

# plot(mcmc_clean[, "theta"],
#   type = "l",
#   ylab = "theta", xlab = "Iteration", main = "Trace of theta", density = FALSE, auto.layout = FALSE,
# )

# # Reset par
# par(mfrow = c(1, 1))

# pairs(as.matrix(mcmc_clean))

# acf(as.matrix(mcmc_clean))
