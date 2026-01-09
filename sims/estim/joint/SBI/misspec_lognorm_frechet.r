# =============================================================================
# Test script: small MCMC run with Gumbel copula and Lognormal margin
# =============================================================================

# Load packages
source("libs/packages.R")

# Load models
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")
source("libs/models/margins/frechet.R")
source("libs/models/builders/maxima_distribution.r")

# Load builders
source("libs/models/builders/simulator.R")
source("libs/models/builders/bsl_logposterior.R")
source("libs/models/builders/synthetic_semibsl.r")
source("libs/models/builders/synthetic_gaussian_bsl.r")

# Load MCMC machinery
source("libs/mcmc/run_chain.R")
source("libs/mcmc/run_parallel_chain.R")
source("libs/mcmc/engines/metropolis_hastings.R")
source("libs/mcmc/proposals/gaussian_rw.R")
source("libs/mcmc/adaptation/none.R")
source("libs/mcmc/adaptation/haario.R")
source("libs/mcmc/adaptation/robbins_monro.R")
source("libs/mcmc/adaptation/hybrid_haario_rm.R")

#=============================================================================
# 1. Simulate data using the simulator
# =============================================================================
param_map <- list(margin = c("mu", "sigma"), copula = "theta")
param_map_true <- list(margin = c("shape", "scale"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)
simulator_true <- build_simulator(
  copula = copula_gumbel,
  margin = margin_frechet,
  param_map = param_map_true
)

set.seed(123)
true_param <- c(mu = 0, sigma = 1, theta = 2)
true_param_frechet <- c(shape = 2, scale = 1, theta = 2)
n_obs <- 1000
X <- simulator_true(true_param_frechet, n_obs)

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

med_mad_max <- function(x) {
  m <- median(x)
  md <- mad(x)
  if (md == 0) md <- 1e-6
  raw_max <- max(x)
  c(m, md, raw_max)
}

log_jacobian <- function(phi) phi[2] + phi[3]

logpost <- build_bsl_logposterior(
  copula = copula_gumbel, margin = margin_lognormal,
  param_map = param_map,
  data = X,
  simulator = simulator,
  sum_stats = med_mad_max,
  n_sim = 200,
  synthetic_loglik = synthetic_semibsl,
  inverse_transform = g_inv,
  log_jacobian = log_jacobian
)

# =============================================================================
# 3. Define proposal and run chain
# =============================================================================
phi_init <- g(c(mu = 0, sigma = 1, theta = 2))
p <- 3
# Sigma0 <- diag(0.001, p, p)
Sigma0 <- matrix(
  c(0.56, -0.12, 0.038,
    -0.12, 0.08, 0.09,
    0.038, 0.09, 0.23), nrow = p, ncol = p, byrow = TRUE)

proposal <- proposal_gaussian_rw(Sigma0 = Sigma0)

res <- run_chain(
  log_target = logpost,
  init = phi_init,
  n_iter = 2000,
  proposal = proposal,
  adapt = adapt_haario()
)

# init_list <- list(
#   init_1 = g(c(mu = 0, sigma = 1, theta = 2)),
#   init_2 = g(c(mu = 0.5, sigma = 0.5, theta = 1.5)),
#   init_3 = g(c(mu = -0.5, sigma = 1.5, theta = 2.5)),
#   init_4 = g(c(mu = 1, sigma = 2, theta = 3))
# )

# res_par <- run_parallel_chains(
#   log_target = logpost,
#   init_values = init_list,
#   n_iter = 20000,
#   proposal = proposal,
#   n_cores = 4,
#   adapt = adapt_none(),
#   transform = g_inv,
#   export = c("g_inv", "margin_lognormal", "copula_gumbel",
#   "simulator", "synthetic_semibsl", "log_jacobian", "Sigma0", "mh_step")
# )

sbi_dir <- here("sims", "estim", "joint", "SBI")
sbi_res_dir <- here(sbi_dir, "res")

save(res, Sigma0, file = here(sbi_res_dir, "semibsl_20kruns_200sims_1kobs_adaptnone_sigma0finalhybrid.Rdata"))

# save(res_par, Sigma0, init_list, file = here(sbi_res_dir, "semibsl_4chains_20kruns_200sims_1kobs_adaptnone_sigma0finalhybrid.Rdata"))

load(here(sbi_res_dir, "semibsl_20kruns_200sims_1kobs_adaptnone_sigma0finalhybrid.Rdata"))

# Convert samples back to natural space
samples_natural <- t(apply(res$samples, 1, g_inv))

mcmc_samples <- mcmc(samples_natural)

# burn_in <- nrow(samples_natural)/3
burn_in <- 0
mcmc_clean <- window(mcmc_samples, start = burn_in + 1, thin = 1)

# mcmc_clean_par <- window(res_par, start = burn_in + 1, thin = 1)
# =============================================================================
# 4. Quick summaries
# =============================================================================
res$conv

cat("Acceptance rate:", res$accept_rate, "\n")

summary(mcmc_clean)

# Traceplot
plot(mcmc_clean)

# =============================================================================
# 4. Recover Limit Distribution G
# =============================================================================

G_dist <- build_maxima_distribution(
  margin = margin_lognormal,
  copula = copula_gumbel,
  param_map = param_map
)

G_dist_true <- build_maxima_distribution(
  margin = margin_frechet,
  copula = copula_gumbel,
  param_map = param_map_true
)

y_grid <- seq(0.01, 40, length.out = 200)

G <- t(
  apply(
    as.matrix(mcmc_clean),
    1,
    G_dist,
    y = y_grid,
    n = n_obs
  )
)

G_true <- G_dist_true(
  true_param_frechet,
  y = y_grid,
  n = n_obs
)

summarize_G <- function(G_mat) {
  list(
    mean   = colMeans(G_mat),
    median = apply(G_mat, 2, median),
    lower  = apply(G_mat, 2, quantile, 0.1),
    upper  = apply(G_mat, 2, quantile, 0.9)
  )
}

stats <- summarize_G(G)

plot(y_grid, G_true, type = "l", lwd = 2, lty = 2, col = "red",
     ylim = c(0, 1), xlab = "y", ylab = expression(P(M[n] <= y)),
     main = "Maxima Distribution G, theta = 2")

# Plot Dependent Model Posterior
lines(y_grid, stats$mean, col = "blue", lwd = 2)

lines(y_grid, stats$median, col = "orange", lwd = 2)

polygon(c(y_grid, rev(y_grid)), c(stats$lower, rev(stats$upper)),
        col = rgb(0, 0, 1, 0.1), border = NA)

legend("bottomright",
       legend = c("True G",
                  "Posterior Mean",
                  "Posterior Median",
                  "80% Credible Interval"),
       col = c("red", "blue", "orange", rgb(0, 0, 1, 0.1)),
       lwd = c(2, 2, 2, NA),
       lty = c(2, 1, 1, NA),
       pch = c(NA, NA, NA, 15),
       pt.cex = 2,
       bty = "n")

# ----------------------------
# Full-grid ISE
# ----------------------------
G_mean <- stats$mean
G_med  <- stats$median

squared_diff_mean <- (G_mean - G_true)^2
squared_diff_med  <- (G_med - G_true)^2

dy <- diff(y_grid)[1]      # uniform grid spacing
ISE_mean <- sum(squared_diff_mean) * dy
ISE_med  <- sum(squared_diff_med)  * dy

cat("Full-grid ISE (Mean):", ISE_mean, "\n")
cat("Full-grid ISE (Median):", ISE_med, "\n")

# ----------------------------
# Tail ISE (e.g., G_true > 0.9)
# ----------------------------
tail_idx <- which(G_true > 0.9)          # indices for upper tail

squared_diff_mean_tail <- squared_diff_mean[tail_idx]
squared_diff_med_tail  <- squared_diff_med[tail_idx]

# compute tail-grid spacing
dy_tail <- diff(y_grid[tail_idx])[1]

ISE_mean_tail <- sum(squared_diff_mean_tail) * dy_tail
ISE_med_tail  <- sum(squared_diff_med_tail) * dy_tail

cat("Tail ISE (Mean):", ISE_mean_tail, "\n")
cat("Tail ISE (Median):", ISE_med_tail, "\n")