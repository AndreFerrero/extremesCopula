# ==============================================================================
# Semi-Parametric BSL from Scratch (An et al., 2019)
# ==============================================================================

library(parallel)
library(stabledist) # For Gumbel Copula simulation
library(mvtnorm)    # For multivariate normal proposal
library(MASS)       # For ginv (pseudo-inverse if correlation matrix is singular)

# ==============================================================================
# 1. Model Definitions (Simulator & Statistics)
# ==============================================================================

# 1.1 Gumbel Generator (Stable distribution parameterization)
rGumbV <- function(n, theta) {
  # theta must be >= 1. If theta=1, it's independence.
  if (theta < 1) return(rep(NA, n))
  
  # Sampling V from Positive Stable Distribution
  # alpha = 1/theta, beta = 1, gamma = (cos(pi/(2*theta)))^theta, delta = 0
  val_gamma <- (cos(pi / (2 * theta)))^theta
  V <- stabledist::rstable(n = 1, alpha = 1 / theta, beta = 1,
                           gamma = val_gamma, delta = 0, pm = 1)
  
  if (V <= 0) V <- .Machine$double.eps # Numerical safety
  
  E <- rexp(n)
  U <- exp(-(E / V)^(1 / theta))
  return(U)
}

# 1.2 Full Simulator
# Params: theta[1] = mu, theta[2] = sigma, theta[3] = copula_theta
fn_sim_gumbel <- function(theta, n_obs) {
  mu <- theta[1]
  sigma <- theta[2]
  cop_theta <- theta[3]
  
  # Safety checks
  if (sigma <= 0 || cop_theta < 1) return(NA)
  
  # Simulate Copula (Dependence)
  U <- tryCatch(rGumbV(n_obs, cop_theta), error = function(e) return(NULL))
  if (is.null(U) || any(is.na(U))) return(NA)
  
  # Inverse Probability Integral Transform for Marginals (LogNormal)
  # Q(u) of LogNorm(mu, sigma)
  X <- qlnorm(U, meanlog = mu, sdlog = sigma)
  
  # Sanity check for infinite values
  if (any(is.infinite(X))) return(NA)
  
  return(X)
}

# 1.3 Summary Statistics

fn_sum_stats <- function(x) {
  if (any(is.na(x))) return(rep(NA, 5)) # Increased dimension to 5
  
  # 1. Location
  m <- median(x)
  
  # 2. Spread (Robust)
  md <- mad(x)
  
  # 3. Spread (Tail-sensitive) - Helps separate Sigma from Theta
  # The difference between 95th and 5th percentile captures the width of the main body
  q_low <- quantile(x, 0.05)
  q_high <- quantile(x, 0.95)
  q_diff <- as.numeric(q_high - q_low)
  
  # 4. Skewness/Shape check (Log-Normal specific)
  # Ratio of upper tail to lower tail spread
  tail_ratio <- as.numeric((q_high - m) / (m - q_low))
  
  # 5. Extreme dependency (Standardized Max)
  if (md == 0) md <- 1e-6 
  s_max <- (max(x) - m) / md
  
  return(c(m, md, q_diff, tail_ratio, s_max))
}

# 1.4 Prior Distribution (Weakly Informative)
# Input: theta in Natural Scale
fn_log_prior <- function(theta) {
  mu <- theta[1]
  sigma <- theta[2]
  cop_theta <- theta[3]
  
  # Hard Bounds
  if (sigma <= 0 || cop_theta <= 1) return(-Inf)
  
  # Priors:
  # mu ~ N(0, 10^2)
  # sigma ~ Cauchy(0, 2.5) truncated > 0 (Weakly informative for scale)
  # cop_theta ~ Gamma(shape=2, rate=0.5) shifted by 1 (Since theta > 1)
  
  lp_mu <- dnorm(mu, mean = 0, sd = 10, log = TRUE)
  lp_sigma <- dcauchy(sigma, location = 0, scale = 2.5, log = TRUE)
  lp_theta <- dgamma(cop_theta - 1, shape = 2, rate = 0.5, log = TRUE)
  
  return(lp_mu + lp_sigma + lp_theta)
}

# ==============================================================================
# 2. Semi-Parametric Likelihood Core (The Math)
# ==============================================================================

# Helper: Gaussian Kernel Density Estimation at specific points
# Returns density estimate and CDF estimate
get_kde_estimates <- function(sim_data, obs_val) {
  # Silverman's rule of thumb for bandwidth
  n <- length(sim_data)
  sd_sim <- sd(sim_data)
  iqr_sim <- IQR(sim_data)
  
  # Handle degenerate case where variance is 0
  min_spread <- min(sd_sim, iqr_sim / 1.34)
  if (min_spread == 0) min_spread <- sd_sim
  if (min_spread == 0) min_spread <- 1e-6
  
  h <- 0.9 * min_spread * n^(-0.2)
  
  # 1. Marginal Log-Density of Observed Stat (using Gaussian Kernel)
  # f_hat(y) = (1/nh) * sum K((y - x_i)/h)
  # We work in log space for stability
  # K is standard normal PDF
  z_scores <- (obs_val - sim_data) / h
  # LogSumExp trick for stability: log(mean(exp(vals)))
  log_vals <- dnorm(z_scores, log = TRUE)
  max_val <- max(log_vals)
  log_pdf <- max_val + log(sum(exp(log_vals - max_val))) - log(n) - log(h)
  
  # 2. PIT (Probability Integral Transform) for Observed Stat
  # F_hat(y) = (1/n) * sum Phi((y - x_i)/h)
  u_obs <- mean(pnorm(z_scores))
  
  # 3. PIT for Simulated Stats (Leave-one-out not strictly necessary for MCMC but cleaner)
  # Here we just use the empirical CDF smoothened by the kernel for the simulations
  # For the correlation matrix, we need the "Normal Scores" of the simulations.
  # An et al. suggest using empirical ranks or kernel CDF. 
  # We use Empirical CDF scaled by n/(n+1) for robustness in converting to Z-scores.
  u_sim <- rank(sim_data) / (n + 1)
  
  return(list(log_pdf = log_pdf, u_obs = u_obs, u_sim = u_sim))
}

# The SemiBSL Log-Likelihood Function
semi_bsl_loglik <- function(theta, y_sum, n_sim, n_obs) {
  
  # 1. Simulate Data
  # We generate n_sim datasets, compute summary stats for each
  sim_stats_mat <- matrix(NA, nrow = n_sim, ncol = length(y_sum))
  
  for(i in 1:n_sim) {
    x_sim <- fn_sim_gumbel(theta, n_obs)
    if (any(is.na(x_sim))) return(-Inf) # Reject if simulation fails
    sim_stats_mat[i, ] <- fn_sum_stats(x_sim)
  }
  
  # Check for NAs in summary stats
  if (any(is.na(sim_stats_mat))) return(-Inf)
  
  d <- length(y_sum)
  total_marginal_loglik <- 0
  eta_obs <- numeric(d)
  eta_sim <- matrix(NA, nrow = n_sim, ncol = d)
  
  # 2. Marginals: Fit KDE and Transform
  for(j in 1:d) {
    kde_res <- get_kde_estimates(sim_stats_mat[, j], y_sum[j])
    
    # Add log density of marginal
    total_marginal_loglik <- total_marginal_loglik + kde_res$log_pdf
    
    # Transform Obs to Normal Space (clip to avoid Inf)
    u_o <- max(min(kde_res$u_obs, 1 - 1e-6), 1e-6)
    eta_obs[j] <- qnorm(u_o)
    
    # Transform Sims to Normal Space
    # u_sim is already safe via rank/(n+1)
    eta_sim[, j] <- qnorm(kde_res$u_sim)
  }
  
  # 3. Dependence: Gaussian Copula
  # Estimate Correlation Matrix R from transformed simulations
  R_hat <- cor(eta_sim)
  
  # Regularize R if singular (rare but possible with low n)
  if (det(R_hat) <= 1e-9) {
    R_hat <- R_hat + diag(1e-6, d)
  }
  
  # Calculate Gaussian Copula Density (Eq 3 in Paper)
  # log(c(u)) = -0.5 * log|R| - 0.5 * eta_obs' * (R^-1 - I) * eta_obs
  
  R_inv <- tryCatch(solve(R_hat), error = function(e) MASS::ginv(R_hat))
  
  log_det_R <- as.numeric(determinant(R_hat, logarithm = TRUE)$modulus)
  
  quadratic_term <- as.numeric(t(eta_obs) %*% (R_inv - diag(1, d)) %*% eta_obs)
  
  copula_loglik <- -0.5 * log_det_R - 0.5 * quadratic_term
  
  # Final Likelihood
  return(total_marginal_loglik + copula_loglik)
}

# ==============================================================================
# 3. MCMC Sampler Infrastructure (Adaptive MH)
# ==============================================================================

# Transformations to Unconstrained Space
# mu -> mu (Identity)
# sigma -> log(sigma)
# theta_cop -> log(theta_cop - 1)
to_unconstrained <- function(theta) {
  c(theta[1], log(theta[2]), log(theta[3] - 1))
}

to_natural <- function(phi) {
  c(phi[1], exp(phi[2]), exp(phi[3]) + 1)
}

# Jacobian Adjustment for Prior (if prior is defined on natural scale)
# log J = log(sigma) + log(theta_cop - 1)
log_jacobian <- function(phi) {
  phi[2] + phi[3] 
}

run_adaptive_chain <- function(y_obs_sum, n_mcmc, n_sim, n_obs, start_theta) {
  
  # MCMC Settings
  d_param <- 3
  burn_in <- floor(n_mcmc * 0.2)
  
  # Storage
  chain <- matrix(NA, nrow = n_mcmc, ncol = d_param)
  log_post_store <- numeric(n_mcmc)
  acc_rate <- 0
  
  # Initialize
  curr_phi <- to_unconstrained(start_theta)
  curr_theta <- start_theta
  
  # Initial Likelihood
  curr_ll <- semi_bsl_loglik(curr_theta, y_obs_sum, n_sim, n_obs)
  curr_lp <- fn_log_prior(curr_theta)
  curr_post <- curr_ll + curr_lp + log_jacobian(curr_phi)
  
  # Check initialization
  if(!is.finite(curr_post)) stop("Initial parameter invalid (Inf posterior)")
  
  # Adaptive Covariance Variables
  # Initial proposal covariance (small diagonal)
  cov_prop <- diag(0.1, d_param) 
  epsilon <- 1e-6
  s_d <- 2.38^2 / d_param # Optimal scaling constant
  
  # Loop
  for (i in 1:n_mcmc) {
    
    # 1. Propose new parameter (in unconstrained space)
    prop_phi <- as.vector(mvtnorm::rmvnorm(1, mean = curr_phi, sigma = cov_prop))
    prop_theta <- to_natural(prop_phi)
    
    # 2. Evaluate Prior + Jacobian
    prop_lp <- fn_log_prior(prop_theta)
    prop_jac <- log_jacobian(prop_phi)
    
    if (prop_lp == -Inf) {
      log_alpha <- -Inf
    } else {
      # 3. Evaluate SemiBSL Likelihood
      prop_ll <- semi_bsl_loglik(prop_theta, y_obs_sum, n_sim, n_obs)
      
      prop_post <- prop_ll + prop_lp + prop_jac
      
      # 4. Metropolis Ratio
      log_alpha <- prop_post - curr_post
    }
    
    # 5. Accept/Reject
    if (log(runif(1)) < log_alpha) {
      curr_phi <- prop_phi
      curr_theta <- prop_theta
      curr_post <- prop_post
      curr_ll <- prop_ll
      acc_rate <- acc_rate + 1
    }
    
    # Store
    chain[i, ] <- curr_theta
    log_post_store[i] <- curr_post
    
    # 6. Adaptive Covariance Update (Haario et al. 2001)
    # Start adapting after some samples, strictly during burn-in (or continue if desired)
    if (i > 50 && i <= burn_in) {
      # Recursive formula or simple re-calculation on history
      # Using simple window re-calc for robustness in this script
      hist_window <- chain[1:i, ]
      # Transform history to unconstrained for covariance calculation
      hist_phi <- t(apply(hist_window, 1, to_unconstrained))
      
      cov_mat <- cov(hist_phi)
      cov_prop <- s_d * cov_mat + s_d * epsilon * diag(d_param)
    }
    
    # Print progress
    if (i %% 100 == 0) {
      cat(sprintf("Iter: %d | Acc: %.2f | Post: %.2f | Theta: %.2f, %.2f, %.2f\n", 
                  i, acc_rate/i, curr_post, curr_theta[1], curr_theta[2], curr_theta[3]))
    }
  }
  
  return(list(chain = chain, log_post = log_post_store, acc_rate = acc_rate/n_mcmc))
}

# ==============================================================================
# 4. Master Execution Block
# ==============================================================================

# Setup Data
set.seed(2025)
N_OBS <- 1000      # Data size
TRUE_THETA <- c(0, 1.0, 2.0) # mu=0, sigma=1, cop_theta=2
y_data <- fn_sim_gumbel(TRUE_THETA, N_OBS)
y_sum <- fn_sum_stats(y_data)

cat("Observed Stats:", round(y_sum, 4), "\n")

# Parallel Setup
n_chains <- 4
n_cores <- min(n_chains, parallel::detectCores() - 1)
n_iter <- 8000     # Total iterations
n_sim_bsl <- 1000   # Number of simulations for BSL likelihood estimation (n in paper)

start_points <- list(
  c(0.5, 0.8, 1.5),
  c(-0.5, 1.2, 2.5),
  c(0.1, 1.0, 1.8),
  c(-0.1, 1.5, 3.0)
)

# Create Cluster
cl <- makeCluster(n_cores, outfile = "")

# Load libraries and export functions to cores
clusterEvalQ(cl, {
  library(stabledist)
  library(mvtnorm)
  library(MASS)
})

clusterExport(cl, varlist = c("rGumbV", "fn_sim_gumbel", "fn_sum_stats", 
                              "fn_log_prior", "get_kde_estimates", 
                              "semi_bsl_loglik", "to_unconstrained", 
                              "to_natural", "log_jacobian", 
                              "run_adaptive_chain", "y_sum", 
                              "N_OBS", "n_iter", "n_sim_bsl"))

# Run Parallel Chains
cat("Starting Parallel MCMC...\n")
results <- parLapply(cl, start_points, function(start_val) {
  run_adaptive_chain(y_obs_sum = y_sum, 
                     n_mcmc = n_iter, 
                     n_sim = n_sim_bsl, 
                     n_obs = N_OBS, 
                     start_theta = start_val)
})

stopCluster(cl)

sbi_dir <- here("sims", "estim", "joint", "SBI")
res_dir <- here(sbi_dir, "res")

# save(results, file = here(res_dir, "morestats_semiBSL_lognorm.Rdata"))
# ==============================================================================
# 5. Diagnostics & Visualization
# ==============================================================================
load(here(res_dir, "morelogstats_semiBSL_lognorm.Rdata"))
library(coda)

# Convert to mcmc.list
mcmc_res <- mcmc.list(lapply(results, function(x) {
  # Remove Burn-in 
  burn <- n_iter/3
  mcmc(x$chain[-(1:burn), ])
}))

summary(mcmc_res)
plot(mcmc_res)

# Gelman-Rubin Convergence Diagnostic
gelman.diag(mcmc_res)

cat("True Theta:", TRUE_THETA, "\n")

##### PPC
# 1. Sample 100 parameters from your posterior
post_samples <- as.matrix(mcmc_res) # Combine all chains
idxs <- sample(nrow(post_samples), 100)
thetas_post <- post_samples[idxs, ]

# 2. Simulate data for each posterior sample
# We will plot the density of the OBSERVED data vs POSTERIOR data
plot(density(y_data), lwd=2, main="Posterior Predictive Check", ylim=c(0, 0.6))

for(i in 1:100) {
  # Simulate using the parameters from the chain
  y_rep <- fn_sim_gumbel(thetas_post[i, ], N_OBS)
  lines(density(y_rep), col=rgb(0, 0, 1, 0.1)) # Blue, transparent lines
}

# Redraw observed data on top
lines(density(y_data), lwd=2, col="black")
legend("topright", legend=c("Observed", "Posterior Predictive"), 
       col=c("black", "blue"), lwd=2)
