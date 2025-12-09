library(BSL)
library(stabledist)

# ==============================================================================
# 1. Helper: Gumbel V Generator
# ==============================================================================
rGumbV <- function(n, theta) {
    # Constraint check handled in prior/wrapper, but safety here:
    if (theta < 1) {
        return(NULL)
    }

    # Map Gumbel theta to stabledist parameters for Archimedean generator
    val_gamma <- (cos(pi / (2 * theta)))^theta

    V <- rstable(
        n = 1, alpha = 1 / theta, beta = 1,
        gamma = val_gamma, delta = 0, pm = 1
    )

    if (V <= 0) V <- .Machine$double.eps

    E <- rexp(n)
    U <- exp(-(E / V)^(1 / theta))
    return(U)
}

# ==============================================================================
# 2. Simulator Function (fnSim)
# ==============================================================================
fn_sim_gumbel <- function(theta, n_obs) {
    mu <- theta[1]
    sigma <- theta[2]
    theta_cop <- theta[3]

    # Generate Uniforms via Gumbel Copula
    U <- rGumbV(n_obs, theta_cop)

    # If stable generation failed (e.g. numerical issues), return NA
    if (is.null(U) || any(is.nan(U))) {
        return(NA)
    }

    # Inverse CDF (Lognormal)
    X <- qlnorm(U, meanlog = mu, sdlog = sigma)
    return(X)
}

# ==============================================================================
# 3. Summary Statistic Function (fnSum)
# ==============================================================================
fn_sum_stats <- function(x) {
    m <- median(x)
    md <- mad(x)
    if (md == 0) md <- 1e-6

    # Standardized Max (Skewed statistic)
    s_max <- (max(x) - m) / md
    return(c(m, md, s_max))
}

# ==============================================================================
# 4. Log-Prior Function (fnLogPrior)
# ==============================================================================
fn_log_prior <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    cop_theta <- theta[3]

    # --- Hard Constraints ---
    # Sigma must be positive, Theta must be >= 1
    if (sigma <= 0 || cop_theta < 1) {
        return(-Inf)
    }

    # --- Weakly Informative Priors ---

    # 1. Mu: Normal(0, sd=10)
    # Weakly informative for the location parameter
    lp_mu <- dnorm(mu, mean = 0, sd = 10, log = TRUE)

    # 2. Sigma: Gamma(shape=2, scale=2)
    # Ensures positivity, mode around 2, heavy tail allowed
    lp_sigma <- dgamma(sigma, shape = 2, scale = 2, log = TRUE)

    # 3. Theta: Shifted Gamma (Theta - 1 ~ Gamma)
    # Ensures Theta >= 1.
    # We calculate density of (theta - 1)
    lp_theta <- dgamma(cop_theta - 1, shape = 2, scale = 2, log = TRUE)

    # Sum the log probabilities
    return(lp_mu + lp_sigma + lp_theta)
}

# ==============================================================================
# 5. Setup Data and Model Object
# ==============================================================================

set.seed(123)
N_OBS <- 1000
true_params <- c(0, 1, 2) # mu=0, sigma=1, theta=2

# Generate "Observed" Data
y_obs <- fn_sim_gumbel(true_params, n_obs = N_OBS)

# ==============================================================================
# 6. Run semiBSL
# ==============================================================================

# Proposal Covariance (needs tuning in practice)
cov_rw <- diag(c(0.05, 0.05, 0.1))

library(parallel)

# 1. Detect cores (leave one free for OS)
# 2. Prepare Starting Values
start_points <- list(
    c(0.1, 1.0, 1.1),
    c(0.8, 2.0, 1.8),
    c(0.5, 1.2, 1.4),
    c(0.3, 1.8, 2.0)
)

n_chains <- length(start_points)

n_cores <- n_chains
cl <- makeCluster(n_cores)


# Ensure the list length matches the number of cores/tasks you want to run
if (length(start_points) > n_cores) start_points <- start_points[1:n_cores]

# 3. Export Libraries to the cluster
clusterEvalQ(cl, {
    library(BSL)
    library(stabledist)
})

# 4. Export Data and Functions to the cluster
# This is CRITICAL. The workers start blank.
clusterExport(cl, varlist = c(
    "y_obs", "cov_rw", # Data objects
    "rGumbV", "fn_sim_gumbel", # Simulation functions
    "fn_sum_stats", "fn_log_prior", # Summary and Prior
    "N_OBS"
))

# 5. Define the Wrapper Function
# This function runs on each worker
run_chain_worker <- function(start_val) {
    # --- Define the BSL Model Object ---
    my_model <- newModel(
        fnSim = fn_sim_gumbel,
        fnSum = fn_sum_stats,
        fnLogPrior = fn_log_prior,
        theta0 = start_val, # Initial guess

        # Auxiliary arguments for simulation
        simArgs = list(n_obs = N_OBS),

        # Parameter names for plotting
        thetaNames = c("mu", "sigma", "theta_cop")
    )

    # Run BSL
    res <- bsl(
        y = y_obs,
        n = 300,
        M = 5000,
        model = my_model,
        covRandWalk = cov_rw,
        method = "semiBSL",
        verbose = FALSE
    )
    return(res)
}

# 6. Run in Parallel
cat("Running chains in parallel on", n_cores, "cores...\n")
results_list <- parLapply(cl, start_points, run_chain_worker)

# 7. Stop Cluster
stopCluster(cl)
cat("Parallel execution finished.\n")

library(coda)

# 1. Extract samples and convert to 'mcmc.list' object
mcmc_list <- mcmc.list(lapply(results_list, function(res) {
    # extract samples, discard burn-in (e.g., first 1000)
    samps <- getTheta(res)
    mcmc(samps[-(1:1000), ])
}))

# 2. Visual Traceplot (Check if they mix/overlap)
plot(mcmc_list)

# 3. Gelman-Rubin Diagnostic
# Point est. should be close to 1 (e.g., < 1.05 or 1.1)
gelman.diag(mcmc_list)

# 4. Combine all chains for final inference
all_samples <- do.call(rbind, lapply(results_list, function(res) {
    getTheta(res)[-(1:1000), ]
}))

# Calculate Final Posterior Means
colMeans(all_samples)
