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
fn_sim_gumbel_lognorm <- function(theta, n_obs) {
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
    lp_theta <- dgamma(cop_theta - 1, shape = 2, scale = 1, log = TRUE)

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
y_obs <- fn_sim_gumbel_lognorm(true_params, n_obs = N_OBS)

# ==============================================================================
# 6. Run semiBSL
# ==============================================================================
cov_rw <- diag(c(0.03, 0.03, 0.03))

mod <- newModel(
    fnSim = fn_sim_gumbel_lognorm,
    fnSum = fn_sum_stats,
    theta0 = c(0, 1, 2),
    fnLogPrior = fn_log_prior,
    simArgs = list(n_obs = N_OBS)
)

bounds <- matrix(c(-Inf, Inf, 0, Inf, 1, Inf), 3, 2, byrow = TRUE)

result_semibsl <- bsl(
    y = y_obs,
    n = 300,
    M = 10000,
    model = mod,
    covRandWalk = cov_rw,
    method = "semiBSL",
    logitTransformBound = bounds,
    verbose = FALSE
)


plot(result_semibsl, which = 1, thetaTrue = true_params)

result_semibsl
summary(result_semibsl)

theta_chain <- result_semibsl@theta

# Traceplots
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

plot(theta_chain[,1], type = "l",
     main = expression(mu), xlab = "Iteration", ylab = "")
plot(theta_chain[,2], type = "l",
     main = expression(sigma), xlab = "Iteration", ylab = "")
plot(theta_chain[,3], type = "l",
     main = expression(theta), xlab = "Iteration", ylab = "")
par(mfrow = c(1,1))
