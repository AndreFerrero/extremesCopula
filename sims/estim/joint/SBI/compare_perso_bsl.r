# =============================================================================
# Full comparison script:
# Modular semiBSL vs BSL::bsl (adaptive vs non-adaptive)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
library(coda)
library(BSL)
library(stabledist)
library(here)

# -----------------------------------------------------------------------------
# 1. Load modular implementation
# -----------------------------------------------------------------------------
source("libs/packages.R")

source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")

source("libs/models/builders/simulator.R")
source("libs/models/builders/bsl_logposterior.R")
source("libs/models/builders/synthetic_semibsl.r")

source("libs/mcmc/run_chain.R")
source("libs/mcmc/engines/metropolis_hastings.R")
source("libs/mcmc/proposals/gaussian_rw.R")
source("libs/mcmc/adaptation/none.R")
source("libs/mcmc/adaptation/haario.R")

# -----------------------------------------------------------------------------
# 2. Data generation (shared)
# -----------------------------------------------------------------------------
set.seed(123)

true_param <- c(mu = 0, sigma = 1, theta = 2)
n_obs <- 1000

param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

X <- simulator(true_param, n_obs)

# -----------------------------------------------------------------------------
# 3. Summary statistics
# -----------------------------------------------------------------------------
med_mad_max <- function(x) {

  # Required for BSL compatibility
  if (length(x) == 1 && is.na(x)) {
    return(rep(NA_real_, 3))
  }

  if (any(!is.finite(x))) {
    return(rep(NA_real_, 3))
  }

  m <- median(x)
  md <- mad(x)

  if (!is.finite(md) || md == 0) {
    md <- 1e-6
  }

  raw_max <- max(x)

  c(m, md, raw_max)
}


# -----------------------------------------------------------------------------
# 4. Parameter transforms (modular code)
# -----------------------------------------------------------------------------
g <- function(param) {
  c(param["mu"], log(param["sigma"]), log(param["theta"] - 1))
}

g_inv <- function(phi) {
  param <- c(
    mu    = phi[1],
    sigma = exp(phi[2]),
    theta = exp(phi[3]) + 1
  )
  names(param) <- c("mu", "sigma", "theta")
  param
}

log_jacobian <- function(phi) phi[2] + phi[3]

phi_init <- g(true_param)

# -----------------------------------------------------------------------------
# 5. Build modular semiBSL log-posterior
# -----------------------------------------------------------------------------
logpost <- build_bsl_logposterior(
  copula = copula_gumbel,
  margin = margin_lognormal,
  param_map = param_map,
  data = X,
  simulator = simulator,
  sum_stats = med_mad_max,
  n_sim = 150,
  synthetic_loglik = synthetic_semibsl,
  inverse_transform = g_inv,
  log_jacobian = log_jacobian
)

# -----------------------------------------------------------------------------
# 6. Proposal
# -----------------------------------------------------------------------------
p <- length(phi_init)
Sigma0 <- diag(1, p) * 0.03
proposal <- proposal_gaussian_rw(Sigma0 = Sigma0)

# =============================================================================
# 7. Run modular chains
# =============================================================================

# --- Non-adaptive ---
set.seed(2024)

res_mod_none <- run_chain(
  log_target = logpost,
  init = phi_init,
  n_iter = 20000,
  proposal = proposal,
  adapt = adapt_none()
)

samples_mod_none <- t(apply(res_mod_none$samples, 1, g_inv))
mcmc_mod_none <- window(
  mcmc(samples_mod_none),
  start = 1
)

# --- Haario adaptive ---
set.seed(2024)

res_mod_haario <- run_chain(
  log_target = logpost,
  init = phi_init,
  n_iter = 20000,
  proposal = proposal,
  adapt = adapt_haario()
)

samples_mod_haario <- t(apply(res_mod_haario$samples, 1, g_inv))
mcmc_mod_haario <- window(
  mcmc(samples_mod_haario),
  start = 1
)

# =============================================================================
# 8. BSL::bsl implementation (non-adaptive)
# =============================================================================
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

fn_sim_gumbel_lognorm <- function(theta, n_obs) {
    mu <- theta[1]
    sigma <- theta[2]
    theta_cop <- theta[3]

    # Generate Uniforms via Gumbel Copula
    U <- rGumbV(n_obs, theta_cop)

    # If stable generation failed (e.g. numerical issues), return NA
    # if (is.null(U) || any(is.nan(U))) {
    #     return(NA)
    # }

    # Inverse CDF (Lognormal)
    X <- qlnorm(U, meanlog = mu, sdlog = sigma)
    return(X)
}

fn_log_prior <- function(theta) {
  if (theta[2] <= 0 || theta[3] < 1) return(-Inf)
  dnorm(theta[1], 0, 10, log = TRUE) +
    dgamma(theta[2], 2, 2, log = TRUE) +
    dgamma(theta[3] - 1, 2, 1, log = TRUE)
}

mod <- newModel(
  fnSim = fn_sim_gumbel_lognorm,
  fnSum = med_mad_max,
  theta0 = true_param,
  fnLogPrior = fn_log_prior,
  simArgs = list(n_obs = n_obs)
)

bounds <- matrix(c(-Inf, Inf, 0, Inf, 1, Inf), 3, 2, byrow = TRUE)

set.seed(2024)

result_bsl <- bsl(
  y = X,
  n = 150,
  M = 20000,
  model = mod,
  covRandWalk = Sigma0,
  method = "semiBSL",
  logitTransformBound = bounds,
  verbose = FALSE
)

theta_bsl <- result_bsl@theta
mcmc_bsl <- window(
  mcmc(theta_bsl),
  start = 1
)

# =============================================================================
# 9. Diagnostics
# =============================================================================

cat("\n================ ACCEPTANCE RATES ================\n")
cat("Modular (non-adapt):", res_mod_none$accept_rate, "\n")
cat("Modular (Haario)  :", res_mod_haario$accept_rate, "\n")
cat("BSL (fixed RW)    :", result_bsl@acceptRate, "\n")

cat("\n================ ESS =============================\n")
print(rbind(
  Modular_None   = effectiveSize(mcmc_mod_none),
  Modular_Haario = effectiveSize(mcmc_mod_haario),
  BSL            = effectiveSize(mcmc_bsl)
))

cat("\n================ POSTERIOR MEANS =================\n")
print(rbind(
  Modular_None   = colMeans(as.matrix(mcmc_mod_none)),
  Modular_Haario = colMeans(as.matrix(mcmc_mod_haario)),
  BSL            = colMeans(as.matrix(mcmc_bsl)),
  True           = true_param
))

# =============================================================================
# 10. Traceplots
# =============================================================================

par(mfrow = c(3, 3), mar = c(3, 4, 2, 1))

for (j in 1:3) {
  plot(as.matrix(mcmc_mod_none)[, j], type = "l",
       main = paste(colnames(mcmc_mod_none)[j], "(Mod None)"),
       xlab = "", ylab = "")
  plot(as.matrix(mcmc_mod_haario)[, j], type = "l",
       main = paste(colnames(mcmc_mod_haario)[j], "(Mod Haario)"),
       xlab = "", ylab = "")
  plot(as.matrix(mcmc_bsl)[, j], type = "l",
       main = paste(colnames(mcmc_bsl)[j], "(BSL)"),
       xlab = "", ylab = "")
}

par(mfrow = c(1,1))
