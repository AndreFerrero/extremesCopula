source("libs/packages.R")

copula_gumbel <- list(

  name = "gumbel",
# --------------------------
  # Generator (phi) and inverse generator (psi)
  # --------------------------
  inv_psi = function(u, theta) (-log(u))^theta,
  psi = function(t, theta) exp(-t^(1/theta)),

  # --------------------------
  # 1. Simulate uniforms using latent variable method
  # --------------------------
  simulate_u = function(theta, n) {
    if (theta == 1) return(runif(n))  # independence
    if (theta < 1) return(NULL)       # invalid

    # --- simulate latent variable V ---
    val_gamma <- (cos(pi / (2 * theta)))^theta
    V <- stabledist::rstable(
      n = 1,
      alpha = 1 / theta,
      beta  = 1,
      gamma = val_gamma,
      delta = 0,
      pm    = 1
    )

    if (!is.finite(V) || V <= 0) return(NULL)

    # --- exponential variates ---
    E <- rexp(n)

    # --- apply inverse generator psi ---
    U <- copula_gumbel$psi(E / V, theta)
    U
  },

  log_density = function(u, theta) {
    copula::dCopula(u, copula::gumbelCopula(theta, dim = length(u)), log = TRUE)
  },

  log_prior = function(theta, a = 2, b = 1) {
    if (theta <= 1) return(-Inf)
    dgamma(theta - 1, a, b, log = TRUE)
  },

  diag = function(u, theta, n = length(u)) {
    
    if (theta < 1) return(NA)             # invalid parameter
    if (theta == 1) return(u^n)           # independence case
    
    u^(n^(1/theta))
}
)
