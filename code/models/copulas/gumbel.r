source("code/packages.R")

copula_gumbel <- list(

  name = "gumbel",
# --------------------------
  # Generator (phi) and inverse generator (psi)
  # --------------------------
  inv_psi = function(u, theta) (-log(u))^theta,
  psi = function(t, theta) exp(-t^(1/theta)),

  # Gumbel h-function (Conditional CDF) for Markov transitions

  h_dist = function(u, v, theta) {
    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    ln_u <- -log(u)
    ln_v <- -log(v)
    term <- ln_u^theta + ln_v^theta
    C <- exp(-term^(1 / theta))
    h <- C * (ln_v^(theta - 1)) * (term^(1 / theta - 1)) / v
    return(h)
  }
)
