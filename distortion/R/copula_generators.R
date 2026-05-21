# ==============================================================================
# COPULA GENERATOR FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# GUMBEL COPULA
# ------------------------------------------------------------------------------

psi_gumbel <- function(t, theta) {

  exp(-(t^(1 / theta)))
}

psi_inv_gumbel <- function(u, theta) {
  (-log(u))^theta
}

psi_prime_gumbel <- function(t, theta) {
  -(1 / theta) * t^(1 / theta - 1) * exp(-t^(1 / theta))
}

psi_inv_prime_gumbel <- function(u, theta) {
  -theta * (-log(u))^(theta - 1) / u
}

# ------------------------------------------------------------------------------
# JOE COPULA
# ------------------------------------------------------------------------------

psi_joe <- function(t, theta) {
  1 - (1 - exp(-t))^(1 / theta)
}

psi_inv_joe <- function(u, theta) {
  -log(1 - (1 - u)^theta)
}

psi_prime_joe <- function(t, theta) {
  -(1 / theta) * (1 - exp(-t))^(1 / theta - 1) * exp(-t)
}

psi_inv_prime_joe <- function(u, theta) {
  -theta * (1 - u)^(theta - 1) / (1 - (1 - u)^theta)
}

# ------------------------------------------------------------------------------
# DISPATCHER FUNCTIONS
# ------------------------------------------------------------------------------

get_psi <- function(family) {
  switch(family,
    gumbel = psi_gumbel,
    joe = psi_joe,
    stop("Unknown copula family: ", family)
  )
}

get_psi_inv <- function(family) {
  switch(family,
    gumbel = psi_inv_gumbel,
    joe = psi_inv_joe,
    stop("Unknown copula family: ", family)
  )
}

get_psi_prime <- function(family) {
  switch(family,
    gumbel = psi_prime_gumbel,
    joe = psi_prime_joe,
    stop("Unknown copula family: ", family)
  )
}

get_psi_inv_prime <- function(family) {
  switch(family,
    gumbel = psi_inv_prime_gumbel,
    joe = psi_inv_prime_joe,
    stop("Unknown copula family: ", family)
  )
}
