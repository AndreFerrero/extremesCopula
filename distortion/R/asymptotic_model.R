# ==============================================================================
# ASYMPTOTIC DISTORTED LIMIT MODEL
# ==============================================================================

source("R/copula_generators.R")
source("R/gev_functions.R")

#' Asymptotic distorted CDF
#'
#' @param x Numeric vector of quantiles
#' @param mu GEV location
#' @param sigma GEV scale
#' @param xi GEV shape
#' @param theta Copula dependence parameter
#' @param family Copula family ("gumbel" or "joe")
#' @return Vector of probabilities
G_asymptotic <- function(x, mu, sigma, xi, theta, family) {

  H <- gev_cdf(x, mu, sigma, xi)
  H <- pmin(pmax(H, 1e-12), 1 - 1e-12)

  V <- -log(H)
  rho <- 1 / theta

  psi <- get_psi(family)
  psi(V^(1 / rho), theta)
}

#' Asymptotic distorted density
#'
#' @param x Numeric vector of quantiles
#' @param mu GEV location
#' @param sigma GEV scale
#' @param xi GEV shape
#' @param theta Copula dependence parameter
#' @param family Copula family ("gumbel" or "joe")
#' @return Vector of densities
g_asymptotic <- function(x, mu, sigma, xi, theta, family) {

  H <- gev_cdf(x, mu, sigma, xi)
  H <- pmin(pmax(H, 1e-12), 1 - 1e-12)

  h <- gev_pdf(x, mu, sigma, xi)
  V <- -log(H)

  # Gumbel: no distortion in density

if (family == "gumbel") {
    return(h)
  }

  # Joe (and other Archimedean families)
  rho <- 1 / theta
  z <- V^(1 / rho)

  psi_prime <- get_psi_prime(family)
  psi_prime_val <- psi_prime(z, theta)

  dens <- -psi_prime_val * (1 / rho) * V^(1 / rho - 1) * h / H
  dens
}

#' Asymptotic distorted quantile function
#'
#' @param p Numeric vector of probabilities
#' @param mu GEV location
#' @param sigma GEV scale
#' @param xi GEV shape
#' @param theta Copula dependence parameter
#' @param family Copula family ("gumbel" or "joe")
#' @return Vector of quantiles
model_quantile <- function(p, mu, sigma, xi, theta, family) {

  eps <- 1e-10
  p <- pmin(pmax(p, eps), 1 - eps)

  if (family == "gumbel") {
    return(gev_quantile(p, mu, sigma, xi))
  }

  rho <- 1 / theta
  psi_inv <- get_psi_inv(family)
  t <- psi_inv(p, theta)

  gev_quantile(exp(-t^rho), mu, sigma, xi)
}
