# ==============================================================================
# GEV DISTRIBUTION FUNCTIONS
# ==============================================================================

#' GEV cumulative distribution function
#'
#' @param x Numeric vector of quantiles
#' @param mu Location parameter
#' @param sigma Scale parameter (must be > 0)
#' @param xi Shape parameter
#' @return Vector of probabilities
gev_cdf <- function(x, mu, sigma, xi) {

  if (sigma <= 0) {
    return(rep(NA_real_, length(x)))
  }

  z <- (x - mu) / sigma

  # Gumbel limit (xi -> 0)
  if (abs(xi) < 1e-8) {
    return(exp(-exp(-z)))
  }

  # General case

  support <- 1 + xi * z
  out <- rep(0, length(x))
  valid <- support > 0

  out[valid] <- exp(-support[valid]^(-1 / xi))
  out
}

#' GEV log probability density function
#'
#' @param x Numeric vector of quantiles
#' @param mu Location parameter
#' @param sigma Scale parameter (must be > 0)
#' @param xi Shape parameter
#' @return Vector of log-densities
gev_logpdf <- function(x, mu, sigma, xi) {

  if (sigma <= 0) {
    return(rep(-Inf, length(x)))
  }

  z <- (x - mu) / sigma

  # Gumbel limit
  if (abs(xi) < 1e-8) {
    t <- exp(-z)
    return(-log(sigma) - z - t)
  }

  # General case
  support <- 1 + xi * z
  out <- rep(-Inf, length(x))
  valid <- support > 0

  t <- support[valid]^(-1 / xi)
  out[valid] <- -log(sigma) - (1 / xi + 1) * log(support[valid]) - t

  out
}

#' GEV probability density function
#'
#' @param x Numeric vector of quantiles
#' @param mu Location parameter
#' @param sigma Scale parameter (must be > 0)
#' @param xi Shape parameter
#' @return Vector of densities
gev_pdf <- function(x, mu, sigma, xi) {
  exp(gev_logpdf(x, mu, sigma, xi))
}

#' GEV quantile function
#'
#' @param p Numeric vector of probabilities
#' @param mu Location parameter
#' @param sigma Scale parameter (must be > 0)
#' @param xi Shape parameter
#' @return Vector of quantiles
gev_quantile <- function(p, mu, sigma, xi) {

  if (sigma <= 0) {
    return(rep(NA_real_, length(p)))
  }

  # Gumbel limit
  if (abs(xi) < 1e-8) {
    return(mu - sigma * log(-log(p)))
  }

  # General case
  mu + (sigma / xi) * ((-log(p))^(-xi) - 1)
}
