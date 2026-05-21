# ==============================================================================
# DIAGONAL MLE ESTIMATION
# ==============================================================================

source("R/copula_generators.R")

#' Diagonal density for Joe copula
#'
#' @param u Vector of pseudo-observations (row-wise maxima)
#' @param theta Copula parameter
#' @param n Dimension (number of columns)
#' @return Vector of densities
d_delta_joe <- function(u, theta, n) {

  t <- psi_inv_joe(u, theta)
  n * psi_prime_joe(n * t, theta) * psi_inv_prime_joe(u, theta)
}

#' Negative log-likelihood for Joe DMLE
#'
#' @param theta Copula parameter
#' @param Y Vector of pseudo-observation maxima
#' @param n Dimension
#' @return Negative log-likelihood value
negloglik_dmle_joe <- function(theta, Y, n) {

  if (theta <= 1) return(1e10)

  dens <- d_delta_joe(Y, theta, n)

  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e10)

  -sum(log(dens))
}

#' DMLE for Joe copula
#'
#' @param Yhat Vector of pseudo-observation maxima
#' @param n Dimension
#' @param init Initial theta value
#' @return Estimated theta
dmle_joe <- function(Yhat, n, init = 2) {

  fit <- optim(
    par = init,
    fn = negloglik_dmle_joe,
    Y = Yhat,
    n = n,
    method = "L-BFGS-B",
    lower = 1.001,
    upper = 20
  )

  fit$par
}

#' DMLE for Gumbel copula (closed form)
#'
#' @param Yhat Vector of pseudo-observation maxima
#' @param n Dimension
#' @param k Number of observations
#' @return Estimated theta
dmle_gumbel <- function(Yhat, n, k) {

  S <- sum(-log(Yhat))
  log(n) / (log(k) - log(S))
}

#' DMLE dispatcher
#'
#' @param Yhat Vector of pseudo-observation maxima
#' @param n Dimension
#' @param k Number of observations
#' @param family Copula family
#' @return Estimated theta
estimate_theta_dmle <- function(Yhat, n, k, family) {

  switch(family,
    joe = dmle_joe(Yhat, n),
    gumbel = dmle_gumbel(Yhat, n, k),
    stop("Unknown copula family: ", family)
  )
}
