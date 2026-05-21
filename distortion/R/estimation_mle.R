# ==============================================================================
# MAXIMUM LIKELIHOOD ESTIMATION
# ==============================================================================

source("R/asymptotic_model.R")

#' Negative log-likelihood for joint estimation (all parameters)
#'
#' @param par Parameter vector: c(mu, log(sigma), xi, theta)
#' @param M Vector of block maxima
#' @param family Copula family
#' @return Negative log-likelihood value
negloglik_joint <- function(par, M, family) {

  mu <- par[1]
  sigma <- exp(par[2])
  xi <- par[3]
  theta <- par[4]

  if (sigma <= 0 || theta <= 1) return(1e10)

  dens <- g_asymptotic(M, mu, sigma, xi, theta, family)

  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e10)

  -sum(log(dens))
}

#' Joint MLE estimation
#'
#' @param M Vector of block maxima
#' @param family Copula family
#' @param theta_init Initial theta (e.g., from DMLE)
#' @param maxit Maximum iterations
#' @return List with estimates
fit_joint <- function(M, family, theta_init = 2, optim_method = "Nelder-Mead", maxit = 5000) {

  mu0 <- mean(M)
  sigma0 <- sd(M)
  xi0 <- 0.1
  theta0 <- max(theta_init, 1.05)

  par0 <- c(mu0, log(sigma0), xi0, theta0)

  fit <- optim(
    par = par0,
    fn = negloglik_joint,
    M = M,
    family = family,
    method = optim_method,
    control = list(maxit = maxit)
  )

  list(
    mu = fit$par[1],
    sigma = exp(fit$par[2]),
    xi = fit$par[3],
    theta = fit$par[4],
    family = family,
    convergence = fit$convergence,
    value = fit$value
  )
}

#' Negative log-likelihood for second stage (GEV only, theta fixed)
#'
#' @param par Parameter vector: c(mu, log(sigma), xi)
#' @param M Vector of block maxima
#' @param theta_fixed Fixed copula parameter
#' @param family Copula family
#' @return Negative log-likelihood value
negloglik_stage2 <- function(par, M, theta_fixed, family) {

  mu <- par[1]
  sigma <- exp(par[2])
  xi <- par[3]

  if (sigma <= 0) return(1e10)

  dens <- g_asymptotic(M, mu, sigma, xi, theta_fixed, family)

  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e10)

  -sum(log(dens))
}

#' Two-stage estimation: fix theta, estimate GEV parameters
#'
#' @param M Vector of block maxima
#' @param theta_hat Estimated theta from DMLE
#' @param family Copula family
#' @param maxit Maximum iterations
#' @return List with estimates
fit_twostage_gev <- function(M, theta_hat, family, optim_method = "Nelder-Mead", maxit = 5000) {

  mu0 <- mean(M)
  sigma0 <- sd(M)
  xi0 <- 0.1

  par0 <- c(mu0, log(sigma0), xi0)

  fit <- optim(
    par = par0,
    fn = negloglik_stage2,
    M = M,
    theta_fixed = theta_hat,
    family = family,
    method = optim_method,
    control = list(maxit = maxit)
  )

  list(
    mu = fit$par[1],
    sigma = exp(fit$par[2]),
    xi = fit$par[3],
    theta = theta_hat,
    family = family,
    convergence = fit$convergence,
    value = fit$value
  )
}
