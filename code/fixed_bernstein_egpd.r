bernstein_nll_fixed <- function(theta, x, m, fix.arg = NULL) {
  fixed_names <- names(fix.arg)

  pos <- 1

  if (!("kappa" %in% fixed_names)) {
    kappa <- exp(theta[pos])
    pos <- pos + 1
  } else {
    kappa <- fix.arg$kappa
  }

  if (!("sigma" %in% fixed_names)) {
    sigma <- exp(theta[pos])
    pos <- pos + 1
  } else {
    sigma <- fix.arg$sigma
  }

  if (!("xi" %in% fixed_names)) {
    xi <- theta[pos]
    pos <- pos + 1
  } else {
    xi <- fix.arg$xi
  }

  alpha <- theta[pos:(pos + m - 1)]

  alpha_shifted <- alpha - max(alpha)
  weights <- exp(alpha_shifted)
  weights <- weights / sum(weights)

  u <- egpd:::.pgpd_std(x, sigma, xi)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)

  uk <- u^kappa
  uk <- pmin(pmax(uk, 1e-10), 1 - 1e-10)

  b_dens <- egpd:::.bernstein_density(uk, weights, m)
  b_dens <- pmax(b_dens, .Machine$double.xmin)

  gpd_logd <- egpd:::.dgpd_std(x, sigma, xi, log = TRUE)

  ll <- sum(
    log(b_dens) +
      log(kappa) +
      (kappa - 1) * log(u) +
      gpd_logd
  )

  if (!is.finite(ll)) {
    return(1e20)
  }

  -ll
}

fitegpd_bernstein_fixed <- function(
  x,
  type = 1,
  m = 8,
  start = NULL,
  fix.arg = NULL,
  optim.method = "Nelder-Mead",
  hessian = TRUE,
  ...
) {
  n <- length(x)

  fixed_names <- names(fix.arg)

  #################################################
  ## starting values
  #################################################

  xpos <- x[x > 0]

  if (length(xpos) < 2) {
    xpos <- abs(x) + 0.01
  }

  mx <- mean(xpos)
  vx <- var(xpos)

  xi0 <- max(min(0.5 * (mx^2 / vx - 1), 0.5), -0.5)
  sigma0 <- max(mx * (1 - xi0), 0.01)
  kappa0 <- 1

  if (!is.null(start)) {
    if ("sigma" %in% names(start)) {
      sigma0 <- start$sigma
    }

    if ("xi" %in% names(start)) {
      xi0 <- start$xi
    }

    if ("kappa" %in% names(start)) {
      kappa0 <- start$kappa
    }
  }

  #################################################
  ## build free parameter vector
  #################################################

  theta0 <- numeric()

  if (!("kappa" %in% fixed_names)) {
    theta0 <- c(theta0, log(kappa0))
  }

  if (!("sigma" %in% fixed_names)) {
    theta0 <- c(theta0, log(sigma0))
  }

  if (!("xi" %in% fixed_names)) {
    theta0 <- c(theta0, xi0)
  }

  theta0 <- c(theta0, rep(0, m))

  #################################################
  ## optimize
  #################################################

  opt <- optim(
    par = theta0,
    fn = bernstein_nll_fixed,
    x = x,
    m = m,
    fix.arg = fix.arg,
    method = optim.method,
    hessian = hessian,
    control = list(maxit = 10000),
    ...
  )

  #################################################
  ## reconstruct estimates
  #################################################

  pos <- 1

  if (!("kappa" %in% fixed_names)) {
    kappa_hat <- exp(opt$par[pos])
    pos <- pos + 1
  } else {
    kappa_hat <- fix.arg$kappa
  }

  if (!("sigma" %in% fixed_names)) {
    sigma_hat <- exp(opt$par[pos])
    pos <- pos + 1
  } else {
    sigma_hat <- fix.arg$sigma
  }

  if (!("xi" %in% fixed_names)) {
    xi_hat <- opt$par[pos]
    pos <- pos + 1
  } else {
    xi_hat <- fix.arg$xi
  }

  alpha_hat <- opt$par[pos:(pos + m - 1)]

  alpha_shifted <- alpha_hat - max(alpha_hat)

  weights <- exp(alpha_shifted)
  weights <- weights / sum(weights)

  estimate <- c(
    sigma = sigma_hat,
    xi = xi_hat,
    kappa = kappa_hat
  )

  #################################################
  ## parameter counting
  #################################################

  n_free_main <- 3 - length(intersect(
    c("sigma", "xi", "kappa"),
    fixed_names
  ))

  npar <- n_free_main + (m - 1)

  ll <- -opt$value

  #################################################
  ## output
  #################################################

  structure(
    list(
      estimate = estimate,
      loglik = ll,
      aic = -2 * ll + 2 * npar,
      bic = -2 * ll + log(n) * npar,
      n = n,
      npar = npar,
      family = "egpd",
      method = "bernstein",
      convergence = opt$convergence,
      optim = opt,
      bernstein.m = m,
      bernstein.weights = weights,
      fix.arg = fix.arg
    ),
    class = "fitegpd"
  )
}
