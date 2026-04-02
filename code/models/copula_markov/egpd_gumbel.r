# Gumbel Copula Log-Density
# u, v are the PIT values from the marginal CDF
log_dgumbel_copula <- function(u, v, theta) {
  # Ensure u, v are in (0, 1)
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  v <- pmin(pmax(v, 1e-6), 1 - 1e-6)

  x <- -log(u)
  y <- -log(v)
  t <- x^theta + y^theta

  # Archimedean formula for Gumbel
  log_c <- -t^(1 / theta) + (theta - 1) * log(x) + (theta - 1) * log(y) +
    (1 / theta - 2) * log(t) + log(t^(1 / theta) + theta - 1) -
    (log(u) + log(v))
  return(log_c)
}

# Parametric EGPD + Gumbel Copula NLL
# theta_vec: [log_kappa, log_sigma, xi, log_theta_minus_1]
egpd_gumbel_nll <- function(theta_vec, x) {
  # 1. Parameter Extraction
  kappa <- exp(theta_vec[1])
  sigma <- exp(theta_vec[2])
  xi <- theta_vec[3]
  theta_c <- exp(theta_vec[4]) + 1

  # 2. Marginal Density (EGPD Power Type 1)
  log_f <- egpd:::degpd_density(
    x,
    sigma = sigma,
    xi = xi, kappa = kappa,
    type = 1, log = TRUE
  )
  if (any(1 + xi * x / sigma <= 0)) {
    return(1e20)
  }

  # 3. Copula Dependence
  # U_t = H(x)^kappa
  x_u <- egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1)
  x_u <- pmin(pmax(x_u, 1e-10), 1 - 1e-10)

  n <- length(x)
  # log_dgumbel_copula from previous code block
  log_c <- log_dgumbel_copula(x_u[2:n], x_u[1:(n - 1)], theta_c)

  total_ll <- sum(log_f) + sum(log_c)

  if (!is.finite(total_ll)) {
    return(1e+20)
  }
  return(-total_ll)
}

fit_egpd_gumbel <- function(x, init_par) {
  optim(
    par = init_par,
    fn = egpd_gumbel_nll,
    x = x,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
}
