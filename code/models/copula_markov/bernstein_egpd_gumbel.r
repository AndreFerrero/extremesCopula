# Gumbel Copula Log-Density
# u, v are the PIT values from the marginal CDF
log_dgumbel_copula <- function(u, v, theta) {
  # Ensure u, v are in (0, 1)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)

  x <- -log(u)
  y <- -log(v)
  t <- x^theta + y^theta

  # Archimedean formula for Gumbel
  log_c <- -t^(1 / theta) + (theta - 1) * log(x) + (theta - 1) * log(y) +
    (1 / theta - 2) * log(t) + log(t^(1 / theta) + theta - 1) -
    (log(u) + log(v))
  return(log_c)
}

bernstein_copula_nll <- function(theta_vec, x, m) {
  # 1. Parameter Extraction (un-transforming to ensure positivity)
  sigma <- exp(theta_vec[1])
  xi <- theta_vec[2]
  kappa <- exp(theta_vec[3])
  # Bernstein log-weights (softmax approach as used in the package)
  alpha <- theta_vec[4:(3 + m)]
  alpha_shifted <- alpha - max(alpha)
  weights <- exp(alpha_shifted) / sum(exp(alpha_shifted))

  # Gumbel parameter (theta > 1)
  theta_copula <- exp(theta_vec[4 + m]) + 1

  # 2. Marginal Calculations (The Berstein EGPD)
  log_f <- log(egpd:::.bernstein_full_density(x, sigma, xi, kappa, weights, m))

  # 3. Copula Calculations (The Dependence part)
  # Get the marginal CDF values (U_t)
  U_pit <- egpd:::.bernstein_full_cdf(x, sigma, xi, kappa, weights, m)

  # Calculate log copula density for pairs (t, t-1)
  n <- length(x)
  log_c <- log_dgumbel_copula(U_pit[2:n], U_pit[1:(n - 1)], theta_copula)

  # 4. Total NLL
  total_ll <- sum(log_f) + sum(log_c)

  if (!is.finite(total_ll)) {
    return(1e+20)
  }
  return(-total_ll)
}

get_init_bern_egpd_gumbel <- function(x, m) {
  # Marginal fit for initialization
  init_marg <- egpd::fitegpd(x, method = "mle", model = 1)

  # Copula Theta: Use Kendall's Tau between x_t and x_{t-1}
  # theta = 1 / (1 - tau)
  n <- length(x)
  tau <- cor(x[1:(n - 1)], x[2:n], method = "kendall", use = "complete.obs")

  # Safety bounds for Tau (Gumbel needs positive dependence)
  tau_safe <- max(0.01, min(0.99, tau))
  theta_tau <- 1 / (1 - tau_safe)

  init_par <- c(
    log(init_marg$estimate["sigma"]),
    init_marg$estimate["xi"],
    log(init_marg$estimate["kappa"]),
    rep(0, m),
    log(theta_tau - 1)
  )

  return(init_par)
}

fit_egpd_bernstein_gumbel <- function(x, init_par = NULL, m = 5) {
  # Automatically initialize if no parameters provided
  # provide parameters in Log scale

  if (is.null(init_par)) {
    init_par <- get_init_bern_egpd_gumbel(x, m)
  } else {
    if (length(init_par) != (4 + m) || any(is.na(init_par))) {
      stop("Provided init_par must be a numeric vector of length 4 + m with no NAs.")
    }
  }

  opt <- optim(
    par = init_par,
    fn = bernstein_copula_nll,
    x = x,
    m = m,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )

  estimate <- c(
    exp(opt$par[2]),
    opt$par[3],
    exp(opt$par[1]),
    theta = exp(opt$par[4 + m]) + 1
  )

  alpha_raw <- opt$par[4:(3 + m)]
  alpha_shifted <- alpha_raw - max(alpha_raw)
  weights <- exp(alpha_shifted) / sum(exp(alpha_shifted))
  
  return(list(estimate = estimate, weights = weights, opt = opt))
}
