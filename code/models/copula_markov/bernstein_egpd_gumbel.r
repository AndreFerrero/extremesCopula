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
  log_c <- -t^(1/theta) + (theta-1)*log(x) + (theta-1)*log(y) + 
           (1/theta - 2)*log(t) + log(t^(1/theta) + theta - 1) - 
           (log(u) + log(v))
  return(log_c)
}

bernstein_copula_nll <- function(theta_vec, x, m) {
  # 1. Parameter Extraction (un-transforming to ensure positivity)
  kappa <- exp(theta_vec[1])
  sigma <- exp(theta_vec[2])
  xi    <- theta_vec[3]
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
  log_c <- log_dgumbel_copula(U_pit[2:n], U_pit[1:(n-1)], theta_copula)
  
  # 4. Total NLL
  total_ll <- sum(log_f) + sum(log_c)
  
  if (!is.finite(total_ll)) return(1e+20)
  return(-total_ll)
}

fitegpd_berstein_gumbel <- function(x, init_par = NULL, m = 5) {
  # Automatically initialize if no parameters provided
  if (is.null(init_par)) {
    message("No initial parameters provided. Calculating automatic starts...")
    init_marg <- egpd::fitegpd(x, method = "bernstein", bernstein.m = m)

    init_par <- c(
      log(init_marg["kappa"]),
      log(init_marg["sigma"]),
      init_marg["xi"],
      rep(0, m), # flat weights to start
      log(2) # starting theta_copula at 3 (exp(log(2))+1)
    )
  }

  optim(
    par = init_par,
    fn = bernstein_copula_nll,
    x = x,
    m = m,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )
}