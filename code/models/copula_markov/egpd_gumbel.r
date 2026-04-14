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

get_init_egpd_gumbel <- function(x) {
  # Marginal fit for initialization
  init_marg <- egpd::fitegpd(x, method = "mle", model = 1)

  # Copula Theta: Use Kendall's Tau between x_t and x_{t-1}
  # theta = 1 / (1 - tau)
  n <- length(x)
  tau <- cor(x[1:(n - 1)], x[2:n], method = "kendall", use = "complete.obs")

  # Safety bounds for Tau (Gumbel needs positive dependence)
  tau_safe <- max(0.01, min(0.99, tau))
  theta_tau <- 1 / (1 - tau_safe)

  # 5. Construct the 4-parameter vector explicitly
  # Ensure these are single numeric values, not lists or named vectors
  # Initial theta_vec
  # [log_kappa, log_sigma, xi, log(theta_copula - 1)]
  init_par <- c(
    log(init_marg$estimate["kappa"]),
    log(init_marg$estimate["sigma"]),
    init_marg$estimate["xi"],
    log(theta_tau - 1)
  )

  # Double check length and NAs
  if (length(init_par) != 4 || any(is.na(init_par))) {
    stop("Initialization failed: parameters contain NAs or wrong length.")
  }

  return(init_par)
}

fit_egpd_gumbel <- function(x, init_par = NULL, method = "L-BFGS-B") {
  # Automatically initialize if no parameters provided
  # Provide Parameters in Log scale
  
  if (is.null(init_par)) {
    init_par <- get_init_egpd_gumbel(x)
  } else {
    # Validate provided init_par
    if (length(init_par) != 4 || any(is.na(init_par))) {
      stop("Provided init_par must be a numeric vector of length 4 with no NAs.")
    }
  }

  opt <- optim(
    par = init_par,
    fn = egpd_gumbel_nll,
    x = x,
    method = method,
    control = list(maxit = 10000)
  )

  estimate <- c(
    exp(opt$par[1]),
    exp(opt$par[2]),
    opt$par[3],
    theta = exp(opt$par[4]) + 1
  )

  return(list(estimate = estimate, opt = opt))
}

#' Calculate Marginal and Conditional PIT values
#' @param data The observed vector x
#' @param theta_vec The estimated parameter vector [log_kappa, log_sigma, xi, log_theta_minus_1]
#' @param h_dist_fn The h-function (conditional CDF) from your copula object
get_pit_values <- function(x, theta_vec, h_dist_fn) {
  # 1. Parameter Extraction
  kappa <- exp(theta_vec[1])
  sigma <- exp(theta_vec[2])
  xi <- theta_vec[3]
  theta_c <- exp(theta_vec[4]) + 1

  # 2. Marginal PIT (u_t)
  u <- egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10) # Numerical safety

  # 3. Conditional PIT (w_t)
  n <- length(u)
  w <- h_dist_fn(u[2:n], u[1:(n - 1)], theta_c)

  return(list(u = u, w = w))
}

#' Generate Side-by-Side QQ-plots
#' @param pit_list A list containing 'u' and 'w' from get_pit_values
plot_diag_plots <- function(pit_list) {
  u <- pit_list$u
  w <- pit_list$w

  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par)) # Restores plot settings after function finishes

  # Plot A: Marginal Check
  plot(stats::ppoints(length(u)), sort(u),
    main = "Marginal QQ-plot (EGPD)",
    xlab = "Theoretical Uniform", ylab = "Empirical u_t",
    pch = 20, col = "grey60"
  )
  abline(0, 1, col = "firebrick", lwd = 2)

  # Plot B: Conditional Check
  plot(stats::ppoints(length(w)), sort(w),
    main = "Conditional QQ-plot (Gumbel)",
    xlab = "Theoretical Uniform", ylab = "Empirical w_t",
    pch = 20, col = "royalblue"
  )
  abline(0, 1, col = "firebrick", lwd = 2)
}
