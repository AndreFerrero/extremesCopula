# 1. Joe Copula Log-Density
log_djoe_copula <- function(u, v, theta) {
  # Numerical stability
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  
  alpha <- 1 / theta
  
  # Term components
  u_term <- (1 - u)^theta
  v_term <- (1 - v)^theta
  
  # om_hJ is the "outer" part of the Joe survival function logic
  om_hJ <- u_term + v_term - (u_term * v_term)
  # hJ is the interaction
  hJ <- (1 - u_term) * (1 - v_term)
  
  # Polynomial term from the derivative of the generator
  # This corresponds to the log_poly in your Stan code
  log_poly <- log1p((1 - alpha) * (hJ / om_hJ))
  
  log_c <- log(theta) + (theta - 1) * (log1p(-u) + log1p(-v)) - 
           (1 - alpha) * log(om_hJ) + log_poly
  
  return(log_c)
}

# 2. Negative Log-Likelihood for Joe-EGPD
# par_trans: [log(sigma), xi, log(kappa), log(theta - 1)]
egpd_joe_nll <- function(par_trans, x) {
  
  sigma <- exp(par_trans[1])
  xi    <- par_trans[2]
  kappa <- exp(par_trans[3])
  theta <- exp(par_trans[4]) + 1

  # Support check
  if (any(1 + xi * x / sigma <= 0)) return(1e20)

  # Marginal density using egpd package internals
  log_f <- egpd:::degpd_density(x, sigma = sigma, xi = xi, kappa = kappa, type = 1, log = TRUE)

  # Probability Integral Transform (PIT)
  u <- egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)

  n <- length(x)
  # Joe Copula density for transitions
  log_c <- log_djoe_copula(u[2:n], u[1:(n - 1)], theta)

  ll <- sum(log_f) + sum(log_c)
  if (!is.finite(ll)) return(1e20)

  return(-ll)
}

# 3. Initialization Function
get_init_egpd_joe <- function(x) {
  # Fit marginal first to get bulk/tail parameters
  init <- egpd::fitegpd(x, type = 1, family = "egpd")
  
  # Use Kendall's Tau for initial theta
  # For Joe: tau = 1 + 4/theta^2 * integral... approx 1 - 2/(theta+2)
  # A simple heuristic for Joe starting value
  n <- length(x)
  tau <- cor(x[-n], x[-1], method = "kendall", use = "complete.obs")
  tau <- max(0.01, min(0.9, tau))
  # Heuristic start: theta approx 2*tau / (1-tau)
  theta_init <- max(1.1, (2 * tau) / (1 - tau))

  c(
    as.numeric(init$estimate["sigma"]),
    as.numeric(init$estimate["xi"]),
    as.numeric(init$estimate["kappa"]),
    theta_init
  )
}

# 4. Fitting Function (MLE)
fit_egpd_joe <- function(x, init_par = NULL, method = "L-BFGS-B", hessian = TRUE) {

  if (is.null(init_par)) {
    init_par <- get_init_egpd_joe(x)
  }

  init_trans <- c(
    log(init_par[1]),      # log sigma
    init_par[2],           # xi
    log(init_par[3]),      # log kappa
    log(init_par[4] - 1)   # log(theta - 1)
  )

  # Using L-BFGS-B to keep xi in physical bounds if needed
  # Wind xi usually -0.2 to 0.4. Constraints help prevent the "infinite mean" trap.
  opt <- optim(
    par = init_trans,
    fn = egpd_joe_nll,
    x = x,
    method = method,
    hessian = hessian,
    control = list(maxit = 10000)
  )

  # Back-transform
  final_sigma <- exp(opt$par[1])
  final_xi    <- opt$par[2]
  final_kappa <- exp(opt$par[3])
  final_theta <- exp(opt$par[4]) + 1

  estimate <- c(sigma = final_sigma, xi = final_xi, kappa = final_kappa, theta = final_theta)

  # Variance estimation (Delta Method)
  se <- rep(NA_real_, 4)
  names(se) <- names(estimate)
  V_par <- NULL

  if (hessian && !is.null(opt$hessian)) {
    V_theta <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    if (!is.null(V_theta)) {
      J <- diag(c(final_sigma, 1, final_kappa, final_theta - 1))
      V_par <- J %*% V_theta %*% t(J)
      rownames(V_par) <- colnames(V_par) <- names(estimate)
      se <- sqrt(pmax(diag(V_par), 0))
    }
  }

  loglik <- -opt$value
  return(structure(list(
    estimate = estimate,
    sd = se,
    vcov = V_par,
    loglik = loglik,
    aic = -2 * loglik + 2 * length(estimate),
    bic = -2 * loglik + log(length(x)) * length(estimate),
    convergence = opt$convergence
  ), class = "fit_egpd_joe"))
}