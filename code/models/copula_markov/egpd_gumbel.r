# 1. Gumbel Copula Log-Density (Unchanged)
log_dgumbel_copula <- function(u, v, theta) {
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  v <- pmin(pmax(v, 1e-6), 1 - 1e-6)

  x <- -log(u)
  y <- -log(v)
  t <- x^theta + y^theta

  log_c <- -t^(1 / theta) + (theta - 1) * log(x) + (theta - 1) * log(y) +
    (1 / theta - 2) * log(t) + log(t^(1 / theta) + theta - 1) -
    (log(u) + log(v))
  log_c
}

# 2. Negative Log-Likelihood
# par_trans contains transformed parameters for unconstrained optimization:
# [log(sigma), xi, log(kappa), log(theta - 1)]
egpd_gumbel_nll <- function(par_trans, x) {
  
  # Back-transform for calculation
  sigma <- exp(par_trans[1])
  xi    <- par_trans[2]
  kappa <- exp(par_trans[3])
  theta <- exp(par_trans[4]) + 1

  # Safety check for xi/sigma support
  if (any(1 + xi * x / sigma <= 0)) return(1e20)

  # Marginal density
  log_f <- egpd:::degpd_density(
    x,
    sigma = sigma,
    xi = xi,
    kappa = kappa,
    type = 1,
    log = TRUE
  )

  # PIT values
  u <- egpd:::pegpd(
    x,
    sigma = sigma,
    xi = xi,
    kappa = kappa,
    type = 1
  )
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)

  n <- length(x)
  # Copula density for the dependency
  log_c <- log_dgumbel_copula(
    u[2:n],
    u[1:(n - 1)],
    theta
  )

  ll <- sum(log_f) + sum(log_c)
  if (!is.finite(ll)) return(1e20)

  -ll
}

# 3. Initialization Function
get_init_egpd_gumbel <- function(x) {
  # Fit marginal first
  init <- egpd::fitegpd(x, type = 1, family = "egpd")
  
  # Kendall's tau for copula parameter theta
  n <- length(x)
  tau <- cor(x[-n], x[-1], method = "kendall", use = "complete.obs")
  tau <- max(0.01, min(0.9, tau))
  theta_init <- 1 / (1 - tau)

  # Return in order: sigma, xi, kappa, theta
  c(
    sigma = as.numeric(init$estimate["sigma"]),
    xi    = as.numeric(init$estimate["xi"]),
    kappa = as.numeric(init$estimate["kappa"]),
    theta = theta_init
  )
}

# 4. Fitting Function
fit_egpd_gumbel <- function(x, init_par = NULL, method = "L-BFGS-B", hessian = TRUE) {

  if (is.null(init_par)) {
    init_par <- get_init_egpd_gumbel(x)
  }

  # Transform parameters for optim (to ensure sigma > 0, kappa > 0, theta > 1)
  # Order: 1:sigma, 2:xi, 3:kappa, 4:theta
  init_trans <- c(
    log(init_par[1]),      # sigma
    init_par[2],           # xi
    log(init_par[3]),      # kappa
    log(init_par[4] - 1)   # theta
  )

  opt <- optim(
    par = init_trans,
    fn = egpd_gumbel_nll,
    x = x,
    method = method,
    hessian = hessian,
    control = list(maxit = 10000)
  )

  # ---- Back-transform estimates ----
  theta_hat <- opt$par
  
  final_sigma <- exp(theta_hat[1])
  final_xi    <- theta_hat[2]
  final_kappa <- exp(theta_hat[3])
  final_theta <- exp(theta_hat[4]) + 1

  estimate <- c(
    sigma = final_sigma,
    xi    = final_xi,
    kappa = final_kappa,
    theta = final_theta
  )

  # ---- Variance estimation (Delta Method) ----
  se <- rep(NA_real_, 4)
  names(se) <- names(estimate)
  V_par <- NULL

  if (hessian && !is.null(opt$hessian)) {
    V_theta <- tryCatch(solve(opt$hessian), error = function(e) NULL)

    if (!is.null(V_theta)) {
      # Jacobian Matrix (Partial derivatives of back-transformations)
      # d(exp(h1)) = sigma
      # d(h2)      = 1
      # d(exp(h3)) = kappa
      # d(exp(h4)+1) = theta - 1
      J <- diag(c(
        final_sigma,
        1,
        final_kappa,
        final_theta - 1
      ))

      V_par <- J %*% V_theta %*% t(J)
      rownames(V_par) <- colnames(V_par) <- names(estimate)
      se <- sqrt(pmax(diag(V_par), 0))
    }
  }

  loglik <- -opt$value
  n <- length(x)
  npar <- length(estimate)

  return(structure(list(
    estimate = estimate,
    sd = se,
    vcov = V_par,
    loglik = loglik,
    aic = -2 * loglik + 2 * npar,
    bic = -2 * loglik + log(n) * npar,
    n = n,
    npar = npar,
    convergence = opt$convergence,
    optim = opt,
    call = match.call()
  ), class = "fit_egpd_gumbel"))
}