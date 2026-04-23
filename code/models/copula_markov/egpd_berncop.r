# 1. Bernstein Copula Density using Beta PDFs
log_dbernstein_copula <- function(u, v, alpha_vec, m) {
  # Standard clamping for safety
  u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
  v <- pmin(pmax(v, 1e-8), 1 - 1e-8)
  
  # Map alpha_vec to a symmetric weight matrix (m+1 x m+1)
  W <- matrix(0, m + 1, m + 1)
  W[lower.tri(W, diag = TRUE)] <- exp(alpha_vec)
  W <- W + t(W) - diag(diag(W))

  # Normalize so weights sum to 1 (Mixture weights)
  W <- W / sum(W)
  
  # Create Beta density matrices: rows are observations, columns are mixture components
  # Degree m involves components 1 to m+1
  bu <- sapply(1:(m + 1), function(i) dbeta(u, shape1 = i, shape2 = m - i + 2))
  bv <- sapply(1:(m + 1), function(j) dbeta(v, shape1 = j, shape2 = m - j + 2))
  
  # Mixture Density: sum_{i,j} w_ij * Beta_i(u) * Beta_j(v)
  # Efficiently: diag( bu %*% W %*% t(bv) )
  dens <- rowSums((bu %*% W) * bv)
  
  log(pmax(dens, 1e-12))
}

get_init_egpd_bernstein <- function(x, m) {
  # Fit marginal first
  init <- egpd::fitegpd(x, type = 1, family = "egpd")

  # Initial alpha_vec for independence (all weights equal)
  num_cop_par <- (m + 1) * (m + 2) / 2
  alpha_init <- rep(0, num_cop_par) # exp(0) = 1, leads to uniform weights

  # Return in order: sigma, xi, kappa, alpha_vec
  c(
    as.numeric(init$estimate["sigma"]),
    as.numeric(init$estimate["xi"]),
    as.numeric(init$estimate["kappa"]),
    alpha_init
  )
}

# 2. stabilized NLL
egpd_bernstein_nll <- function(par_trans, x, m) {
  sigma <- exp(par_trans[1])
  xi    <- par_trans[2]
  kappa <- exp(par_trans[3])
  alpha_vec <- par_trans[4:length(par_trans)]

  # Support constraint: x must be within (0, sigma/(-xi)) if xi < 0
  if (any(1 + xi * x / sigma <= 1e-7)) return(1e15)

  # Marginal Density (Protected)
  log_f <- try(egpd:::degpd_density(x,
    sigma = sigma,
    xi = xi,
    kappa = kappa,
    type = 1,
    log = TRUE),
    silent = TRUE)
  if (inherits(log_f, "try-error") || any(!is.finite(log_f))) return(1e15)

  # PIT values (Protected)
  u <- try(egpd:::pegpd(x,
    sigma = sigma,
    xi = xi,
    kappa = kappa,
    type = 1),
    silent = TRUE)
  if (inherits(u, "try-error") || any(!is.finite(u))) return(1e15)
  u <- pmin(pmax(u, 1e-9), 1 - 1e-9)

  # Copula Density
  n <- length(x)
  log_c <- log_dbernstein_copula(u[2:n], u[1:(n-1)], alpha_vec, m = m)
  
  val <- -(sum(log_f) + sum(log_c))
  if (!is.finite(val) || is.na(val)) return(1e15)
  val
}

# 3. Robust Fitting Function
fit_egpd_bernstein <- function(x, init_par = NULL, m = 3, method = "Nelder-Mead", hessian = TRUE) {
  if (is.null(init_par)) {
    init_par <- get_init_egpd_bernstein(x, m)
  }
  
  init_trans <- c(
    log(init_par[1]),      # sigma
    init_par[2],           # xi
    log(init_par[3]),      # kappa
    init_par[4:length(init_par)]
  )

  # Final safety check before starting
  test_ll <- egpd_bernstein_nll(init_trans, x, m)
  if(test_ll > 1e14) stop("Initial NLL is non-finite.")
  message(paste("Starting NLL:", round(test_ll, 2)))
  
  # Step 3: Joint Optimization
  opt <- optim(
    par = init_trans,
    fn = egpd_bernstein_nll,
    x = x, m = m,
    method = method,
    hessian = hessian,
    control = list(maxit = 5000)
  )

  # ---- Back-transform estimates ----
  theta_hat <- opt$par
  
  final_sigma <- exp(theta_hat[1])
  final_xi    <- theta_hat[2]
  final_kappa <- exp(theta_hat[3])
  final_alphas <- theta_hat[4:length(theta_hat)]

  estimate <- c(
    sigma = final_sigma,
    xi    = final_xi,
    kappa = final_kappa
  )
  # Name the alpha parameters
  names(final_alphas) <- paste0("alpha", 1:length(final_alphas))
  estimate <- c(estimate, final_alphas)

  # ---- Variance estimation (Delta Method) ----
  npar <- length(estimate)
  se <- rep(NA_real_, npar)
  names(se) <- names(estimate)
  V_par <- NULL

  if (hessian && !is.null(opt$hessian)) {
    # We use Moore-Penrose inverse (ginv) or tryCatch for stability 
    # since Bernstein parameters can occasionally cause singular Hessians
    V_theta <- tryCatch(solve(opt$hessian), error = function(e) NULL)

    if (!is.null(V_theta)) {
      # Jacobian Matrix (Partial derivatives of back-transformations)
      # d(exp(h1))   = sigma
      # d(h2)        = 1
      # d(exp(h3))   = kappa
      # d(h4...hp)   = 1 (for all alphas)
      
      diag_jacobian <- c(
        final_sigma,
        1,
        final_kappa,
        rep(1, length(final_alphas))
      )
      
      J <- diag(diag_jacobian)

      V_par <- J %*% V_theta %*% t(J)
      rownames(V_par) <- colnames(V_par) <- names(estimate)
      
      # Ensure non-negative variances
      se <- sqrt(pmax(diag(V_par), 0))
    }
  }

  loglik <- -opt$value
  n <- length(x)

  # Return structured object
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
    m = m, # Store degree used
    call = match.call()
  ), class = "fit_egpd_bernstein"))
}

# 1. Updated NLL with Weighted Likelihood and Regularization
egpd_bernstein_nll_robust <- function(par_trans, x, m, weight_copula = 1.0) {
  sigma <- exp(par_trans[1])
  xi    <- par_trans[2] # Will be constrained by L-BFGS-B
  kappa <- exp(par_trans[3])
  alpha_vec <- par_trans[4:length(par_trans)]

  # 1. Marginal Part
  log_f <- try(egpd:::degpd_density(x, sigma = sigma, xi = xi, kappa = kappa, type = 1, log = TRUE), silent = TRUE)
  if (inherits(log_f, "try-error") || any(!is.finite(log_f))) return(1e15)

  # 2. Copula Part
  u <- try(egpd:::pegpd(x, sigma=sigma, xi=xi, kappa=kappa, type = 1), silent = TRUE)
  if (inherits(u, "try-error")) return(1e15)
  u <- pmin(pmax(u, 1e-9), 1 - 1e-9)

  n <- length(x)
  log_c <- log_dbernstein_copula(u[2:n], u[1:(n - 1)], alpha_vec, m = m)
  
  # Penalty to prevent Alpha Explosion (Tikhonov Regularization)
  # This keeps weights from becoming "infinite spikes"
  penalty <- 0.01 * sum(alpha_vec^2)

  # Total Weighted Likelihood
  # We divide the copula part by (n-1) and multiply by n to put them on the same scale
  val <- -sum(log_f) - (weight_copula * sum(log_c)) + penalty
  
  if (!is.finite(val)) return(1e15)
  val
}

# 2. Updated Fitting Function with Box Constraints
fit_egpd_bernstein_robust <- function(x, m = 2, weight_copula = 1.0) {
  init_m <- egpd::fitegpd(x, type = 1, family = "egpd")
  num_cop_par <- (m + 1) * (m + 2) / 2
  
  init_trans <- c(
    log(as.numeric(init_m$estimate["sigma"])),
    as.numeric(init_m$estimate["xi"]),
    log(as.numeric(init_m$estimate["kappa"])),
    rep(0, num_cop_par)
  )

  # Use L-BFGS-B to enforce physical boundaries
  # xi should stay between -0.25 and 0.25 for wind
  # kappa should stay > 1 to preserve the mode (the "bump")
  lower_bounds <- c(-Inf, -0.25, log(1.1), rep(-10, num_cop_par))
  upper_bounds <- c(Inf,  0.5,  Inf,      rep( 10, num_cop_par))

  opt <- optim(
    par = init_trans,
    fn = egpd_bernstein_nll_robust,
    x = x, m = m, weight_copula = weight_copula,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(maxit = 5000, pgtol = 1e-6)
  )

  # Back-transform and structure output as before...
  # (Standard transform code here)
  return(opt)
}
