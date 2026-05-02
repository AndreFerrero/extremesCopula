library(egpd)
library(copula)
library(numDeriv)

# --- Internal Utility: Student-t Log-Density ---
log_dt_copula <- function(u, v, rho, nu) {
  # Numerical safety boundaries
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  
  # t-Copula density using the copula package
  log_density <- dCopula(cbind(u, v), tCopula(param = rho, df = nu), log = TRUE)
  return(log_density)
}

# --- Internal Utility: Joint NLL for Student-t ---
egpd_t_nll <- function(par_trans, x) {
  sigma <- exp(par_trans[1])
  xi    <- par_trans[2]
  kappa <- exp(par_trans[3])
  
  # Transformations for t-copula parameters
  rho   <- tanh(par_trans[4])         # Map (-inf, inf) to (-1, 1)
  nu    <- exp(par_trans[5]) + 2.01    # Map (-inf, inf) to (2.01, inf)

  if (any(1 + xi * x / sigma <= 0)) {
    return(1e20)
  }

  # Marginal Density
  log_f <- egpd:::degpd_density(x, sigma = sigma, xi = xi, kappa = kappa, type = 1, log = TRUE)
  
  # Probability Integral Transform (PIT)
  u <- pmin(pmax(egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1), 1e-10), 1 - 1e-10)

  n <- length(x)
  # Copula Density for the Markov chain
  log_c <- log_dt_copula(u[2:n], u[1:(n - 1)], rho, nu)

  ll <- sum(log_f) + sum(log_c)
  
  if (!is.finite(ll)) {
    return(1e20)
  }
  -ll
}

# --- THE UNIFIED STUDENT-T FITTING FUNCTION ---
fit_egpd_t_copula <- function(x, method = c("mle", "iterative", "ifm"),
                              init_par = NULL, hessian = TRUE, optim.method = "Nelder-Mead") {
  method <- match.arg(method)
  n <- length(x)

  # 1. Initialization
  if (is.null(init_par)) {
    # Marginal starts
    init_marg <- egpd::fitegpd(x, type = 1, family = "egpd")$estimate
    
    # Copula starts (Pearson correlation for rho, 4 for nu as a standard heavy-tail start)
    rho_init <- cor(x[-n], x[-1], method = "pearson")
    init_par <- c(init_marg, rho = rho_init, nu = 4)
  }

  # 2. Optimization Paths
  if (method == "mle") {
    # Transform: log(sigma), xi, log(kappa), atanh(rho), log(nu-2.01)
    init_trans <- c(log(init_par[1]), init_par[2], log(init_par[3]), 
                    atanh(init_par[4]), log(init_par[5] - 2.01))
    
    opt <- optim(init_trans, egpd_t_nll, x = x, method = optim.method, 
                 control = list(maxit = 3000, parscale = c(1, 0.01, 1, 1, 1)))
    final_trans <- opt$par
    conv <- opt$convergence
    
  } else if (method == "iterative") {
    cur_t <- c(log(init_par[1]), init_par[2], log(init_par[3]), 
               atanh(init_par[4]), log(init_par[5] - 2.01))
    
    for (i in 1:50) {
      # Step A: Margins (Fix rho, nu)
      opt_m <- optim(cur_t[1:3], function(p) egpd_t_nll(c(p, cur_t[4:5]), x), method = optim.method)
      # Step B: Copula (Fix margins)
      opt_c <- optim(cur_t[4:5], function(p) egpd_t_nll(c(opt_m$par, p), x), method = optim.method)
      
      new_t <- c(opt_m$par, opt_c$par)
      if (sum((new_t - cur_t)^2) < 1e-6) break
      cur_t <- new_t
    }
    final_trans <- cur_t
    conv <- 0
    
  } else if (method == "ifm") {
    # Stage 1: IID Margins
    fit_m <- egpd::fitegpd(x, type = 1, family = "egpd")
    m_p <- fit_m$estimate
    
    # Stage 2: Copula only
    u <- egpd:::pegpd(x, sigma = m_p[1], xi = m_p[2], kappa = m_p[3], type = 1)
    u_pairs <- cbind(u[-n], u[-1])
    
    # Fit t-copula using pseudo-observations
    fit_c <- fitCopula(tCopula(dim = 2), pobs(u_pairs), method = "ml")
    
    final_trans <- c(log(m_p[1]), m_p[2], log(m_p[3]), 
                     atanh(fit_c@estimate[1]), log(fit_c@estimate[2] - 2.01))
    conv <- 0
  }

  # 3. Back-transform Results
  est <- c(
    exp(final_trans[1]), 
    final_trans[2],
    exp(final_trans[3]), 
    tanh(final_trans[4]),
    exp(final_trans[5]) + 2.01
  )

  loglik <- -egpd_t_nll(final_trans, x)

  # Standard Errors via Hessian
  se <- rep(NA, 5)
  vcov_mat <- NULL
  if (hessian) {
    h <- tryCatch(numDeriv::hessian(function(p) egpd_t_nll(p, x), final_trans), error = function(e) NULL)
    if (!is.null(h)) {
      v_t <- tryCatch(solve(h), error = function(e) NULL)
      if (!is.null(v_t)) {
        # Jacobian for the back-transformation
        jac <- diag(c(est[1], 1, est[3], 1 - est[4]^2, est[5] - 2.01))
        vcov_mat <- jac %*% v_t %*% t(jac)
        se <- sqrt(pmax(diag(vcov_mat), 0))
      }
    }
  }

  structure(list(
    estimate = est, sd = se, vcov = vcov_mat, loglik = loglik,
    aic = -2 * loglik + 10, method = method, n = n
  ), class = "fit_egpd_copula")
}