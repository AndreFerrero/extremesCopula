library(egpd)
library(copula)
library(numDeriv)

# --- Internal Utility: t-EV Log-Density ---
log_dtev_copula <- function(u, v, rho, df) {
  # Numerical safety boundaries
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)

  # tevCopula in R takes two parameters: 
  # rho (correlation) and df (degrees of freedom)
  log_density <- dCopula(cbind(u, v), tevCopula(param = rho, df = df), log = TRUE)
  
  return(log_density)
}

# --- Internal Utility: Joint NLL for t-EV-EGPD ---
egpd_tev_nll <- function(par_trans, x) {
  # Back-transform parameters
  sigma  <- exp(par_trans[1])
  xi     <- par_trans[2]           # xi is real
  kappa  <- exp(par_trans[3])
  # rho must be in [0, 1]. Use sigmoid (plogis) for stability
  rho    <- plogis(par_trans[4])   
  df     <- exp(par_trans[5])      # df > 0
  
  # Domain check for EGPD
  if (any(1 + xi * x / sigma <= 0)) return(1e20)
  
  # Marginal Log-Likelihood (EGPD)
  log_f <- egpd:::degpd_density(x, sigma = sigma, xi = xi, kappa = kappa, type = 1, log = TRUE)
  
  # Probability Integral Transform (PIT)
  u <- pmin(pmax(egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1), 1e-10), 1 - 1e-10)
  
  # Dependence Log-Likelihood (t-EV)
  n <- length(x)
  # Model x_t and x_{t+1} as a Markov Chain
  log_c <- log_dtev_copula(u[1:(n - 1)], u[2:n], rho, df)
  
  ll <- sum(log_f) + sum(log_c)
  
  if (!is.finite(ll)) return(1e20)
  return(-ll)
}

# --- THE UNIFIED FITTING FUNCTION (t-EV Version) ---
fit_egpd_tev_copula <- function(x, method = c("mle", "ifm"), 
                                init_par = NULL, hessian = TRUE) {
  method <- match.arg(method)
  n <- length(x)
  
  # 1. Initialization
  if (is.null(init_par)) {
    # Fit margins IID first
    init_marg <- egpd::fitegpd(x, type = 1, family = "egpd")$estimate
    
    # Estimate starting rho from Spearman/Kendall
    rho_init <- cor(x[-n], x[-1], method = "kendall")
    
    # Standard starting df for extremal-t is often around 4
    init_par <- c(init_marg, rho = rho_init, df = 4)
  }
  
  # 2. Optimization
  if (method == "mle") {
    # par_trans: log(sigma), xi, log(kappa), qlogis(rho), log(df)
    init_trans <- c(log(init_par[1]), init_par[2], log(init_par[3]), 
                    qlogis(init_par[4]), log(init_par[5]))
    
    opt <- optim(init_trans, egpd_tev_nll, x = x, method = "Nelder-Mead", 
                 control = list(maxit = 4000))
    final_trans <- opt$par
    conv <- opt$convergence
    
  } else if (method == "ifm") {
    # Stage 1: Margins only (IID)
    fit_m <- egpd::fitegpd(x, type = 1, family = "egpd")
    m_p <- fit_m$estimate
    
    # Stage 2: Copula only (fix margins)
    u <- egpd:::pegpd(x, sigma = m_p[1], xi = m_p[2], kappa = m_p[3], type = 1)
    u_pairs <- cbind(u[-n], u[-1])
    
    fit_c <- fitCopula(tevCopula(), pobs(u_pairs), method = "ml")
    
    final_trans <- c(log(m_p[1]), m_p[2], log(m_p[3]), 
                     qlogis(fit_c@estimate[1]), log(fit_c@estimate[2]))
    conv <- 0
  }
  
  # 3. Results Preparation
  est <- c(exp(final_trans[1]), 
           final_trans[2], 
           exp(final_trans[3]), 
           plogis(final_trans[4]), 
           exp(final_trans[5]))
  
  loglik <- -egpd_tev_nll(final_trans, x)
  
  # Standard Errors via Delta Method / Hessian
  se <- rep(NA, 5); vcov_mat <- NULL
  if (hessian) {
    h <- tryCatch(numDeriv::hessian(function(p) egpd_tev_nll(p, x), final_trans), 
                  error = function(e) NULL)
    if (!is.null(h)) {
      v_t <- tryCatch(solve(h), error = function(e) NULL)
      if (!is.null(v_t)) {
        # Jacobian for the back-transformations:
        # d(exp(p1))/dp1 = exp(p1)
        # d(p2)/dp2 = 1
        # d(exp(p3))/dp3 = exp(p3)
        # d(plogis(p4))/dp4 = dlogis(p4)
        # d(exp(p5))/dp5 = exp(p5)
        jac <- diag(c(est[1], 1, est[3], dlogis(final_trans[4]), est[5]))
        vcov_mat <- jac %*% v_t %*% t(jac)
        se <- sqrt(pmax(diag(vcov_mat), 0))
      }
    }
  }
  
  result <- list(estimate = est, sd = se, vcov = vcov_mat, loglik = loglik,
                 aic = -2 * loglik + 10, conv = conv, method = method)
  class(result) <- "fit_egpd_tev"
  return(result)
}