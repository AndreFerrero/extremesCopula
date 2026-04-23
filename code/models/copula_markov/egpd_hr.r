library(egpd)
library(copula)
library(numDeriv)

# --- Internal Utility: Hüsler-Reiss Log-Density ---
# We use the Pickands dependence function approach for numerical stability
log_dhr_copula <- function(u, v, lambda) {
  # Numerical safety boundaries
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  
  # The density of an HR copula is given by:
  # c(u,v) = C(u,v)/(u*v) * [ (Phi(z) - phi(z)/(2*lambda*x)) * (Phi(w) - phi(w)/(2*lambda*y)) + ... ]
  # A more stable way is using the dCopula function from the copula package:
  log_density <- dCopula(cbind(u, v), huslerReissCopula(lambda), log = TRUE)
  
  return(log_density)
}

# --- Internal Utility: Joint NLL for HR-EGPD ---
egpd_hr_nll <- function(par_trans, x) {
  # Back-transform parameters
  sigma <- exp(par_trans[1])
  xi    <- par_trans[2] # xi is real-valued
  kappa <- exp(par_trans[3])
  lambda <- exp(par_trans[4]) # lambda > 0
  
  # Domain check for EGPD
  if (any(1 + xi * x / sigma <= 0)) return(1e20)
  
  # Marginal Log-Likelihood (EGPD)
  log_f <- egpd:::degpd_density(x, sigma = sigma, xi = xi, kappa = kappa, type = 1, log = TRUE)
  
  # Probability Integral Transform (PIT)
  u <- pmin(pmax(egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1), 1e-10), 1 - 1e-10)
  
  # Dependence Log-Likelihood (Hüsler-Reiss)
  n <- length(x)
  # Model x_t and x_{t+1} as a Markov Chain
  log_c <- log_dhr_copula(u[1:(n - 1)], u[2:n], lambda)
  
  ll <- sum(log_f) + sum(log_c)
  
  if (!is.finite(ll)) return(1e20)
  return(-ll)
}

# --- THE UNIFIED FITTING FUNCTION (Hüsler-Reiss Version) ---
fit_egpd_hr_copula <- function(x, method = c("mle", "ifm"), 
                               init_par = NULL, hessian = TRUE) {
  method <- match.arg(method)
  n <- length(x)
  
  # 1. Initialization
  if (is.null(init_par)) {
    # Fit margins IID first for a good starting point
    init_marg <- egpd::fitegpd(x, type = 1, family = "egpd")$estimate
    
    # Estimate lambda using Kendall's Tau relationship for HR
    # Note: tau = 1 - integral... no closed form, so we use a common heuristic
    tau <- cor(x[-n], x[-1], method = "kendall", use = "complete.obs")
    # Heuristic for starting lambda: Higher tau -> Smaller lambda
    lambda_init <- sqrt(2) * qnorm((tau + 1) / 2) # Inverse of approx tau logic
    lambda_init <- pmax(0.1, pmin(5, 1/lambda_init)) 
    
    init_par <- c(init_marg, lambda = lambda_init)
  }
  
  # 2. Optimization
  if (method == "mle") {
    # par_trans: log(sigma), xi, log(kappa), log(lambda)
    init_trans <- c(log(init_par[1]), init_par[2], log(init_par[3]), log(init_par[4]))
    
    opt <- optim(init_trans, egpd_hr_nll, x = x, method = "Nelder-Mead", 
                 control = list(maxit = 3000))
    final_trans <- opt$par
    conv <- opt$convergence
    
  } else if (method == "ifm") {
    # Stage 1: Margins only (IID)
    fit_m <- egpd::fitegpd(x, type = 1, family = "egpd")
    m_p <- fit_m$estimate
    
    # Stage 2: Copula only (fix margins)
    u <- egpd:::pegpd(x, sigma = m_p[1], xi = m_p[2], kappa = m_p[3], type = 1)
    u_pairs <- cbind(u[-n], u[-1])
    # Use copula package for stage 2
    fit_c <- fitCopula(huslerReissCopula(), pobs(u_pairs), method = "ml")
    
    final_trans <- c(log(m_p[1]), m_p[2], log(m_p[3]), lambda = log(fit_c@estimate))
    conv <- 0
  }
  
  # 3. Results Preparation
  est <- c(exp(final_trans[1]), 
           final_trans[2], 
           exp(final_trans[3]), 
           exp(final_trans[4]))
  
  loglik <- -egpd_hr_nll(final_trans, x)
  
  # Standard Errors via Delta Method / Hessian
  se <- rep(NA, 4); vcov_mat <- NULL
  if (hessian) {
    h <- tryCatch(numDeriv::hessian(function(p) egpd_hr_nll(p, x), final_trans), 
                  error = function(e) NULL)
    if (!is.null(h)) {
      v_t <- tryCatch(solve(h), error = function(e) NULL)
      if (!is.null(v_t)) {
        # Jacobian for log-transformed parameters
        jac <- diag(c(est[1], 1, est[3], est[4]))
        vcov_mat <- jac %*% v_t %*% t(jac)
        se <- sqrt(pmax(diag(vcov_mat), 0))
      }
    }
  }
  
  result <- list(estimate = est, sd = se, vcov = vcov_mat, loglik = loglik,
                 aic = -2 * loglik + 8, conv = conv, method = method)
  class(result) <- "fit_egpd_hr"
  return(result)
}
