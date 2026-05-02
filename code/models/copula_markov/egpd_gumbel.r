library(egpd)
library(copula)
library(numDeriv)

# --- Internal Utility: Gumbel Log-Density ---
# log_dgumbel_copula <- function(u, v, theta) {
#   u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
#   v <- pmin(pmax(v, 1e-6), 1 - 1e-6)
#   x <- -log(u); y <- -log(v)
#   t <- x^theta + y^theta
#   log_c <- -t^(1 / theta) + (theta - 1) * log(x) + (theta - 1) * log(y) +
#     (1 / theta - 2) * log(t) + log(t^(1 / theta) + theta - 1) - (log(u) + log(v))
#   return(log_c)
# }

log_dgumbel_copula <- function(u, v, theta) {
  # Numerical safety boundaries
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)

  log_density <- dCopula(cbind(u, v), gumbelCopula(theta), log = TRUE)

  return(log_density)
}


# --- Internal Utility: Joint NLL ---
egpd_gumbel_nll <- function(par_trans, x) {
  sigma <- exp(par_trans[1])
  xi <- par_trans[2]
  kappa <- exp(par_trans[3])
  theta <- exp(par_trans[4]) + 1

  if (any(1 + xi * x / sigma <= 0)) {
    return(1e20)
  }

  log_f <- egpd:::degpd_density(x, sigma = sigma, xi = xi, kappa = kappa, type = 1, log = TRUE)
  u <- pmin(pmax(egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1), 1e-10), 1 - 1e-10)

  n <- length(x)
  log_c <- log_dgumbel_copula(u[2:n], u[1:(n - 1)], theta)

  ll <- sum(log_f) + sum(log_c)
  if (!is.finite(ll)) {
    return(1e20)
  }
  -ll
}

# --- THE UNIFIED FITTING FUNCTION ---
fit_egpd_gumbel_copula <- function(x, method = c("mle", "iterative", "ifm"),
                                   init_par = NULL, hessian = TRUE, optim.method = "Nelder-Mead") {
  method <- match.arg(method)
  n <- length(x)

  # 1. Initialization
  if (is.null(init_par)) {
    init_marg <- egpd::fitegpd(x, type = 1, family = "egpd")$estimate
    tau <- cor(x[-n], x[-1], method = "kendall", use = "complete.obs")
    theta_init <- 1 / (1 - max(0.01, min(0.8, tau)))
    init_par <- c(init_marg, theta = theta_init)
  }

  # 2. Optimization Paths
  if (method == "mle") {
    init_trans <- c(log(init_par[1]), init_par[2], log(init_par[3]), log(init_par[4] - 1))
    opt <- optim(init_trans, egpd_gumbel_nll, x = x, method = optim.method, control = list(maxit = 2000, parscale = c(1, 0.1, 1, 1)))
    final_trans <- opt$par
    conv <- opt$convergence
  } else if (method == "iterative") {
    cur_t <- c(log(init_par[1]), init_par[2], log(init_par[3]), log(init_par[4] - 1))
    for (i in 1:50) {
      # Step A: Margins (Fix theta)
      opt_m <- optim(cur_t[1:3], function(p) egpd_gumbel_nll(c(p, cur_t[4]), x), method = optim.method)
      # Step B: Copula (Fix margins)
      opt_c <- optim(cur_t[4], function(p) egpd_gumbel_nll(c(opt_m$par, p), x), method = optim.method)
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
    fit_c <- fitCopula(gumbelCopula(dim = 2), pobs(u_pairs), method = "ml")
    final_trans <- c(log(m_p[1]), m_p[2], log(m_p[3]), theta = log(fit_c@estimate - 1))
    conv <- 0
  }

  # 3. Back-transform and Results
  est <- c(
    exp(final_trans[1]), final_trans[2],
    exp(final_trans[3]), exp(final_trans[4]) + 1
  )

  loglik <- - egpd_gumbel_nll(final_trans, x)

  # Standard Errors via Hessian of the joint likelihood at the solution
  se <- rep(NA, 4)
  vcov_mat <- NULL
  if (hessian) {
    h <- tryCatch(numDeriv::hessian(function(p) egpd_gumbel_nll(p, x), final_trans), error = function(e) NULL)
    if (!is.null(h)) {
      v_t <- tryCatch(solve(h), error = function(e) NULL)
      if (!is.null(v_t)) {
        jac <- diag(c(est[1], 1, est[3], est[4] - 1))
        vcov_mat <- jac %*% v_t %*% t(jac)
        se <- sqrt(pmax(diag(vcov_mat), 0))
      }
    }
  }

  structure(list(
    estimate = est, sd = se, vcov = vcov_mat, loglik = loglik,
    aic = -2 * loglik + 8, method = method, n = n
  ), class = "fit_egpd_copula")
}
