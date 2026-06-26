
source("markov_dyn_mix_copula.r")

library(copula)
# library(egpd)   # load if installed

# ============================================================
# 1. LATENT MIXTURE MODEL (same as your sampler)
# ============================================================

# Gaussian tail copula density
d_tail_cop <- function(u_star, v_star, rho) {
  cop <- normalCopula(param = rho, dim = 2)
  dCopula(cbind(u_star, v_star), copula = cop)
}

# Gumbel bulk copula density
d_bulk_cop <- function(u_star, v_star, alpha) {
  cop <- gumbelCopula(param = alpha, dim = 2)
  dCopula(cbind(u_star, v_star), copula = cop)
}

# Weight function
pi_weight <- function(u_star, v_star, theta) {
  (u_star * v_star)^theta
}

# Unnormalized latent mixture density
c_star_unnorm <- function(u_star, v_star, theta, rho, alpha) {
  w <- pi_weight(u_star, v_star, theta)
  w * d_tail_cop(u_star, v_star, rho) +
    (1 - w) * d_bulk_cop(u_star, v_star, alpha)
}

# ============================================================
# 2. BUILD THE INDUCED COPULA OBJECT FOR A GIVEN psi
#    psi = (theta, rho, alpha)
# ============================================================

build_markov_copula <- function(theta, rho, alpha,
                                grid_size = 120,
                                eps = 1e-6) {
  # ----------------------------------------------------------
  # 2.1 Normalization constant K
  # ----------------------------------------------------------
  K_val <- integrate(
    function(v) {
      sapply(v, function(vv) {
        integrate(
          function(u) c_star_unnorm(u, vv, theta, rho, alpha),
          lower = eps, upper = 1 - eps,
          rel.tol = 1e-5
        )$value
      })
    },
    lower = eps, upper = 1 - eps,
    rel.tol = 1e-5
  )$value

  # Normalized latent joint density c*(u*,v*)
  c_star <- function(u_star, v_star) {
    c_star_unnorm(u_star, v_star, theta, rho, alpha) / K_val
  }

  # ----------------------------------------------------------
  # 2.2 Latent marginal density f_{U*}(u*)
  # ----------------------------------------------------------
  f_u_star <- function(s) {
    sapply(s, function(ss) {
      integrate(
        function(v) c_star(ss, v),
        lower = eps, upper = 1 - eps,
        rel.tol = 1e-5
      )$value
    })
  }

  # ----------------------------------------------------------
  # 2.3 Build monotone spline for F_{U*} and inverse
  # ----------------------------------------------------------
  u_star_grid <- seq(eps, 1 - eps, length.out = grid_size)

  f_grid <- f_u_star(u_star_grid)

  # Numerically integrate marginal density over the grid
  # cumulative trapezoid
  du <- diff(u_star_grid)
  mid_area <- 0.5 * (f_grid[-1] + f_grid[-length(f_grid)]) * du
  F_grid <- c(0, cumsum(mid_area))

  # rescale to [0,1] to absorb numerical integration drift
  F_grid <- F_grid / max(F_grid)

  # strictly increasing safeguard
  F_grid <- pmin(pmax(F_grid, 0), 1)
  for (i in 2:length(F_grid)) {
    if (F_grid[i] <= F_grid[i - 1]) {
      F_grid[i] <- min(1, F_grid[i - 1] + 1e-10)
    }
  }

  F_u_star_spline <- splinefun(u_star_grid, F_grid, method = "monoH.FC")
  F_u_star_inv_spline <- splinefun(F_grid, u_star_grid, method = "monoH.FC")

  # Safe wrappers
  F_u_star <- function(u_star) {
    out <- F_u_star_spline(pmin(pmax(u_star, eps), 1 - eps))
    pmin(pmax(out, eps), 1 - eps)
  }

  F_u_star_inv <- function(u) {
    out <- F_u_star_inv_spline(pmin(pmax(u, eps), 1 - eps))
    pmin(pmax(out, eps), 1 - eps)
  }

  # ----------------------------------------------------------
  # 2.4 Copula density on the observed uniform scale
  #
  # c(u,v) = c*(u*,v*) / [f_{U*}(u*) f_{U*}(v*)]
  # ----------------------------------------------------------
  d_copula <- function(u, v, log = FALSE) {
    u <- pmin(pmax(u, eps), 1 - eps)
    v <- pmin(pmax(v, eps), 1 - eps)

    u_star <- F_u_star_inv(u)
    v_star <- F_u_star_inv(v)

    fu <- f_u_star(u_star)
    fv <- f_u_star(v_star)

    dens <- c_star(u_star, v_star) / (fu * fv)

    dens <- pmax(dens, .Machine$double.xmin)
    if (log) log(dens) else dens
  }

  list(
    theta = theta,
    rho = rho,
    alpha = alpha,
    K = K_val,
    c_star = c_star,
    f_u_star = f_u_star,
    F_u_star = F_u_star,
    F_u_star_inv = F_u_star_inv,
    d_copula = d_copula
  )
}

# ============================================================
# 3. NEGATIVE LOG-LIKELIHOOD FOR FULL MARKOV MODEL
# ============================================================

# Parameter vector layout:
# par = c(theta, z_rho, z_alpha, margin_par...)
#
# We optimize on unconstrained scales and transform:
#   theta > 0          via exp()
#   rho in (-1,1)      via tanh()
#   alpha > 1          via 1 + exp()
#
# You can adapt bounds if needed.

negloglik_markov_egpd <- function(par, x,
                                  margin_npar,
                                  grid_size = 120,
                                  eps = 1e-6,
                                  penalty = 1e12) {
  # ----------------------------------------------------------
  # 3.1 unpack dependence parameters
  # ----------------------------------------------------------
  theta <- exp(par[1]) # > 0
  rho <- tanh(par[2]) # in (-1,1)
  alpha <- 1 + exp(par[3]) # > 1 for Gumbel

  margin_par <- par[4:(3 + margin_npar)]

  # ----------------------------------------------------------
  # 3.2 evaluate marginal CDF / density
  # ----------------------------------------------------------
  # If EGPD parameters are invalid, the wrappers should error or
  # return non-finite values; we catch that below.
  out <- try(
    {
      u <- egpd::pegpd(x,
        kappa = margin_par[1], sigma = margin_par[2], xi = margin_par[3]
      )
      logf <- egpd::degpd_density(x,
        kappa = margin_par[1], sigma = margin_par[2],
        xi = margin_par[3], log = TRUE
      )

      # protect against exact 0/1
      u <- pmin(pmax(u, eps), 1 - eps)

      # build copula object for current dependence parameters
      cop_obj <- build_markov_copula(theta, rho, alpha,
        grid_size = grid_size,
        eps = eps
      )

      # transition log copula density
      logc <- cop_obj$d_copula(u[-length(u)], u[-1], log = TRUE)

      # full Markov log-likelihood:
      # log f(x1) + sum_{t=1}^{n-1}[ log c(u_t,u_{t+1}) + log f(x_{t+1}) ]
      ll <- sum(logf) + sum(logc)

      if (!is.finite(ll)) stop("Non-finite log-likelihood")
      ll
    },
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    return(penalty)
  }
  return(-out)
}

# ============================================================
# 4. FITTING ROUTINE
# ============================================================

fit_markov_egpd <- function(x,
                            start_dep,
                            start_margin,
                            grid_size = 120,
                            eps = 1e-6,
                            method = "L-BFGS-B",
                            hessian = TRUE,
                            control = list(maxit = 300, trace = 3)) {
  # start_dep should be on natural scale:
  #   list(theta = ..., rho = ..., alpha = ...)
  #
  # start_margin: vector of EGPD starting values on natural scale

  # unconstrained parametrization for optimizer
  par0 <- c(
    log(start_dep$theta),
    atanh(pmin(pmax(start_dep$rho, -0.999), 0.999)),
    log(start_dep$alpha - 1),
    start_margin
  )

  n_margin <- length(start_margin)

  opt <- optim(
    par = par0,
    fn = negloglik_markov_egpd,
    x = x,
    margin_npar = n_margin,
    grid_size = grid_size,
    eps = eps,
    method = method,
    hessian = hessian,
    control = control
  )

  # back-transform estimates
  est_theta <- exp(opt$par[1])
  est_rho <- tanh(opt$par[2])
  est_alpha <- 1 + exp(opt$par[3])
  est_margin <- opt$par[4:(3 + n_margin)]

  # optional covariance on transformed scale:
  vcov_unconstrained <- NULL
  se_unconstrained <- NULL
  if (!is.null(opt$hessian)) {
    invH <- try(solve(opt$hessian), silent = TRUE)
    if (!inherits(invH, "try-error")) {
      vcov_unconstrained <- invH
      se_unconstrained <- sqrt(diag(invH))
    }
  }

  list(
    convergence = opt$convergence,
    message = opt$message,
    logLik = -opt$value,
    par_unconstrained = opt$par,
    estimates = list(
      theta = est_theta,
      rho = est_rho,
      alpha = est_alpha,
      margin = est_margin
    ),
    hessian = opt$hessian,
    vcov_unconstrained = vcov_unconstrained,
    se_unconstrained = se_unconstrained,
    optim = opt
  )
}

fit <- fit_markov_egpd(
  x = x_chain[1:20],
  start_dep = list(theta = 1.2, rho = 0.3, alpha = 1.4),
  start_margin = c(1, 1, 0.4),   # example EGPD starts
  grid_size = 30,
  method = "BFGS",
  control = list(maxit = 10, trace = 6)
)

fit$estimates
fit$logLik
