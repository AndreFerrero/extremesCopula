############################################################
# FULL ASYMPTOTIC DISTORTION EVT SCRIPT
# TWO-STAGE ESTIMATION:
#
#   STEP 1:
#       Estimate copula dependence theta via DMLE
#
#   STEP 2:
#       Condition on theta_hat
#       Fit only GEV parameters inside distorted limit law
#
############################################################

library(copula)

set.seed(123)

############################################################
# SETTINGS
############################################################

n <- 200
k <- 1000

theta_true <- 2.5

family_choice <- "joe"
# family_choice <- "gumbel"

############################################################
# SIMULATE COPULA
############################################################

if (family_choice == "gumbel") {
  cop <- gumbelCopula(theta_true, dim = n)
}

if (family_choice == "joe") {
  cop <- joeCopula(theta_true, dim = n)
}

############################################################
# SIMULATE LATENT U
############################################################

U <- rCopula(k, cop)

############################################################
# MARGINS
############################################################

X <- qt(U, df = 4)

############################################################
# BLOCK MAXIMA
############################################################

M <- apply(X, 1, max)

############################################################
# PSEUDO OBSERVATIONS
############################################################

Uhat <- matrix(
  rank(as.vector(X)) / (k * n + 1),
  nrow = k,
  ncol = n
)

Yhat <- apply(Uhat, 1, max)

############################################################
# =========================================================
# DMLE FOR THETA
# =========================================================
############################################################

psi_joe <- function(t, theta) {
  1 - (1 - exp(-t))^(1 / theta)
}

psi_inv_joe <- function(u, theta) {
  -log(1 - (1 - u)^theta)
}

psi_prime_joe <- function(t, theta) {

  -(1 / theta) *
    (1 - exp(-t))^(1 / theta - 1) *
    exp(-t)
}

psi_inv_prime_joe <- function(u, theta) {

  -theta *
    (1 - u)^(theta - 1) /
    (1 - (1 - u)^theta)
}

############################################################
# DIAGONAL DENSITY
############################################################

d_delta_joe <- function(u, theta, n) {

  t <- psi_inv_joe(u, theta)

  n *
    psi_prime_joe(n * t, theta) *
    psi_inv_prime_joe(u, theta)
}

############################################################
# DMLE OBJECTIVE
############################################################

negloglik_dmle_joe <- function(theta, Y, n) {

  if (theta <= 1) return(1e10)

  dens <- d_delta_joe(Y, theta, n)

  if (any(!is.finite(dens))) return(1e10)
  if (any(dens <= 0)) return(1e10)

  -sum(log(dens))
}

############################################################
# ESTIMATE THETA
############################################################

if (family_choice == "joe") {

  fit_dmle <- optim(
    par = 2,
    fn = negloglik_dmle_joe,
    Y = Yhat,
    n = n,
    method = "L-BFGS-B",
    lower = 1.001,
    upper = 20
  )

  theta_hat <- fit_dmle$par

  cat("\n====================================\n")
  cat("JOE DMLE ESTIMATE\n")
  cat("====================================\n")

  print(theta_hat)

} else {

  ##########################################################
  # GUMBEL CLOSED FORM DMLE
  ##########################################################

  S <- sum(-log(Yhat))

  theta_hat <- log(n) /
    (log(k) - log(S))

  cat("\n====================================\n")
  cat("GUMBEL DMLE ESTIMATE\n")
  cat("====================================\n")

  print(theta_hat)
}

############################################################
# =========================================================
# GEV FUNCTIONS
# =========================================================
############################################################

gev_cdf <- function(x, mu, sigma, xi) {

  if (sigma <= 0)
    return(rep(NA, length(x)))

  z <- (x - mu) / sigma

  ##########################################################
  # GUMBEL LIMIT
  ##########################################################

  if (abs(xi) < 1e-8) {

    return(exp(-exp(-z)))

  }

  ##########################################################
  # GENERAL CASE
  ##########################################################

  support <- 1 + xi * z

  out <- rep(0, length(x))

  valid <- support > 0

  out[valid] <-
    exp(-(support[valid])^(-1 / xi))

  out
}

############################################################
# GEV LOG DENSITY
############################################################

gev_logpdf <- function(x, mu, sigma, xi) {

  if (sigma <= 0)
    return(rep(-Inf, length(x)))

  z <- (x - mu) / sigma

  ##########################################################
  # GUMBEL LIMIT
  ##########################################################

  if (abs(xi) < 1e-8) {

    t <- exp(-z)

    return(
      -log(sigma) - z - t
    )

  }

  ##########################################################
  # GENERAL CASE
  ##########################################################

  support <- 1 + xi * z

  out <- rep(-Inf, length(x))

  valid <- support > 0

  t <- support[valid]^(-1 / xi)

  out[valid] <-
    -log(sigma) -
    (1 / xi + 1) *
    log(support[valid]) -
    t

  out
}

############################################################
# GUMBEL GENERATOR
############################################################

psi_gumbel <- function(t, theta) {
  exp(-(t^(1 / theta)))
}

############################################################
# =========================================================
# ASYMPTOTIC DISTORTED LIMIT
# =========================================================
############################################################

G_asymptotic <- function(x,
                         mu,
                         sigma,
                         xi,
                         theta,
                         family) {

  H <- gev_cdf(x, mu, sigma, xi)

  H <- pmin(
    pmax(H, 1e-12),
    1 - 1e-12
  )

  V <- -log(H)

  ##########################################################
  # GUMBEL
  ##########################################################

  if (family == "gumbel") {

    rho <- 1 / theta

    return(
      psi_gumbel(
        V^(1 / rho),
        theta
      )
    )
  }

  ##########################################################
  # JOE
  ##########################################################

  if (family == "joe") {

    rho <- 1 / theta

    return(
      psi_joe(
        V^(1 / rho),
        theta
      )
    )
  }
}

############################################################
# DENSITY
############################################################

g_asymptotic <- function(x,
                         mu,
                         sigma,
                         xi,
                         theta,
                         family) {

  H <- gev_cdf(x, mu, sigma, xi)

  H <- pmin(
    pmax(H, 1e-12),
    1 - 1e-12
  )

  h <- exp(
    gev_logpdf(
      x, mu, sigma, xi
    )
  )

  V <- -log(H)

  ##########################################################
  # GUMBEL
  ##########################################################

  if (family == "gumbel") {

    return(h)

  }

  ##########################################################
  # JOE
  ##########################################################

  if (family == "joe") {

    rho <- 1 / theta

    z <- V^(1 / rho)

    psi_prime <-
      -(1 / theta) *
      (1 - exp(-z))^(1 / theta - 1) *
      exp(-z)

    dens <-
      -psi_prime *
      (1 / rho) *
      V^(1 / rho - 1) *
      h / H

    return(dens)
  }
}

############################################################
# =========================================================
# SECOND STAGE:
# FIX THETA = DMLE ESTIMATE
# FIT ONLY (mu, sigma, xi)
# =========================================================
############################################################

negloglik_stage2 <- function(par,
                             M,
                             theta_fixed,
                             family) {

  mu <- par[1]

  sigma <- exp(par[2])

  xi <- par[3]

  if (sigma <= 0)
    return(1e10)

  dens <- g_asymptotic(
    x = M,
    mu = mu,
    sigma = sigma,
    xi = xi,
    theta = theta_fixed,
    family = family
  )

  if (any(!is.finite(dens)))
    return(1e10)

  if (any(dens <= 0))
    return(1e10)

  -sum(log(dens))
}

############################################################
# INITIAL VALUES
############################################################

mu0 <- mean(M)

sigma0 <- sd(M)

xi0 <- 0.1

par0 <- c(
  mu0,
  log(sigma0),
  xi0
)

############################################################
# OPTIMIZATION
############################################################

fit_stage2 <- optim(
  par = par0,
  fn = negloglik_stage2,
  M = M,
  theta_fixed = theta_hat,
  family = family_choice,
  method = "Nelder-Mead",
  control = list(maxit = 5000)
)

############################################################
# FINAL ESTIMATES
############################################################

mu_hat <- fit_stage2$par[1]

sigma_hat <- exp(
  fit_stage2$par[2]
)

xi_hat <- fit_stage2$par[3]

############################################################
# OUTPUT
############################################################

cat("\n====================================\n")
cat("TRUE PARAMETERS\n")
cat("====================================\n")

cat("theta_true =", theta_true, "\n")

cat("\n====================================\n")
cat("TWO-STAGE ESTIMATES\n")
cat("====================================\n")

cat("theta_hat =", theta_hat, "\n")

cat("mu_hat    =", mu_hat, "\n")

cat("sigma_hat =", sigma_hat, "\n")

cat("xi_hat    =", xi_hat, "\n")

############################################################
# =========================================================
# QUANTILE FUNCTION
# =========================================================
############################################################

gev_quantile <- function(p,
                         mu,
                         sigma,
                         xi) {

  if (sigma <= 0)
    return(rep(NA, length(p)))

  ##########################################################
  # GUMBEL LIMIT
  ##########################################################

  if (abs(xi) < 1e-8) {

    return(
      mu - sigma * log(-log(p))
    )

  }

  ##########################################################
  # GENERAL CASE
  ##########################################################

  mu +
    (sigma / xi) *
    ((-log(p))^(-xi) - 1)
}

############################################################
# DISTORTED QUANTILE
############################################################

model_quantile <- function(p,
                           mu,
                           sigma,
                           xi,
                           theta,
                           family) {

  eps <- 1e-10

  p <- pmin(
    pmax(p, eps),
    1 - eps
  )

  ##########################################################
  # GUMBEL
  ##########################################################

  if (family == "gumbel") {

    return(
      gev_quantile(
        p,
        mu,
        sigma,
        xi
      )
    )
  }

  ##########################################################
  # JOE
  ##########################################################

  if (family == "joe") {

    rho <- 1 / theta

    t <- psi_inv_joe(
      p,
      theta
    )

    return(
      gev_quantile(
        exp(-(t^rho)),
        mu,
        sigma,
        xi
      )
    )
  }
}

############################################################
# =========================================================
# DIAGNOSTICS
# =========================================================
############################################################

n_obs <- length(M)

prob_grid <-
  (1:n_obs) / (n_obs + 1)

############################################################
# MODEL QUANTILES
############################################################

q_theory <- model_quantile(
  prob_grid,
  mu_hat,
  sigma_hat,
  xi_hat,
  theta_hat,
  family_choice
)

q_emp <- sort(M)

############################################################
# PLOTS
############################################################

par(mfrow = c(1,3))

############################################################
# QQ PLOT
############################################################

plot(
  q_theory,
  q_emp,
  pch = 16,
  cex = 0.5,
  main = "QQ Plot"
)

abline(
  0, 1,
  col = "red"
)

############################################################
# BOOTSTRAP QQ CONFIDENCE BANDS
############################################################

B <- 200  # bootstrap replications
n_obs <- length(M)
prob_grid <- (1:n_obs) / (n_obs + 1)

q_boot <- matrix(NA, nrow = B, ncol = n_obs)

for (b in 1:B) {

  # simulate maxima from fitted model
  if (family_choice == "joe") {
    cop_b <- joeCopula(theta_hat, dim = n)
  } else {
    cop_b <- gumbelCopula(theta_hat, dim = n)
  }

  U_b <- rCopula(k, cop_b)
  X_b <- qt(U_b, df = 4)
  M_b <- apply(X_b, 1, max)

  # recompute fitted quantiles (same parameters, so only randomness)
  q_boot[b, ] <- sort(M_b)
}

# pointwise bands
q_low  <- apply(q_boot, 2, quantile, 0.025)
q_high <- apply(q_boot, 2, quantile, 0.975)
q_med  <- apply(q_boot, 2, quantile, 0.5)

plot(q_theory, q_emp,
     pch = 16, cex = 0.5,
     main = "QQ plot with 95% bands",
     xlab = "Theoretical quantiles",
     ylab = "Empirical maxima")

# band region
lines(q_theory, q_low, col = "gray", lwd = 2)
lines(q_theory, q_high, col = "gray", lwd = 2)

# median bootstrap line
lines(q_theory, q_med, col = "blue", lwd = 2)

# diagonal
abline(0,1,col="red")
