############################################################
# FULL ASYMPTOTIC DISTORTION EVT SCRIPT (WITH QQ/PP)
############################################################

library(copula)

set.seed(123)

############################################################
# SETTINGS
############################################################

n <- 1000
k <- 5000

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
# 1. DMLE (DIAGONAL MLE)
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
  -theta * (1 - u)^(theta - 1) /
    (1 - (1 - u)^theta)
}

d_delta_joe <- function(u, theta, n) {
  t <- psi_inv_joe(u, theta)

  n * psi_prime_joe(n * t, theta) * psi_inv_prime_joe(u, theta)
}

negloglik_dmle_joe <- function(theta, Y, n) {

  if (theta <= 1) return(1e10)

  dens <- d_delta_joe(Y, theta, n)

  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e10)

  -sum(log(dens))
}

############################################################
# DMLE ESTIMATION
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

  theta_dmle_joe <- fit_dmle$par

  cat("\nJoe DMLE:", theta_dmle_joe, "\n")

} else {

  S <- sum(-log(Yhat))

  theta_dmle_gumbel <- log(n) / (log(k) - log(S))
  cat("\n Gumbel DMLE:", theta_dmle_gumbel, "\n")
}

############################################################
# =========================================================
# GEV FUNCTIONS
# =========================================================
############################################################

gev_cdf <- function(x, mu, sigma, xi) {

  if (sigma <= 0) return(rep(NA, length(x)))

  z <- (x - mu) / sigma

  if (abs(xi) < 1e-8) {
    return(exp(-exp(-z)))
  }

  support <- 1 + xi * z
  out <- rep(0, length(x))
  valid <- support > 0

  out[valid] <- exp(-(support[valid])^(-1 / xi))
  out
}

gev_logpdf <- function(x, mu, sigma, xi) {

  if (sigma <= 0) return(rep(-Inf, length(x)))

  z <- (x - mu) / sigma

  if (abs(xi) < 1e-8) {
    t <- exp(-z)
    return(-log(sigma) - z - t)
  }

  support <- 1 + xi * z

  out <- rep(-Inf, length(x))
  valid <- support > 0

  t <- support[valid]^(-1 / xi)

  out[valid] <-
    -log(sigma) -
    (1 / xi + 1) * log(support[valid]) -
    t

  out
}

############################################################
# GENERATORS
############################################################

psi_gumbel <- function(t, theta) {
  exp(-(t^(1 / theta)))
}

psi_joe <- function(t, theta) {
  1 - (1 - exp(-t))^(1 / theta)
}

psi_inv_joe <- function(p, theta) {
  -log(1 - (1 - p)^theta)
}

############################################################
# ASYMPTOTIC MODEL
############################################################

G_asymptotic <- function(x, mu, sigma, xi, theta, family) {

  H <- gev_cdf(x, mu, sigma, xi)
  H <- pmin(pmax(H, 1e-12), 1 - 1e-12)

  V <- -log(H)

  if (family == "gumbel") {

    rho <- 1 / theta
    return(psi_gumbel(V^(1 / rho), theta))

  }

  if (family == "joe") {

    rho <- 1 / theta
    return(psi_joe(V^(1 / rho), theta))

  }
}

############################################################
# DENSITY
############################################################

g_asymptotic <- function(x, mu, sigma, xi, theta, family) {

  H <- gev_cdf(x, mu, sigma, xi)
  H <- pmin(pmax(H, 1e-12), 1 - 1e-12)

  h <- exp(gev_logpdf(x, mu, sigma, xi))
  V <- -log(H)

  if (family == "gumbel") return(h)

  if (family == "joe") {

    rho <- 1 / theta
    z <- V^(1 / rho)

    psi_prime <- -(1 / theta) *
      (1 - exp(-z))^(1 / theta - 1) *
      exp(-z)

    dens <- -psi_prime *
      (1 / rho) *
      V^(1 / rho - 1) *
      h / H

    return(dens)
  }
}

############################################################
# NEG LOG-LIKELIHOOD
############################################################

negloglik <- function(par, M, family) {

  mu <- par[1]
  sigma <- exp(par[2])
  xi <- par[3]
  theta <- par[4]

  if (sigma <= 0 || theta <= 1) return(1e10)

  dens <- g_asymptotic(M, mu, sigma, xi, theta, family)

  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e10)

  -sum(log(dens))
}

############################################################
# INITIAL VALUES
############################################################

mu0 <- mean(M)
sigma0 <- sd(M)
xi0 <- 0.1

theta0 <- if (family_choice == "joe") {
  max(theta_dmle_joe, 1.05)
} else {
  max(theta_dmle_gumbel, 1.05)
}

par0 <- c(mu0, log(sigma0), xi0, theta0)

############################################################
# OPTIMIZATION
############################################################

fit <- optim(
  par = par0,
  fn = negloglik,
  M = M,
  family = family_choice,
  method = "Nelder-Mead",
  control = list(maxit = 5000)
)

############################################################
# ESTIMATES
############################################################

mu_hat <- fit$par[1]
sigma_hat <- exp(fit$par[2])
xi_hat <- fit$par[3]
theta_hat <- fit$par[4]

############################################################
# OUTPUT
############################################################

cat("\nTRUE theta:\n", theta_true,
    "\nEST theta:\n", theta_hat, "\n")

cat("GEV estimates:\n mu:", mu_hat, "\n",
"sigma:", sigma_hat, "\n", "xi:", xi_hat, "\n")

############################################################
# =========================================================
# CLOSED-FORM QUANTILE
# =========================================================
############################################################

gev_quantile <- function(p, mu, sigma, xi) {

  if (sigma <= 0) return(rep(NA, length(p)))

  if (abs(xi) < 1e-8) {
    return(mu - sigma * log(-log(p)))
  }

  mu + (sigma / xi) * ((-log(p))^(-xi) - 1)
}

model_quantile <- function(p, mu, sigma, xi, theta, family) {

  eps <- 1e-10
  p <- pmin(pmax(p, eps), 1 - eps)

  if (family == "gumbel") {
    return(gev_quantile(p, mu, sigma, xi))
  }

  if (family == "joe") {
    rho <- 1 / theta
    t <- psi_inv_joe(p, theta)
    return(gev_quantile(exp(-(t^rho)), mu, sigma, xi))
  }
}

############################################################
# DIAGNOSTICS
############################################################

n_obs <- length(M)
prob_grid <- (1:n_obs) / (n_obs + 1)

q_theory <- model_quantile(prob_grid,
                           mu_hat, sigma_hat, xi_hat,
                           theta_hat, family_choice)

q_emp <- sort(M)

par(mfrow = c(1,3))

plot(q_theory, q_emp,
     pch = 16, cex = 0.5,
     main = "QQ plot")
abline(0,1,col="red")

plot(sort(G_asymptotic(M, mu_hat, sigma_hat, xi_hat,
                       theta_hat, family_choice)),
     prob_grid,
     pch=16, cex=0.5,
     main="PP plot")
abline(0,1,col="red")

tail_idx <- floor(0.9*n_obs):n_obs

plot(q_theory[tail_idx], q_emp[tail_idx],
     pch=16, cex=0.7,
     main="TAIL QQ")
abline(0,1,col="blue")
par(mfrow = c(1,1))
