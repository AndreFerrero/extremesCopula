library(copula)

# -----------------------------------------------------------
# 1. Choose dimension and true parameter
# -----------------------------------------------------------
n <- 100 # dimension of the sample vector
theta_true <- 2 # true Gumbel parameter

# -----------------------------------------------------------
# 2. Simulate ONE vector from the n-dimensional Gumbel copula
# -----------------------------------------------------------
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop)) # one vector of uniforms

# -----------------------------------------------------------
# 3. Choose margins and generate the observed sample X
# -----------------------------------------------------------
# lognormal margins just for concreteness:
mu <- 0
sigma <- 1
X <- qlnorm(U, meanlog = mu, sdlog = sigma)

# -----------------------------------------------------------
# 4. Estimate margins and obtain pseudo-observations
# -----------------------------------------------------------
mu_hat <- mean(log(X))
sigma_hat <- sd(log(X))

u_hat <- plnorm(X, meanlog = mu_hat, sdlog = sigma_hat)

# -----------------------------------------------------------
# 5. Define log-likelihood for the n-dimensional Gumbel copula
# -----------------------------------------------------------
loglik_gumbel <- function(theta, u) {
  if (theta <= 1) {
    return(-Inf)
  }
  cop <- gumbelCopula(param = theta, dim = length(u))
  ll <- log(dCopula(u, copula = cop))
  return(ll)
}

# -----------------------------------------------------------
# 6. Evaluate likelihood over a grid of theta values
# -----------------------------------------------------------
theta_grid <- seq(1.0, 5, length.out = 500)
ll_values <- sapply(theta_grid, loglik_gumbel, u = u_hat)
true_ll <- sapply(theta_grid, loglik_gumbel, u = U)
# -----------------------------------------------------------
# 7. Plot the likelihood curve
# -----------------------------------------------------------
plot(theta_grid, true_ll,
  type = "l", lwd = 2,
  xlab = expression(theta),
  ylab = "Log-likelihood"
)
lines(theta_grid, ll_values, col = "green", lwd = 2)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("True theta", "U", "Uhat"), col = c("red", "black", "green"), lwd = 2, lty = c(2, 1, 1))

joint_loglik <- function(par, x) {
  mu <- par[1]
  sigma <- par[2]
  theta <- par[3]

  if (sigma <= 0 || theta <= 1) {
    return(-Inf)
  }

  # Margins
  ll_marg <- sum(dlnorm(x, meanlog = mu, sdlog = sigma, log = TRUE))

  # Copula part
  u <- plnorm(x, meanlog = mu, sdlog = sigma)
  cop <- gumbelCopula(param = theta, dim = length(x))
  ll_cop <- log(dCopula(u, copula = cop))

  ll_marg + ll_cop
}

start <- c(mean(log(X)), sd(log(X)), 1.5)

fit <- optim(
  par = start,
  fn = joint_loglik,
  x = X,
  method = "L-BFGS-B",
  lower = c(-Inf, 1e-6, 1 + 1e-6),
  control = list(fnscale = -1)
)

fit$par


## V formulation likelihood
val_gamma <- (cos(pi / (2 * theta_true)))^theta_true
V <- stabledist::rstable(
  n = 1,
  alpha = 1 / theta_true,
  beta = 1,
  gamma = val_gamma,
  delta = 0,
  pm = 1
)

if (!is.finite(V) || V <= 0) {
  return(NULL)
}

# --- exponential variates ---
E <- rexp(n)

# --- apply inverse generator psi ---
psi = function(t, theta) exp(-t^(1/theta))
U_V <- psi(E / V, theta_true)
U_V

logp_U_Vtheta <- function(theta, U, V) {
    n <- length(U)

    n * log(V) + n * log(theta) - sum(log(U)) + (theta - 1) * sum(log(-log(U))) - V * sum((-log(U))^theta)
}

Vlik_values <- sapply(theta_grid, logp_U_Vtheta, U = U_V, V = V)
lines(theta_grid, Vlik_values, col = "orange", lwd = 2)


theta_grid[which.max(Vlik_values)]
theta_grid[which.max(true_ll)]
