# ============================================================
# Gumbel Copula Joint Likelihood Analysis
# ============================================================
# Simulate data, compute pseudo-observations, 
# and explore likelihood surfaces with ridge extraction
# ============================================================

library(copula)

set.seed(1)

# ------------------------------------------------------------
# 1. Dimension and true copula parameter
# ------------------------------------------------------------
n <- 100            # sample size
theta_true <- 4     # true Gumbel copula parameter

# ------------------------------------------------------------
# 2. Draw from Gumbel copula
# ------------------------------------------------------------
cop <- gumbelCopula(theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))  # single draw

# ------------------------------------------------------------
# 3. Margins and observed data
# ------------------------------------------------------------
mu_true <- 0
sigma_true <- 1
X <- qlnorm(U, mu_true, sigma_true)  # lognormal data

# ------------------------------------------------------------
# 4. Plug-in marginal estimation (pseudo-observations)
# ------------------------------------------------------------
mu_hat    <- mean(log(X))
sigma_hat <- sd(log(X))
u_hat     <- plnorm(X, mu_hat, sigma_hat)  # pseudo-observations

# ------------------------------------------------------------
# 5. Copula log-likelihood functions
# ------------------------------------------------------------
loglik_gumbel <- function(theta, u) {
  if (theta <= 1) return(-Inf)
  cop <- gumbelCopula(theta, dim = length(u))
  log(dCopula(u, cop))
}

# ------------------------------------------------------------
# 6. Likelihood curves: Oracle vs pseudo
# ------------------------------------------------------------
theta_grid <- seq(1.05, 5, length.out = 100)

ll_true <- sapply(theta_grid, loglik_gumbel, u = U)
ll_pseudo <- sapply(theta_grid, loglik_gumbel, u = u_hat)

plot(theta_grid, ll_true, type = "l", lwd = 2,
     xlab = expression(theta),
     ylab = "Log-likelihood",
     main = "Oracle vs Pseudo Likelihood")
lines(theta_grid, ll_pseudo, col = "darkgreen", lwd = 2)
abline(v = theta_true, col = "red", lty = 2)
legend("topright",
       legend = c("Oracle (true U)", "Pseudo (Uhat)", "True theta"),
       col = c("black", "darkgreen", "red"),
       lwd = 2, lty = c(1,1,2))

# ------------------------------------------------------------
# 7. Joint likelihood: (mu, sigma, theta)
# ------------------------------------------------------------
joint_loglik <- function(par, x) {
  mu    <- par[1]
  sigma <- par[2]
  theta <- par[3]
  
  if (sigma <= 0 || theta <= 1) return(-Inf)
  
  # Marginal likelihood (lognormal)
  ll_marg <- sum(dlnorm(x, mu, sigma, log = TRUE))
  
  # Copula likelihood (Gumbel)
  u <- plnorm(x, mu, sigma)
  cop <- gumbelCopula(theta, dim = length(x))
  ll_cop <- log(dCopula(u, cop))
  
  ll_marg + ll_cop
}

# ------------------------------------------------------------
# 8. Parameter grids for slicing
# ------------------------------------------------------------
mu_grid    <- seq(-1, 1, length.out = 40)
sigma_grid <- seq(0.05, 1.5, length.out = 40)
theta_grid <- seq(1.05, 5, length.out = 40)

# Fix parameters for 2D slices
mu_fix    <- mu_true
sigma_fix <- sigma_true
theta_fix <- theta_true

# ------------------------------------------------------------
# 9. Likelihood surfaces
# ------------------------------------------------------------
# 9a. theta vs mu | sigma fixed
LL_theta_mu <- outer(theta_grid, mu_grid, Vectorize(function(th, mu) {
  joint_loglik(c(mu, sigma_fix, th), X)
}))

persp(theta_grid, mu_grid, LL_theta_mu,
      xlab = expression(theta),
      ylab = expression(mu),
      zlab = "Log-likelihood",
      main = "Joint log-likelihood surface: (theta, mu) | sigma fixed",
      theta = 40, phi = 25, col = "lightblue", ticktype = "detailed")

# 9b. theta vs sigma | mu fixed
LL_theta_sigma <- outer(theta_grid, sigma_grid, Vectorize(function(th, sg) {
  joint_loglik(c(mu_fix, sg, th), X)
}))

persp(theta_grid, sigma_grid, LL_theta_sigma,
      xlab = expression(theta),
      ylab = expression(sigma),
      zlab = "Log-likelihood",
      main = "Joint log-likelihood surface: (theta, sigma) | mu fixed",
      theta = 40, phi = 25, col = "lightgreen", ticktype = "detailed")

# 9c. mu vs sigma | theta fixed
LL_mu_sigma <- outer(mu_grid, sigma_grid, Vectorize(function(mu, sg) {
  joint_loglik(c(mu, sg, theta_fix), X)
}))

persp(mu_grid, sigma_grid, LL_mu_sigma,
      xlab = expression(mu),
      ylab = expression(sigma),
      zlab = "Log-likelihood",
      main = "Joint log-likelihood surface: (mu, sigma) | theta fixed",
      theta = 40, phi = 25, col = "lightpink", ticktype = "detailed")

# ------------------------------------------------------------
# 10. Contours for 2D slices
# ------------------------------------------------------------
contour(theta_grid, mu_grid, LL_theta_mu,
        xlab = expression(theta), ylab = expression(mu),
        main = "Contour: (theta, mu) | sigma fixed")
contour(theta_grid, sigma_grid, LL_theta_sigma,
        xlab = expression(theta), ylab = expression(sigma),
        main = "Contour: (theta, sigma) | mu fixed")
contour(mu_grid, sigma_grid, LL_mu_sigma,
        xlab = expression(mu), ylab = expression(sigma),
        main = "Contour: (mu, sigma) | theta fixed")

# ------------------------------------------------------------
# 11. Heatmaps for 2D slices
# ------------------------------------------------------------
image(theta_grid, mu_grid, LL_theta_mu,
      xlab = expression(theta), ylab = expression(mu),
      main = "Heatmap: (theta, mu) | sigma fixed",
      col = terrain.colors(50))
image(theta_grid, sigma_grid, LL_theta_sigma,
      xlab = expression(theta), ylab = expression(sigma),
      main = "Heatmap: (theta, sigma) | mu fixed",
      col = terrain.colors(50))
image(mu_grid, sigma_grid, LL_mu_sigma,
      xlab = expression(mu), ylab = expression(sigma),
      main = "Heatmap: (mu, sigma) | theta fixed",
      col = terrain.colors(50))

# ------------------------------------------------------------
# 12. Ridge extraction: theta vs sigma | mu fixed
# ------------------------------------------------------------
LL_clean <- LL_theta_sigma
LL_clean[!is.finite(LL_clean)] <- NA

# Find sigma that maximizes likelihood for each theta
ridge_sigma_idx <- sapply(seq_len(ncol(LL_clean)), function(j) {
  col <- LL_clean[, j]
  if (all(is.na(col))) return(NA_integer_)
  which.max(col)
})

ridge_sigma <- sigma_grid[ridge_sigma_idx]
ridge_theta <- theta_grid
ridge_LL    <- apply(LL_clean, 2, max, na.rm = TRUE)

# Overlay ridge on contour
contour(theta_grid, sigma_grid, LL_clean,
        xlab = expression(theta),
        ylab = expression(sigma),
        main = "Ridge in (theta, sigma) likelihood")
lines(ridge_theta, ridge_sigma, col = "red", lwd = 2)
points(ridge_theta, ridge_sigma, col = "red", pch = 19, cex = 0.6)

# Fisher information
library(numDeriv)

psi_true <- c(mu_true, sigma_true, theta_true)

negloglik <- function(par) {
  -joint_loglik(par, X)
}

H <- hessian(negloglik, psi_true)
I_obs <- H
eigen(I_obs)$values
eigen(I_obs)$vectors

true_H <- hessian(function(th) loglik_gumbel(th, U), theta)
true_H
eigen(true_H)
