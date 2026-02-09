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
n <- 200            # sample size
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
u_par     <- plnorm(X, mu_hat, sigma_hat)  # pseudo-observations
u_ranks <- pobs(X)

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
ll_par <- sapply(theta_grid, loglik_gumbel, u = u_par)
ll_ranks <- sapply(theta_grid, loglik_gumbel, u = u_ranks)

plot(theta_grid, ll_true, type = "l", lwd = 2,
     xlab = expression(theta),
     ylab = "Log-likelihood",
     main = "Oracle vs Pseudo Likelihood")
lines(theta_grid, ll_par, col = "darkgreen", lwd = 2)
lines(theta_grid, ll_ranks, col = "orange", lwd = 2)
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
mu_grid    <- seq(-1, 1, length.out = 100)
sigma_grid <- seq(0.05, 1.5, length.out = 100)
theta_grid <- seq(1.05, 5, length.out = 100)

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

# 1. Setup the 3D Grid
# Note: Keep resolution manageable (e.g., 30-40) to avoid crashing R
N_res <- 35
mu_grid    <- seq(-0.5, 0.5, length.out = N_res)
sigma_grid <- seq(0.5, 1.5, length.out = N_res)
theta_grid <- seq(1.5, 6,   length.out = N_res)

# Create the 3D coordinate space
grid_3d <- expand.grid(mu = mu_grid, sigma = sigma_grid, theta = theta_grid)

# 2. Compute the Joint Log-Likelihood for the entire volume
# Using your existing joint_loglik function
grid_3d$LL <- apply(grid_3d, 1, function(p) {
  joint_loglik(c(p[1], p[2], p[3]), X)
})

save(grid_3d, file = "C:/Users/Andrea Ferrero/extremesCopula/sims/estim/joint/grid_ll.Rdata")
load("C:/Users/Andrea Ferrero/extremesCopula/sims/estim/joint/grid_ll.Rdata")
# 3. Visualization: 3D Isosurface
# We plot "shells" of likelihood. The inner-most shell is the "peak".
# Values far from the max are made transparent.
max_ll <- max(grid_3d$LL, na.rm = TRUE)

library(dplyr)
library(ggplot2)

# Profile Likelihood for Theta
profile_theta <- grid_3d %>%
  group_by(theta) %>%
  summarize(max_LL = max(LL, na.rm = TRUE))

ggplot(profile_theta, aes(x = theta, y = max_LL)) +
  geom_line(size = 1) +
  geom_vline(xintercept = theta_true, color = "red", linetype = "dashed") +
  labs(title = "Profile Log-Likelihood for Theta",
       subtitle = "If this is flat, Theta is not identified",
       y = "Max Log-Likelihood (optimized over mu/sigma)")

# Find the best sigma for each theta
ridge_line <- grid_3d %>%
  group_by(theta) %>%
  filter(LL == max(LL, na.rm = TRUE))

ggplot(ridge_line, aes(x = theta, y = sigma)) +
  geom_line(color = "blue", size = 1) +
  geom_point(aes(x = theta_true, y = sigma_true), color = "red", size = 3) +
  labs(title = "The Identifiability Ridge",
       x = "Copula Parameter (Theta)",
       y = "Best Marginal Scale (Sigma)",
       subtitle = "The line shows how Sigma 'compensates' for Theta")

# Filter for the top 1% of points
top_points <- grid_3d %>%
  filter(LL > quantile(LL, 0.99, na.rm = TRUE))

# Plot Sigma vs Theta for these top points
ggplot(top_points, aes(x = theta, y = sigma, color = mu)) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Top 1% Likelihood Points",
       subtitle = "Shows the 'banana' shape of the confounding")

# Pick 4 specific values of mu from your grid
mu_slices <- mu_grid[seq(1, length(mu_grid), length.out = 4)]

grid_3d %>%
  filter(mu %in% mu_slices) %>%
  ggplot(aes(x = theta, y = sigma, fill = LL)) +
  geom_tile() +
  facet_wrap(~mu, labeller = label_both) +
  scale_fill_viridis_c() +
  labs(title = "Likelihood Slices for different Mu levels",
       x = "Theta", y = "Sigma")
