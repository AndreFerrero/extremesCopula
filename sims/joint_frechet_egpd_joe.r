library(copula)
library(evd)
library(tidyverse)

# Note: Using explicit formulas for EGPD to ensure the Joint Likelihood is transparent
# EGPD CDF: F(x) = [1 - (1 + xi*x/sigma)^(-1/xi)]^kappa
pegpd_custom <- function(x, kappa, sigma, xi) {
  # Use ifelse to ensure vectorization/matrix structure is preserved
  res <- if (abs(xi) < 1e-10) {
    (1 - exp(-x/sigma))^kappa
  } else {
    # Ensure we don't take powers of negative numbers
    inner <- 1 + xi * x / sigma
    (1 - pmax(0, inner)^(-1/xi))^kappa
  }
  # If x was a matrix, res is currently a vector. We must return it to matrix form.
  if (is.matrix(x)) return(matrix(res, nrow = nrow(x), ncol = ncol(x)))
  return(res)
}

degpd_custom <- function(x, kappa, sigma, xi) {
  res <- if (abs(xi) < 1e-10) {
    (kappa/sigma) * exp(-x/sigma) * (1 - exp(-x/sigma))^(kappa - 1)
  } else {
    t <- pmax(1e-10, 1 + xi * x / sigma)
    (kappa/sigma) * t^(-1/xi - 1) * (1 - t^(-1/xi))^(kappa - 1)
  }
  if (is.matrix(x)) return(matrix(res, nrow = nrow(x), ncol = ncol(x)))
  return(res)
}

# ==============================================================================
# 1. SETUP & SIMULATION
# ==============================================================================
set.seed(101)
p <- 1000       # Number of blocks
k <- 50        # Block length 
alpha_true <- 2
theta_true <- 2

# We generate data using Joe Copula + Frechet Margins
generate_blocks_frechet <- function(theta, p, k, alpha) {
  cop <- joeCopula(theta, dim = k)
  u <- rCopula(p, cop)
  x <- (-log(u))^(-1/alpha) 
  return(x)
}

data_matrix <- generate_blocks(theta_true, p, k, alpha_true)

# ==============================================================================
# 2. JOINT LIKELIHOOD FUNCTIONS
# ==============================================================================

# LIKELIHOOD A: Correct Frechet Margin + Joe Copula
lik_frechet_joe <- function(params, data) {
  alpha <- params[1]; theta <- params[2]
  if(alpha <= 0 || theta <= 1) return(1e10)
  
  u <- exp(-(data)^(-alpha))
  u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
  
  log_f <- log(alpha) - (alpha + 1) * log(data) - (data)^(-alpha)
  log_c <- dCopula(u, joeCopula(theta, dim = ncol(data)), log = TRUE)
  
  return(-(sum(log_f) + sum(log_c)))
}

# LIKELIHOOD B: Flexible EGPD Margin + Joe Copula
lik_egpd_joe <- function(params, data) {
  # params: kappa, sigma, xi, theta
  ka <- params[1]; si <- params[2]; xi <- params[3]; th <- params[4]
  if(ka <= 0 || si <= 0 || th <= 1) return(1e10)
  
  u <- pegpd_custom(data, ka, si, xi)
  u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
  
  log_f <- log(degpd_custom(data, ka, si, xi))
  log_c <- dCopula(u, joeCopula(th, dim = ncol(data)), log = TRUE)
  
  return(-(sum(log_f) + sum(log_c)))
}

# ==============================================================================
# 3. OPTIMIZATION
# ==============================================================================

cat("Fitting Joint Frechet-Joe Model...\n")
fit_frechet <- optim(par = c(2, 2), fn = lik_frechet_joe, data = data_matrix)

cat("Fitting Joint EGPD-Joe Model (Flexible)...\n")
# Starting values: kappa=1 (standard GPD), sigma=1, xi=0.5, theta=2
fit_egpd <- optim(par = c(1, 1, 0.5, 2), fn = lik_egpd_joe, data = data_matrix)

# ==============================================================================
# 4. RECONSTRUCTION & PLOTTING
# ==============================================================================
eval_x <- seq(0.1, 40, length.out = 300)
diag_joe <- function(u, theta, k) 1 - (1 - (1 - (1-u)^theta)^k)^(1/theta)

# 1. Reconstruct from Frechet Fit
F_frechet <- exp(-(eval_x)^(-fit_frechet$par[1]))
G_frechet_recon <- diag_joe(F_frechet, fit_frechet$par[2], k)

# 2. Reconstruct from EGPD Fit
F_egpd <- pegpd_custom(eval_x, fit_egpd$par[1], fit_egpd$par[2], fit_egpd$par[3])
G_egpd_recon <- diag_joe(F_egpd, fit_egpd$par[4], k)

# Empirical Maxima
emp_max <- apply(data_matrix, 1, max)

# Results Comparison Table
results <- data.frame(
  Model = c("True", "Joint Frechet", "Joint EGPD"),
  Theta = c(theta_true, fit_frechet$par[2], fit_egpd$par[4]),
  Effective_Xi = c(1/alpha_true, 1/fit_frechet$par[1], fit_egpd$par[3])
)
print(results)

# Plot
ggplot() +
  stat_ecdf(aes(emp_max, color = "Empirical Block Maxima"), linetype = "dashed", size = 1) +
  geom_line(aes(eval_x, G_frechet_recon, color = "Joint Frechet-Joe Recon"), size = 1) +
  geom_line(aes(eval_x, G_egpd_recon, color = "Joint EGPD-Joe Recon"), size = 1) +
  labs(title = "Joint Likelihood: Parametric vs. Flexible Margins",
       subtitle = "Comparing block maxima reconstruction under complex Joe dependence",
       x = "Max Value", y = "Probability (CDF)") +
  theme_minimal() +
  scale_color_manual(values = c("Empirical Block Maxima" = "black", 
                                "Joint Frechet-Joe Recon" = "blue", 
                                "Joint EGPD-Joe Recon" = "red"))
