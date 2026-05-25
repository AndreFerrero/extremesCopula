# ==============================================================================
# COMPARISON: DISTORTED MODEL VS. STANDARD PACKAGE (extRemes)
# ==============================================================================
# install.packages("extRemes")
library(extRemes)
source("R/utils.R")
load_all()

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
n <- 1000
k <- 10000
theta_true <- 3
family_choice <- "joe" # Joe is non-GEV, so the package fit should be worse
seed <- 46

qfrechet <- function(u, xi = 0.25) {
  # alpha = 1/xi
  return( (-log(u))^(-xi) )
}
qpareto <- function(u, xi = 0.25) {
  # Standard Pareto starting at 1
  # u = 1 - x^(-1/xi)  =>  x = (1-u)^(-xi)
  return( (1 - u)^(-xi) )
}

# How to use it in your script:
sim <- simulate_data(k, n, theta_true, family_choice, 
                     margin_fn = function(u) qpareto(u), 
                     seed = seed)

pseudo <- compute_pseudo_obs(sim$X)

# ------------------------------------------------------------------------------
# 1. Standard Block Maxima using 'extRemes' package
# ------------------------------------------------------------------------------
# fevd = "Fitting Extreme Value Distributions"
# sim$M is the vector of k block maxima
fit_extRemes <- fevd(sim$M, type = "GEV", method = "MLE")

# Extract parameters from package object
# Note: extRemes uses (location, scale, shape). Shape = xi.
params_pkg <- distill(fit_extRemes) 
mu_pkg    <- params_pkg["location"]
sigma_pkg <- params_pkg["scale"]
xi_pkg    <- params_pkg["shape"]

# ------------------------------------------------------------------------------
# 2. Your Two-Stage Approach (Stage 1: DMLE, Stage 2: Distorted GEV)
# ------------------------------------------------------------------------------
# Step 1: Estimate dependency component
theta_dmle <- estimate_theta_dmle(pseudo$Yhat, n, k, family_choice)

# Step 2: Estimate GEV components given theta
fit_distorted <- fit_joint(sim$M, family_choice, theta_init = theta_dmle)

# ------------------------------------------------------------------------------
# Comparison Results
# ------------------------------------------------------------------------------
cat("\n================================================")
cat("\n           COMPARISON OF ESTIMATORS             ")
cat("\n================================================")
cat("\nTrue Dependency (Theta):", theta_true)
cat("\nEstimated Dependency:   ", theta_dmle)

cat("\n\n--- Standard GEV (extRemes Package) ---")
cat("\nLocation (mu):    ", mu_pkg)
cat("\nScale (sigma):    ", sigma_pkg)
cat("\nShape (xi):       ", xi_pkg)
summary(fit_extRemes)

cat("\n\n--- Distorted Model ---")
cat("\nLocation (mu):    ", fit_distorted$mu)
cat("\nScale (sigma):    ", fit_distorted$sigma)
cat("\nShape (xi):       ", fit_distorted$xi)
# Calculate a rough AIC for comparison
lik_val <- -fit_distorted$value
aic_val <- 2*4 - 2*lik_val # 4 parameters: mu, sigma, xi, theta
bic_val <- log(k)*4 - 2*lik_val
cat("\nApprox AIC:       ", aic_val)
cat("\nApprox BIC:       ", bic_val)
cat("\n================================================\n")

# ------------------------------------------------------------------------------
# Visual Comparison
# ------------------------------------------------------------------------------
par(mfrow=c(1,2))

# 1. Package Quantiles
# We create a theoretical grid based on standard GEV
prob_grid <- (1:k)/(k+1)
q_theory_pkg <- qevd(prob_grid, loc=mu_pkg, scale=sigma_pkg, shape=xi_pkg)
plot(q_theory_pkg, sort(sim$M), main="extRemes (Standard GEV)", 
     xlab="Theoretical", ylab="Empirical", pch=16, cex=0.5)
abline(0,1, col="red")

# 2. Your Model Quantiles
plot_qq_bootstrap(sim$M, fit_distorted, n, k, B = 200)
# (Added sub-title to the plot function output)


# B = 200
# boot_quantiles <- numeric(B)
# for(b in 1:B) {
#   # 1. Resample the blocks (rows of your data matrix X)
#   idx <- sample(1:k, k, replace = TRUE)
#   X_boot <- sim$X[idx, ]
#   M_boot <- sim$M[idx] # Maxima of those rows
  
#   # 2. Stage 1: DMLE
#   pseudo_boot <- compute_pseudo_obs(X_boot)
#   theta_b <- estimate_theta_dmle(pseudo_boot$Yhat, n, k, family_choice)
  
#   # 3. Stage 2: MLE
#   fit_b <- fit_twostage_gev(M_boot, theta_b, family_choice)
  
#   # 4. Calculate Quantile of interest (e.g. 0.99)
#   boot_quantiles[b] <- model_quantile(0.99, fit_b$mu, fit_b$sigma, fit_b$xi, fit_b$theta, family_choice)
# }

# # CI for the 0.99 quantile
# quantile(boot_quantiles, probs = c(0.025, 0.975))