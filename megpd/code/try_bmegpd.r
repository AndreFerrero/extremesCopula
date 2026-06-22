library(copula)
library(egpd)
library(splines)
library(rstan)
library(tidyverse)
library(bayesplot)

# --- 1. DATA SIMULATION ---
set.seed(46)
n <- 1000
gumbelC <- gumbelCopula(param = 5, dim = 2)
u <- rCopula(n, gumbelC)

# Generate observations with different tail weights
x <- egpd::qegpd(u[,1], kappa = 2, sigma = 1, xi = 0.2)
y <- egpd::qegpd(u[,2], kappa = 2, sigma = 1, xi = 0.2)

# Calculate Radius for basis expansion
R_obs <- x + y

# --- 2. BASIS EXPANSION (The "Fixed" Part) ---
# Create a B-spline basis matrix for the observed radii
K_basis <- 20
# We use a log-transform on R for the basis to give more knots to the bulk
B_matrix <- bs(R_obs, df = K_basis, degree = 3, intercept = TRUE)

# --- 3. RUN STAN MODEL ---
stan_data <- list(
  N = n,
  x = x,
  y = y,
  K = K_basis,
  B = B_matrix
)

# Compile and sample
# Note: Increase 'iter' for production runs
fit_bayesian <- stan(
  file = 'megpd/bmegpd.stan',
  data = stan_data,
  iter = 2000,
  chains = 4,
  cores = 4
)

mcmc_trace(fit_bayesian, pars = c("kappa", "sigma", "xi", "tau"))
mcmc_trace(fit_bayesian, pars = "delta")

# --- 4. POSTERIOR ANALYSIS ---
# Extract parameters
print(fit_bayesian, pars = c("kappa", "sigma", "xi", "tau"))
print(fit_bayesian, pars = "beta")

# Visualize the fitted delta(R)
post_beta <- extract(fit_bayesian)$beta
# Calculate the mean function delta(R)
# Note: This is an average over all posterior spline shapes
mean_beta <- colMeans(post_beta)
log_delta_fit <- B_matrix %*% mean_beta
delta_fit <- exp(log_delta_fit) + 0.01

plot(R_obs, delta_fit, col = "blue", pch = 16, cex = 0.5,
     xlab = "Radius R", ylab = "delta(R)", main = "Bayesian Heteroscedastic Fit")

# --- 1. Prepare a Smooth Grid for Plotting ---
# Create a sequence of R values from min to max observed
R_grid <- seq(min(R_obs), max(R_obs), length.out = 500)

# IMPORTANT: The basis functions must be identical to the ones used in fitting.
# We use the 'predict' style for splines by using the original B_matrix attributes.
# Or simpler: re-generate the basis on the log scale to match the Stan code.
B_grid <- bs(log(R_grid), df = K_basis, degree = 3, intercept = TRUE)

# --- 2. Calculate delta(R) for every Posterior Sample ---
# post_beta is a matrix [samples x K_basis]
# B_grid is a matrix [500 x K_basis]
# We use matrix multiplication to get [samples x 500 points]
log_delta_samples <- post_beta %*% t(B_grid)
delta_samples <- exp(log_delta_samples) + 0.01

# --- 3. Calculate Summary Statistics (Mean and Quantiles) ---
delta_mean  <- colMeans(delta_samples)
delta_lower <- apply(delta_samples, 2, quantile, probs = 0.025)
delta_upper <- apply(delta_samples, 2, quantile, probs = 0.975)

# --- 4. Plot with Credible Bands ---
# Initialize the plot
plot(R_grid, delta_mean, type = "n", # "n" means don't plot yet
     ylim = c(0, max(delta_upper) * 1.1),
     xlab = "Radius R", ylab = expression(delta(R)),
     main = "Bayesian delta(R) with 95% Credible Band")

# Add the shaded credible band
polygon(c(R_grid, rev(R_grid)), 
        c(delta_lower, rev(delta_upper)), 
        col = rgb(0, 0, 1, 0.2), border = NA)

# Add the posterior mean line
lines(R_grid, delta_mean, col = "blue", lwd = 2)

# Optional: Add the original data points at the bottom for reference (rug plot)
rug(R_obs, col = rgb(0,0,0,0.1))

