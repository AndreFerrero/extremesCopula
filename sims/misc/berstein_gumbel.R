# ==============================================================================
# SETUP: Libraries, Sources, and Shared Parameters
# ==============================================================================
library(ExtremalDep)
library(copula)
library(rstan)

# Load custom model components
source("code/models/margins/egp.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/sim_copula_markov.R")
source("code/models/copula_markov/copula_markov_model.R")

rstan_options(auto_write = TRUE)

# Simulation & Model Constants
EGPD_PARAM  <- c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2)
THETA_TRUE  <- 2
HYPERPARAM  <- list(mu.nbinom = 3.2, var.nbinom = 4.48)
NSIM_MCMC   <- 30000
BURN        <- 5000

# Reusable Plotting Function for Spaghetti Plots
plot_nonparam_h <- function(summary_obj, theta_true, main_title) {
  # 1. Define truth
  # Formula: h(w) = 0.5 * (theta-1) * (w^theta + (1-w)^theta)^(1/theta - 2)
  # Simplified for plotting on the summary object grid
  h_true_vals <- 0.5 * (theta_true - 1) * 
                 (summary_obj$w^2 + (1 - summary_obj$w)^2)^(1 / theta_true - 2)
  
  # 2. Setup Plot
  plot(summary_obj$w, summary_obj$h.mean, type = "n", ylim = c(0, 2),
       xlab = expression(omega), ylab = expression(h(omega)), main = main_title)
  
  # 3. Spaghetti lines
  thin_idx <- seq(1, nrow(summary_obj$h_post), length.out = 300)
  for (i in thin_idx) {
    lines(summary_obj$w, summary_obj$h_post[i, ], col = rgb(0.7, 0.7, 0.7, 0.1))
  }
  
  # 4. Credibility Interval Band
  polygon(c(summary_obj$w, rev(summary_obj$w)),
          c(summary_obj$h.low, rev(summary_obj$h.up)),
          col = rgb(1, 0, 0, 0.2), border = NA)
  
  # 5. Mean and Truth Overlays
  lines(summary_obj$w, summary_obj$h.mean, col = "red", lwd = 2)
  lines(summary_obj$w, h_true_vals, col = "blue", lwd = 2, lty = 2)
  
  legend("top", legend = c("Post. Samples", "95% Cred. Int.", "Post. Mean", "True Gumbel"),
         col = c("grey", rgb(1, 0, 0, 0.2), "red", "blue"),
         lty = c(1, 1, 1, 2), lwd = c(1, 10, 2, 2), bty = "n")
}

# ==============================================================================
# CASE 1: THINNED MARKOV CHAIN DATA
# ==============================================================================
cat("\n--- Running Case 1: Markov Chain Analysis ---\n")

# Compile and Simulate Markov Model
gumbel_stan <- rstan::stan_model("code/stan/gumbel_egpd.stan")
gumbel_model <- make_copula_markov_model(margin = margin_egp, copula = copula_gumbel, stan_mod = gumbel_stan)

gumbel_sim <- gumbel_model$simulate(n = 10000, copula_param = THETA_TRUE, 
                                    margin_param = EGPD_PARAM, seed = 46)

# Extract Thinned Pairs (Skip-lag k=20)
k_skip <- 20
idx <- seq(1, length(gumbel_sim$x) - 1, by = k_skip)
pairs_markov <- cbind(gumbel_sim$x[idx], gumbel_sim$x[idx + 1])

# Fit Non-parametric Model
res_markov <- fExtDep.np(x = pairs_markov, method = "Bayesian", mar.fit = TRUE, 
                         u = TRUE, par10 = rep(0.1, 3), par20 = rep(0.1, 3),
                         sig10 = 0.0001, sig20 = 0.0001, k0 = 5, 
                         hyperparam = HYPERPARAM, nsim = NSIM_MCMC)

# Summarize and Plot
sum_markov <- summary(res_markov, burn = BURN)
plot_nonparam_h(sum_markov, THETA_TRUE, "Angular Density: Thinned Markov Pairs")

# ==============================================================================
# CASE 2: ACTUAL i.i.d. GUMBEL SAMPLES
# ==============================================================================
cat("\n--- Running Case 2: i.i.d. Gumbel Analysis ---\n")

# Generate True i.i.d. Gumbel realisations
gumbel_cop <- gumbelCopula(param = THETA_TRUE, dim = 2)
u_iid <- rCopula(gumbel_cop, n = 1000)

# Quantile transform to EGPD margins
x_iid <- apply(u_iid, 2, margin_egp$quantile, param = EGPD_PARAM)

# Fit Non-parametric Model
res_iid <- fExtDep.np(x = x_iid, method = "Bayesian", mar.fit = TRUE, 
                      u = TRUE, par10 = rep(0.1, 3), par20 = rep(0.1, 3),
                      sig10 = 0.0001, sig20 = 0.0001, k0 = 5, 
                      hyperparam = HYPERPARAM, nsim = NSIM_MCMC)

# Summarize and Plot
sum_iid <- summary(res_iid, burn = BURN)
plot_nonparam_h(sum_iid, THETA_TRUE, "Angular Density: i.i.d. Gumbel Samples")

# Compare diagnostics if needed
# diagnostics(res_markov)
# diagnostics(res_iid)