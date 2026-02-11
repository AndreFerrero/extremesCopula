###############################################################################
# PhD Research: Principled Bayesian Calibration (SBC)
# Model: EGPD-Gumbel Markov Chain
###############################################################################

# --- 1. LIBRARIES & ENGINES ---
library(SBC)
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(evd)

# Load your custom models
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/margins/egp.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/copulas/gumbel.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/builders/markov/sim_gumbel_markov.R")

# Settings for computation
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(future)
plan(multisession)
options(SBC.min_chunk_size = 5)

# --- 2. DEFINE THE GENERATOR FUNCTION ---
# This matches your priors and simulation logic exactly
egpd_gumbel_single_gen <- function(n_obs) {
  
  # A. Sample "True" Parameters from the Priors
  true_mu      <- rnorm(1, 5, 2.5)
  true_kappa   <- rlnorm(1, 2, 1)
  true_sigma   <- rexp(1, 0.1)
  
  # Focused xi sampling (Gumbel-Frechet domain, truncated at 0.5)
  true_xi <- 1
  while(true_xi > 0.5){
    true_xi <- rgamma(1, 2, 10)
  }
  
  true_thetam1 <- rgamma(1, 2, 1)
  true_theta   <- true_thetam1 + 1
  
  # B. Simulate a Markov Chain dataset from these parameters

  sim <- simulate_gumbel_markov(
    n = n_obs, 
    copula_h = copula_gumbel$h_dist, 
    theta = true_theta, 
    margin = margin_egp, 
    margin_param = c(kappa = true_kappa, sigma = true_sigma, xi = true_xi)
  )
  
  # C. Return the structure required by the SBC package
  return(list(
    # 'variables' must match the parameters block in Stan
    variables = list(
      mu = true_mu,
      kappa = true_kappa,
      sigma = true_sigma,
      xi = true_xi,
      thetam1 = true_thetam1
    ),
    # 'data' must match the data block in Stan
    generated = list(
      T = n_obs,
      x = as.vector(sim$X + true_mu),
      unif_prior = 0,
      prior_check = 0,
      run_ppc = 0,
      I = 25
    )
  ))
}

# --- 3. INITIALIZE SBC OBJECTS ---

# Create the generator object (N=100 iterations is robust for a PhD)
set.seed(46)
sbc_sims <- 50
sbc_gen <- SBC_generator_function(egpd_gumbel_single_gen, n_obs = 500)
sbc_data <- generate_datasets(sbc_gen, sbc_sims)

# 1. Compile the model once
stan_mod <- stan_model("C:/Users/Andrea Ferrero/extremesCopula/code/stan/gumbel_egpd_ppc.stan")

# 2. Define the backend specifically for rstan
sbc_backend <- SBC_backend_rstan_sample(
  stan_mod,
  warmup = 1000,
  iter = 2000,
  chains = 2,
  control = list(adapt_delta = 0.9)
)

# --- 4. EXECUTE CALIBRATION ---
sbc_results <- compute_SBC(sbc_data, sbc_backend, cache_mode = "results",
    cache_location = "C:/Users/Andrea Ferrero/extremesCopula/sims/estim/copula_markov/egpd_gumbel/res/sbc_egpd_gumbel_cache.rds")

# --- 5. ANALYSIS & VISUALIZATION ---

# A. Rank Histograms (The standard check)
# It automatically adds the red dashed line and the binomial uncertainty band
plot_rank_hist(sbc_results) +
  theme_minimal() +
  labs(title = "SBC Rank Histograms")

# B. ECDF Difference Plot (The 'Advanced' check)
# If the line stays within the grey region, the model is calibrated.
# This is much more sensitive than a histogram for small biases.
plot_ecdf(sbc_results) +
  theme_minimal() +
  labs(title = "ECDF Plot")

# C. Diagnostic Summary (Divergences and R-hat)
# High R-hat or many divergences will be flagged here
diag_summary <- summary(sbc_results)
print(diag_summary)

