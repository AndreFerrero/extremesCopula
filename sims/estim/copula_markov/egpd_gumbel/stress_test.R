# --- 1. LIBRARIES & GLOBAL SETTINGS ---
library(rstan)
library(ggplot2)
library(bayesplot)
library(copula)
library(tidyr)
library(dplyr)
library(purrr)
library(evd)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

# --- 2. MATHEMATICAL ENGINE: EGPD & COPULA FUNCTIONS ---
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/margins/egp.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/margins/frechet.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/copulas/gumbel.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/builders/markov/sim_gumbel_markov.R")

gumbel_h <- copula_gumbel$h_dist

mod_sim <- simulate_gumbel_markov

# --- 3. GROUND TRUTH & DATA GENERATION ---

set.seed(46)
frechet_param <- c(scale = 2.0, shape = 4)
theta_true <- 2
n_obs <- 1000

# Generate simulated "observed" data
sim_frechet <- mod_sim(n = n_obs, copula_h = gumbel_h, theta = theta_true, margin = margin_frechet, margin_param = frechet_param)
X_frechet <- sim_frechet$X

plot(1:n_obs, X_frechet, type = "l", main = "Simulated Extremal Time Series", col = "darkblue")
hist(X_frechet)

stan_mod <- stan_model("C:/Users/Andrea Ferrero/extremesCopula/code/stan/gumbel_egpd_ppc.stan")

stan_data_frechet <- list(T = length(X_frechet), x = X_frechet, unif_prior = 0, prior_check = 0, run_ppc = 1, I = 25)

stan_fit_frechet <- sampling(
    stan_mod,
    data = stan_data_frechet,
    iter = 2000, chains = 4, cores = 4, seed = 42
)

print(stan_fit_frechet, pars = c("mu", "kappa", "sigma", "xi", "theta"))

mcmc_trace(stan_fit_frechet, pars = c("mu", "kappa", "sigma", "xi", "theta"))
mcmc_acf(stan_fit_frechet, pars = c("mu", "kappa", "sigma", "xi", "theta"))
pairs(stan_fit_frechet, pars = c("mu", "kappa", "sigma", "xi", "theta"))
mcmc_neff_hist(neff_ratio(stan_fit_frechet))
mcmc_rhat_hist(rhat(stan_fit_frechet))

y_rep <- rstan::extract(stan_fit_frechet)$x_rep

ppc_dens_overlay(X_frechet, y_rep[1:200, ]) +
    ggtitle("Posterior Predictive Check: Densities")

ppc_stat(X_frechet, y_rep, stat = "max") +
    ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(X_frechet, y_rep, stat = "min") +
    ggtitle("Posterior Predictive Check: Minimum Values")

# Extremogram Calculation (Lag-dependent Tail dependence)
calc_extremogram <- function(x, prob = 0.95, max_lag = 10) {
  u <- quantile(x, prob)
  n <- length(x)
  ext_vec <- numeric(max_lag)
  is_ext <- x > u
  denom <- sum(is_ext)
  for (h in 1:max_lag) {
    num <- sum(is_ext[1:(n - h)] & is_ext[(h + 1):n])
    ext_vec[h] <- num / denom
  }
  return(ext_vec)
}
# Extremogram PPC (Temporal Dependence Validation)
obs_ext <- calc_extremogram(X_frechet)
rep_exts <- t(apply(y_rep[1:100, ], 1, calc_extremogram))

plot_df <- data.frame(
  lag = 1:10, obs = obs_ext,
  low = apply(rep_exts, 2, quantile, 0.025),
  high = apply(rep_exts, 2, quantile, 0.975)
)

ggplot(plot_df, aes(x = lag)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = obs), color = "red", size = 1) +
  labs(title = "Extremogram PPC", subtitle = "95% Model Credible Interval") +
  theme_minimal()

# --- 11. RETURN LEVEL ASSESSMENT (STRESS TEST: EGPD FIT TO FRECHET DATA) ---

# Define return periods (m observations)
return_periods <- c(100, 500, 1000, 5000, 10000)

# Extract posterior draws from Stan
# Note: We now include 'mu' because your Stan model estimated a location parameter
draws_df <- as.data.frame(stan_fit_frechet, pars = c("mu", "kappa", "sigma", "xi"))

# Calculate Return Levels for every MCMC draw
rl_results <- map_df(return_periods, function(m) {
  p_target <- 1 - (1 / m)

  # Compute level for each posterior draw using the margin_egp object
  # row_params will contain c(mu, kappa, sigma, xi)
  levels <- apply(draws_df, 1, function(row_params) {
    # margin_egp$quantile calculates the excess; we add the estimated mu
    row_params["mu"] + margin_egp$quantile(p_target, row_params)
  })

  # Compute the GROUND TRUTH level using the TRUE Frechet parameters
  # This uses the margin_frechet object you sourced earlier
  # true_params was list(scale = 2.0, shape = 4, theta = 2.0)
  true_level <- margin_frechet$quantile(p_target, frechet_param)

  data.frame(
    m = m,
    mean_rl = mean(levels),
    low_rl = quantile(levels, 0.025),
    high_rl = quantile(levels, 0.975),
    true_rl = true_level
  )
})

# --- 12. VISUALIZE RETURN LEVEL PLOT ---

ggplot(rl_results, aes(x = m)) +
  # Credible Interval Ribbon
  geom_ribbon(aes(ymin = low_rl, ymax = high_rl, fill = "95% Credible Interval"), alpha = 0.2) +
  
  # Posterior Mean Line (The model's estimate)
  geom_line(aes(y = mean_rl, color = "EGPD Posterior Mean"), size = 1) +
  
  # True Value Line (The Fréchet reality)
  geom_line(aes(y = true_rl, color = "True Frechet Value"), linetype = "dashed", size = 1) +
  
  # Logarithmic scale for rarity
  scale_x_log10(breaks = return_periods, labels = scales::comma) +
  
  # Styling the legend and colors
  scale_color_manual(values = c("EGPD Posterior Mean" = "red", "True Frechet Value" = "blue")) +
  scale_fill_manual(values = c("95% Credible Interval" = "blue")) +
  
  labs(
    title = "Stress Test: Marginal Robustness Check",
    subtitle = "Fitting a Bayesian EGPD model to Fréchet-generated data",
    x = "Return Period (m observations, Log Scale)",
    y = "Return Level (z_m)",
    color = "Legend",
    fill = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")



########################
#Lognormal

source("C:/Users/Andrea Ferrero/extremesCopula/code/models/margins/lognormal.r")

# --- RUN TEST ---
set.seed(46)
lognorm_param <- c(mu = 1, sigma = 1)
sim_lognorm <- simulate_gumbel_markov(n = n_obs, copula_h = gumbel_h, theta = 2.0, 
                                   margin = margin_lognormal, margin_param = lognorm_param)
X_lognorm <- sim_lognorm$X
plot(1:length(X_lognorm), X_lognorm, type = "l")
hist(X_lognorm)
# Fit your EGPD Stan model to sim_lognorm$X

stan_data_lognorm <- list(T = length(X_lognorm), x = X_lognorm, unif_prior = 0, prior_check = 0, run_ppc = 1, I = 25)

stan_fit_lognorm <- sampling(
    stan_mod,
    data = stan_data_lognorm,
    iter = 2000, chains = 4, cores = 4, seed = 42
)

print(stan_fit_lognorm, pars = c("mu", "kappa", "sigma", "xi", "theta"))

mcmc_trace(stan_fit_lognorm, pars = c("mu", "kappa", "sigma", "xi", "theta"))
mcmc_acf(stan_fit_lognorm, pars = c("mu", "kappa", "sigma", "xi", "theta"))
pairs(stan_fit_lognorm, pars = c("mu", "kappa", "sigma", "xi", "theta"))
mcmc_neff_hist(neff_ratio(stan_fit_lognorm))
mcmc_rhat_hist(rhat(stan_fit_lognorm))

y_rep <- rstan::extract(stan_fit_lognorm)$x_rep

ppc_dens_overlay(X_frechet, y_rep[1:200, ]) +
    ggtitle("Posterior Predictive Check: Densities")

ppc_stat(X_frechet, y_rep, stat = "max") +
    ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(X_frechet, y_rep, stat = "min") +
    ggtitle("Posterior Predictive Check: Minimum Values")

# Extremogram Calculation (Lag-dependent Tail dependence)
calc_extremogram <- function(x, prob = 0.95, max_lag = 10) {
  u <- quantile(x, prob)
  n <- length(x)
  ext_vec <- numeric(max_lag)
  is_ext <- x > u
  denom <- sum(is_ext)
  for (h in 1:max_lag) {
    num <- sum(is_ext[1:(n - h)] & is_ext[(h + 1):n])
    ext_vec[h] <- num / denom
  }
  return(ext_vec)
}
# Extremogram PPC (Temporal Dependence Validation)
obs_ext <- calc_extremogram(X_frechet)
rep_exts <- t(apply(y_rep[1:100, ], 1, calc_extremogram))

plot_df <- data.frame(
  lag = 1:10, obs = obs_ext,
  low = apply(rep_exts, 2, quantile, 0.025),
  high = apply(rep_exts, 2, quantile, 0.975)
)

ggplot(plot_df, aes(x = lag)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = obs), color = "red", size = 1) +
  labs(title = "Extremogram PPC", subtitle = "95% Model Credible Interval") +
  theme_minimal()
