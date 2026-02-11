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
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/copulas/gaussian.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/builders/markov/sim_copula_markov.R")

mod_sim <- simulate_copula_markov

# --- 3. GROUND TRUTH & DATA GENERATION ---

set.seed(46)

param_egp <- c(kappa = 1.5, sigma = 2.0, xi = 0.2)
rho_true <- 0.7
n_obs <- 1000

# Generate simulated "observed" data
sim_obs <- mod_sim(n = n_obs, copula = copula_gaussian, theta = rho_true, margin = margin_egp, margin_param = param_egp)
X_sim <- sim_obs$X

plot(1:n_obs, X_sim, type = "l", main = "Simulated Extremal Time Series", col = "darkblue")
hist(X_sim, breaks = 50)

stan_mod <- stan_model("C:/Users/Andrea Ferrero/extremesCopula/code/stan/gumbel_egpd.stan")

stan_data <- list(T = length(X_sim), x = X_sim, unif_prior = 0, prior_check = 0, run_ppc = 1, I = 25)

fit_wrong <- sampling(stan_mod, data = stan_data, iter = 2000, chains = 4, seed= 42)

# Check the estimated theta
print(fit_wrong, pars = "theta")

# --- 8. PHASE 4: POSTERIOR PREDICTIVE CHECKS (PPC) ---

y_rep <- rstan::extract(fit_wrong)$x_rep

ppc_dens_overlay(X_sim, y_rep[1:100, ]) +
  ggtitle("Posterior Predictive Check: Densities") +
    coord_cartesian(xlim = c(0, 10))

ppc_stat(X_sim, y_rep, stat = "max") +
  ggtitle("Posterior Predictive Check: Maximum Values")

# Extremogram PPC (Temporal Dependence Validation)
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

obs_ext <- calc_extremogram(X_sim)
rep_exts <- t(apply(y_rep[1:100, ], 1, calc_extremogram))

plot_df <- data.frame(
  lag = 1:10, obs = obs_ext,
  low = apply(rep_exts, 2, quantile, 0.025),
  high = apply(rep_exts, 2, quantile, 0.975)
)

ggplot(plot_df, aes(x = lag)) +
  # Uncertainty Ribbon
  geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  # Observed Data Line
  geom_line(aes(y = obs), color = "red", linewidth = 1) +
  # Points to emphasize each discrete lag
  geom_point(aes(y = obs), color = "red") +
  # Force integer breaks on the x-axis
  scale_x_continuous(breaks = 1:10) + 
  labs(
    title = "Extremogram PPC", 
    subtitle = "95% Model Credible Interval",
    x = "Lag (Steps in Time)",
    y = "Tail Dependence (rho)"
  ) +
  theme_minimal()

  
# --- 11. RETURN LEVEL ASSESSMENT ---

# Define the Return Periods of interest (in number of observations)
# If daily data: 1 year = 365, 10 years = 3650, 100 years = 36500
return_periods <- c(100, 500, 1000, 5000, 10000)

# Extract posterior draws
draws <- as.data.frame(fit_wrong, pars = c("kappa", "sigma", "xi"))

# Function to calculate EGPD Quantile
# (Internal consistency with your margin_egp object)
calc_egpd_quantile <- function(p, kappa, sigma, xi) {
  p_gpd <- p^(1 / kappa)
  # Standard GPD quantile formula
  if (abs(xi) < 1e-10) {
    return(-sigma * log(1 - p_gpd))
  } else {
    return((sigma / xi) * ((1 - p_gpd)^(-xi) - 1))
  }
}

# Calculate Return Levels for every MCMC draw
rl_results <- map_df(return_periods, function(m) {
  p_target <- 1 - (1 / m)

  # Compute level for each draw
  levels <- mapply(calc_egpd_quantile,
    MoreArgs = list(p = p_target),
    kappa = draws$kappa,
    sigma = draws$sigma,
    xi = draws$xi
  )

  # Compute True Return Level (from true_params)
  true_level <- calc_egpd_quantile(
    p_target,
    true_params$kappa,
    true_params$sigma,
    true_params$xi
  )

  data.frame(
    m = m,
    mean_rl = mean(levels),
    low_rl = quantile(levels, 0.025),
    high_rl = quantile(levels, 0.975),
    true_rl = true_level
  )
})

print(rl_results)

# --- 12. VISUALIZE RETURN LEVEL PLOT ---
ggplot(rl_results, aes(x = m)) +
  geom_ribbon(aes(ymin = low_rl, ymax = high_rl), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = mean_rl, color = "Posterior Mean"), size = 1) +
  geom_line(aes(y = true_rl, color = "True Value"), linetype = "dashed", size = 1) +
  scale_x_log10(breaks = return_periods) +
  labs(
    title = "Return Level Plot with 95% Credible Intervals",
    x = "Return Period (m observations)",
    y = "Return Level (z_m)",
    color = "Legend"
  ) +
  theme_minimal()