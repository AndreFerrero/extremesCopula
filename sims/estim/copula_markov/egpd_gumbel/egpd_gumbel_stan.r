###############################################################################
# Copula-based Markov Models with EGPD Margins
# STAGE: Model Validation, Inference, and Calibration Workflow
###############################################################################

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
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/copulas/gumbel.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/builders/markov/sim_gumbel_markov.R")

gumbel_h <- copula_gumbel$h_dist

mod_sim <- simulate_gumbel_markov

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

# --- 3. GROUND TRUTH & DATA GENERATION ---

set.seed(46)
param_egp <- c(kappa = 1.5, sigma = 2.0, xi = 0.2)
theta_true <- 2.0
n_obs <- 1000

# Generate simulated "observed" data
sim_obs <- mod_sim(n = n_obs, copula_h = gumbel_h, theta = theta_true, margin = margin_egp, margin_param = param_egp)
X_sim <- sim_obs$X

plot(1:n_obs, X_sim, type = "l", main = "Simulated Extremal Time Series", col = "darkblue")

# --- 4. STAN MODEL COMPILATION ---

# Combined model for Prior/Posterior and PPC
# Replace the file path with your local path
stan_mod <- stan_model("C:/Users/Andrea Ferrero/extremesCopula/code/stan/gumbel_egpd_ppc.stan")

# --- 5. PHASE 1: PRIOR PREDICTIVE CHECKS (PPC) ---

stan_data_prior <- list(T = n_obs, x = X_sim, prior_check = 1, run_ppc = 0, I = 25)
fit_prior <- sampling(stan_mod, data = stan_data_prior, iter = 1000)
y_prior <- rstan::extract(fit_prior)$x_rep

ppc_dens_overlay(X_sim, y_prior[1:50, ]) +
  scale_x_log10() +
  labs(title = "Prior Predictive Overlay", subtitle = "Log Scale") + theme_minimal()

ppc_stat(X_sim, y_prior, stat = "max") +
  labs(title = "Prior Predictive Check: Maximum Value")

ppc_ribbon(X_sim, y_prior) +
  labs(title = "Prior Predictive: 95% Uncertainty Ribbon") + theme_minimal()

# --- 6. PHASE 2: POSTERIOR INFERENCE ---

stan_data_post_ppc <- list(T = length(X_sim), x = X_sim, unif_prior = 0, prior_check = 0, run_ppc = 1, I = 25)

stan_fit_ppc <- sampling(
  stan_mod,
  data = stan_data_post_ppc,
  iter = 3000, chains = 4, cores = 4, seed = 42,
  control = list(adapt_delta = 0.90, max_treedepth = 10)
)

# save(stan_fit_ppc, X_sim, true_params, param_egp, theta_true, file = "C:/Users/Andrea Ferrero/extremesCopula/sims/estim/copula_markov/egpd_gumbel/seed46_ppcfit_3kmcmc_4chains_xsim3k_unifpriors.Rdata")

# load("C:/Users/Andrea Ferrero/extremesCopula/sims/estim/copula_markov/egpd_gumbel/seed46_ppcfit_2kmcmc_4chains_xsim1k.Rdata")

# --- 7. PHASE 3: MCMC DIAGNOSTICS & RECOVERY ---

mcmc_recover_intervals(
  as.array(stan_fit_ppc, pars = c("kappa", "sigma", "xi", "theta")),
  true = unlist(true_params)
) + ggtitle("Parameter Recovery Check")

print(stan_fit_ppc, pars = c("mu", "kappa", "sigma", "xi", "theta"))

mcmc_trace(stan_fit_ppc, pars = c("mu", "kappa", "sigma", "xi", "theta"))
mcmc_acf(stan_fit_ppc, pars = c("mu", "kappa", "sigma", "xi", "theta"))
pairs(stan_fit_ppc, pars = c("mu", "kappa", "sigma", "xi", "theta"))
mcmc_neff_hist(neff_ratio(stan_fit_ppc))
mcmc_rhat_hist(rhat(stan_fit_ppc))

# --- 8. PHASE 4: POSTERIOR PREDICTIVE CHECKS (PPC) ---

y_rep <- rstan::extract(stan_fit_ppc)$x_rep

ppc_dens_overlay(X_sim, y_rep[1:100, ]) +
  ggtitle("Posterior Predictive Check: Densities")

ppc_stat(X_sim, y_rep, stat = "max") +
  ggtitle("Posterior Predictive Check: Maximum Values")

# Extremogram PPC (Temporal Dependence Validation)
obs_ext <- calc_extremogram(X_sim)
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

# --- 9. PHASE 5: SENSITIVITY - PRIOR VS POSTERIOR ---

post_samples <- as.data.frame(stan_fit_ppc, pars = c("kappa", "sigma", "xi", "theta")) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "posterior")

prior_data <- post_samples %>%
  group_by(parameter) %>%
  # We expand the sequence slightly to ensure we see the prior "box"
  summarise(min_val = min(posterior) * 0.5, max_val = max(posterior) * 1.5) %>%
  rowwise() %>%
  mutate(x_seq = list(seq(min_val, max_val, length.out = 400))) %>%
  unnest(x_seq) %>%
  mutate(prior = case_when(
    parameter == "kappa" ~ dunif(x_seq, 0.1, 10),
    parameter == "sigma" ~ dunif(x_seq, 0.1, 15),
    parameter == "xi"    ~ dunif(x_seq, -0.2, 0.5),
    # Fixed the ifelse syntax and the missing closing parenthesis
    parameter == "theta" ~ ifelse(x_seq >= 1 & x_seq <= 11, dunif(x_seq, 1, 11), 0),
    TRUE ~ NA_real_
  ))

# Plotting
ggplot() +
  # We use geom_line for priors if they are flat to make them more visible
  geom_area(data = prior_data, aes(x = x_seq, y = prior, fill = "Prior"), alpha = 0.3) +
  geom_density(data = post_samples, aes(x = posterior, fill = "Posterior"), alpha = 0.7) +
  geom_vline(
    data = data.frame(
      parameter = c("kappa", "sigma", "xi", "theta"),
      true_val = unlist(true_params)
    ),
    aes(xintercept = true_val), color = "red", linetype = "dashed", size = 1
  ) +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("Prior" = "grey70", "Posterior" = "blue")) +
  theme_minimal() +
  labs(title = "Prior-Posterior Sensitivity Analysis (Uniform Priors)",
       subtitle = "Flat grey line represents the objective Uniform prior")

# --- 11. RETURN LEVEL ASSESSMENT ---

# Define the Return Periods of interest (in number of observations)
# If daily data: 1 year = 365, 10 years = 3650, 100 years = 36500
return_periods <- c(100, 500, 1000, 5000, 10000)

# Extract posterior draws
draws <- as.data.frame(stan_fit_ppc, pars = c("kappa", "sigma", "xi"))

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

# --- 11. ANNUALIZED RETURN LEVEL ASSESSMENT ---

# Define observations per year (Daily = 365.25, Hourly = 8766)
obs_per_year <- 365.25

# Define return periods in YEARS
years_seq <- c(1, 2, 5, 10, 20, 50, 100, 200)
m_periods <- years_seq * obs_per_year
p_targets <- 1 - (1 / m_periods)

# Extract posterior draws
post_draws <- as.data.frame(stan_fit_ppc, pars = c("kappa", "sigma", "xi"))

# Calculate Return Levels for every MCMC draw
rl_matrix <- sapply(p_targets, function(p) {
  apply(post_draws, 1, function(row) {
    margin_egp$quantile(p, row)
  })
})

# Summary table
rl_summary <- data.frame(
  years = years_seq,
  mean_rl = colMeans(rl_matrix),
  median_rl = apply(rl_matrix, 2, median),
  low_rl = apply(rl_matrix, 2, quantile, probs = 0.025),
  high_rl = apply(rl_matrix, 2, quantile, probs = 0.975)
)

# Calculate True Return Level
true_param_vec <- c(kappa = true_params$kappa, sigma = true_params$sigma, xi = true_params$xi)
rl_summary$true_rl <- sapply(p_targets, function(p) {
  margin_egp$quantile(p, true_param_vec)
})

# --- 12. VISUALIZATION WITH YEAR LABELS ---
ggplot(rl_summary, aes(x = years)) +
  geom_ribbon(aes(ymin = low_rl, ymax = high_rl), fill = "royalblue", alpha = 0.2) +
  geom_line(aes(y = mean_rl, color = "Posterior Mean"), size = 1) +
  geom_line(aes(y = true_rl, color = "True Value"), linetype = "dashed", size = 1) +
  scale_x_log10(breaks = years_seq, labels = years_seq) +
  labs(
    title = "Return Level Plot: Annualized Risk",
    subtitle = "Extrapolation from Daily Process-based EGPD Model",
    x = "Return Period (Years)",
    y = "Return Level (Magnitude)",
    color = "Legend"
  ) +
  theme_minimal() +
  annotation_logticks(sides = "b") # Adds those small "ruler" lines to the log scale
