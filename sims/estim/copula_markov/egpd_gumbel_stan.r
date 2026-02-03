library(rstan)
library(ggplot2)
library(bayesplot)
library(copula)

# --- 1. SETTINGS & SIMULATION ---
margin_egp <- list(
  name = "egp",

  G_dist = function(u, param) u^param["kappa"],
  G_inv  = function(u, param) u^(1 / param["kappa"]),
  g_dist = function(u, param) {
    param["kappa"] * u^(param["kappa"] - 1)
  },

  cdf = function(x, param) {
    u <- evd::pgpd(x / param["sigma"], shape = param["xi"])
    margin_egp$G_dist(u, param)
  },

  lpdf = function(x, param) {
    u <- evd::pgpd(x / param["sigma"], shape = param["xi"])

    log(margin_egp$g_dist(u, param)) +
      log(evd::dgpd(x / param["sigma"], shape = param["xi"])) -
      log(param["sigma"])
  },

  quantile = function(p, param) {
    param["sigma"] * evd::qgpd(
      margin_egp$G_inv(p, param),
      shape = param["xi"]
    )
  },

  sample = function(n, param) {
    U <- runif(n)
    margin_egp$quantile(U, param)
  }
)


# 1. The Gumbel h-function (Conditional CDF)
# We use this to solve: h(u | v, theta) = w
gumbel_hfunc <- function(u, v, theta) {
  # Numerical clamping to prevent log(0) or log(1)
  u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
  v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
  
  ln_u <- -log(u)
  ln_v <- -log(v)
  term <- ln_u^theta + ln_v^theta
  C <- exp(-term^(1/theta))
  
  # Conditional probability P(U <= u | V = v)
  h <- C * (ln_v^(theta - 1)) * (term^(1/theta - 1)) / v
  return(h)
}

# 2. Optimized Markov Simulation using Bisection
simulate_gumbel_markov_egpd <- function(n, theta, param_egp, burn = 100, bisec_it = 20) {
  
  total_n <- n + burn
  U <- numeric(total_n)
  
  # Initial state
  U[1] <- runif(1, 1e-5, 1 - 1e-5)
  
  for (t in 2:total_n) {
    v_prev <- U[t - 1]
    w_target <- runif(1) # The random uniform to invert
    
    # Bisection Solver (replaces the old cCopula logic)
    low <- 1e-10
    high <- 1 - 1e-10
    
    for (i in 1:bisec_it) {
      mid <- (low + high) / 2
      if (gumbel_hfunc(mid, v_prev, theta) < w_target) {
        low <- mid
      } else {
        high <- mid
      }
    }
    U[t] <- (low + high) / 2
  }
  
  # Remove burn-in and transform to EGPD scale
  U_final <- U[(burn + 1):total_n]
  X <- margin_egp$quantile(U_final, param_egp)
  
  return(list(
    X = X,
    U = U_final,
    n_failures = 0 # With Bisection, this is always 0
  ))
}

set.seed(123)
true_params <- list(kappa = 1.5, sigma = 2.0, xi = 0.1, theta = 2.0)
param_egp <- c(kappa = 1.5, sigma = 2.0, xi = 0.1)

n_obs <- 1000

sim_obs <- simulate_gumbel_markov_egpd(n_obs, true_params[["theta"]], param_egp)
X_sim <- sim_obs$X

plot(1:n_obs, X_sim, type = "l")

# --- 2. COMPILE & FIT STAN MODEL ---
# mod <- stan_model("C:/Users/Andrea Ferrero/extremesCopula/sims/estim/copula_markov/gumbel_egpd.stan")

# # Running the sampler
# stan_fit <- sampling(
#     mod,
#     data = stan_data,
#     iter = 2000,
#     chains = 4,
#     cores = 4,
#     seed = 42,
#     control = list(adapt_delta = 0.98, max_treedepth = 15)
# )

mod_ppc <- stan_model("C:/Users/Andrea Ferrero/extremesCopula/sims/estim/copula_markov/gumbel_egpd_ppc.stan")

# --- 1. RUN PRIOR CHECK ONLY ---
stan_data_prior <- list(T = n_obs, x = X_sim, prior_check = 1)
fit_prior <- sampling(mod_ppc, data = stan_data_prior, iter = 1000)

y_prior <- rstan::extract(fit_prior)$x_rep

# --- 1. Distribution Overlay ---
# Does the data (X_sim) look like a plausible draw from the prior?
ppc_dens_overlay(X_sim, y_prior[1:50, ]) + 
  scale_x_log10() +
  labs(title = "Prior Predictive Overlay",
       subtitle = "Observed data (red) vs. 50 datasets generated from Priors only",
       x = "Value (Log Scale)") +
  theme_minimal()

# --- 2. Extremes Check ---
# Does the prior allow for the maximum value we actually observed?
ppc_stat(X_sim, y_prior, stat = "max") +
  labs(title = "Prior Predictive: Maximum Value",
       subtitle = "Where does the observed max fall in the 'Prior World'?")

ppc_ribbon(X_sim, y_prior) +
  labs(title = "Prior Predictive: 95% Uncertainty Ribbon",
       subtitle = "Observed data (red) vs. Range of possibilities from Priors",
       x = "Time", y = "Value") +
  theme_minimal()

# RUN POSTERIOR

stan_data <- list(T = length(X_sim), x = X_sim, prior_check = 0)

stan_fit_ppc <- sampling(
    mod_ppc,
    data = stan_data,
    iter = 2000,
    chains = 4,
    cores = 4,
    seed = 42,
    control = list(adapt_delta = 0.98, max_treedepth = 15)
)

# --- 3. POST-ESTIMATION ANALYSIS (The stanfit Object) ---

# A. Standard Print (Check Rhat and ESS)
print(stan_fit_ppc, pars = c("kappa", "sigma", "xi", "theta"))

# B. Parameter Recovery Visualization
mcmc_recover_intervals(
    as.array(stan_fit_ppc, pars = c("kappa", "sigma", "xi", "theta")),
    true = unlist(true_params)
) + ggtitle("Posterior Recovery: Red marks true value")

# C. Traceplots (Check for "hairy caterpillars" and no divergences)
color_scheme_set("brightblue")
mcmc_trace(stan_fit_ppc, pars = c("kappa", "sigma", "xi", "theta"))

# D. Correlation between parameters (Very important for PhD work)
# This shows if kappa and sigma are "fighting" each other
pairs(stan_fit_ppc, pars = c("kappa", "sigma", "xi", "theta"))


## Posterior Predictive Checking

y_rep <- rstan::extract(stan_fit_ppc)$x_rep # This is a matrix [samples x T]

# 2. Density Overlay (Does the model capture the distribution?)
# We take a subset of 50 samples to keep the plot clean
ppc_dens_overlay(X_sim, y_rep[1:50, ]) + 
  ggtitle("PPC: Observed vs Replicated Densities")

# 3. Maximum Value Check (CRITICAL for Extremes)
# Does our model generate maximums as high as the observed data?
ppc_stat(X_sim, y_rep, stat = "max") + 
  ggtitle("PPC: Distribution of Maximum Values")

# --- 4. PRIOR VS POSTERIOR OVERLAP ---
library(tidyr)
library(dplyr)

# 1. Extract posterior samples into a data frame
post_samples <- as.data.frame(stan_fit_ppc, pars = c("kappa", "sigma", "xi", "theta")) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "posterior")

# 2. Define theoretical prior densities based on your Stan model
prior_data <- post_samples %>%
  group_by(parameter) %>%
  # Create a sequence that goes slightly wider than the posterior to see the "tails"
  summarise(min_val = min(posterior) * 0.8, max_val = max(posterior) * 1.2) %>%
  rowwise() %>%
  mutate(x_seq = list(seq(min_val, max_val, length.out = 400))) %>%
  unnest(x_seq) %>%
  mutate(prior = case_when(
    parameter == "kappa" ~ dlnorm(x_seq, 0, 2),
    parameter == "sigma" ~ dlnorm(x_seq, 0, 2),
    parameter == "xi"    ~ dnorm(x_seq, 0, 0.5),
    
    # Corrected Shifted Exponential logic
    parameter == "theta" ~ ifelse(x_seq >= 1, dexp(x_seq - 1, 0.1), 0),
    
    TRUE ~ NA_real_
  ))

# 3. Create the Plot
ggplot() +
  # Prior distribution (Area)
  geom_area(data = prior_data, aes(x = x_seq, y = prior, fill = "Prior"), alpha = 0.3) +
  geom_line(data = prior_data, aes(x = x_seq, y = prior), color = "black", linetype = "dashed") +
  
  # Posterior distribution (Density)
  geom_density(data = post_samples, aes(x = posterior, fill = "Posterior"), alpha = 0.7) +
  
  # True parameter values (Vertical lines)
  geom_vline(data = data.frame(
    parameter = c("kappa", "sigma", "xi", "theta"),
    true_val = unlist(true_params)
  ), aes(xintercept = true_val), color = "red", size = 1) +
  
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("Prior" = "grey70", "Posterior" = "blue")) +
  labs(
    title = "Prior-Posterior Overlap Analysis",
    subtitle = "Red lines indicate the 'True' values used in simulation",
    x = "Parameter Value",
    y = "Density",
    fill = "Distribution"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

calc_extremogram <- function(x, prob = 0.95, max_lag = 10) {
  u <- quantile(x, prob)
  n <- length(x)
  ext_vec <- numeric(max_lag)
  
  # Indicator: 1 if extreme, 0 otherwise
  is_ext <- x > u
  denom <- sum(is_ext)
  
  for (h in 1:max_lag) {
    # Count times both t and t+h are extreme
    num <- sum(is_ext[1:(n-h)] & is_ext[(h+1):n])
    ext_vec[h] <- num / denom
  }
  return(ext_vec)
}

# --- For your PPC Plot ---
obs_ext <- calc_extremogram(X_sim)
# Calculate extremograms for 100 posterior replicates
rep_exts <- t(apply(y_rep[1:100, ], 1, calc_extremogram))

# Plotting logic
library(tidyverse)
plot_df <- data.frame(
  lag = 1:10,
  obs = obs_ext,
  low = apply(rep_exts, 2, quantile, 0.025),
  high = apply(rep_exts, 2, quantile, 0.975)
)

ggplot(plot_df, aes(x = lag)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = obs), color = "red", size = 1) +
  geom_point(aes(y = obs), color = "red") +
  labs(title = "Extremogram PPC", 
       subtitle = "Red: Observed Data | Blue: 95% Model Credible Interval",
       y = "Conditional Prob of Extreme", x = "Lag (Time)") +
  theme_minimal()
