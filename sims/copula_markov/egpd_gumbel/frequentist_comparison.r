# ==============================================================================
# SIMULATION STUDY: THE IMPACT OF DEPENDENCE AND THRESHOLD SELECTION
# ==============================================================================
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")
source("winter2016/marg_dep_mods_funcs.R")

library(dplyr)
library(ggplot2)

sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

lag_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) |>
    na.omit()

  return(df)
}

# 1. Setup True Parameters (Based on ERA5 results)
true_margin_par <- c(mu = 0, kappa = 6, sigma = 1, xi = 0)
true_theta <- 1
n_sim <- 10000

# 2. Generate Synthetic "Dependent" Wind Data
synthetic_data <- sim_model$simulate(
  n = n_sim,
  copula_param = true_theta,
  margin_param = true_margin_par,
  seed = 123
)$x

threshrange.plot(synthetic_data, nint = 50)
quantile(synthetic_data, seq(0.9, 0.99, by = 0.01))

# ------------------------------------------------------------------------------
# POINT 1: PROVING THAT IGNORING DEPENDENCE BIASES XI
# ------------------------------------------------------------------------------

# Fit A: Independent (IID) EGPD Model (Ignoring Copula)
fit_iid_egpd <- egpd::fitegpd(synthetic_data, type = 1)
xi_iid <- fit_iid_egpd$estimate["xi"]

# Fit B: Independent (IID) GPD Model (Standard Threshold Approach)
u_fix <- quantile(synthetic_data, 0.98)
fit_iid_gpd <- extRemes::fevd(synthetic_data, threshold = u_fix, type = "GP")
xi_gpd_iid <- fit_iid_gpd$results$par["shape"]

# Fit C: Your Proposition (Joint EGPD + Gumbel)
synthetic_lag <- lag_df(synthetic_data)

tau <- cor(synthetic_lag$x, synthetic_lag$x_lag, method = "kendall")

theta_tau <- 1 / (1 - tau)
theta_tau

# Initial values from marginal-only fit
init_marg <- fitegpd(synthetic_data, method = "mle", model = 1)

# Initial theta_vec
# [log_kappa, log_sigma, xi, log(theta_copula - 1)]
start_params <- c(
  log(init_marg$estimate["kappa"]),
  log(init_marg$estimate["sigma"]),
  init_marg$estimate["xi"],
  log(theta_tau - 1)
)

fit_res_egpd_gumbel <- fit_egpd_gumbel(synthetic_data, init_par = start_params)

res_kappa <- exp(fit_res_egpd_gumbel$par[1])
res_sigma <- exp(fit_res_egpd_gumbel$par[2])
res_xi <- fit_res_egpd_gumbel$par[3]
res_theta <- exp(fit_res_egpd_gumbel$par[4]) + 1

cat("--- XI Bias Analysis ---\n")
cat("True Xi:      ", true_margin_par["xi"], "\n")
cat("GUMBEL EGPD:  ", res_xi, "\n")
cat("IID EGPD:     ", xi_iid, "\n")
cat("IID GPD:      ", xi_gpd_iid, "\n")

# ------------------------------------------------------------------------------
# POINT 2: PROVING STABILITY (THRESHOLD-FREE vs. CENSORED)
# ------------------------------------------------------------------------------

# We will vary a "cut-off" point from 85% to 98%
u_sequence <- seq(0.9, 0.98, by = 0.01)
stab_results <- data.frame(prob = u_sequence, xi_censored = NA, xi_egpd = NA)


for (i in seq_along(u_sequence)) {
  u_val <- quantile(synthetic_data, u_sequence[i])

  # 1. Censored Gumbel-GPD (Winter2016 style)
  fit_cens <- FitGpd(dat = lag_df(synthetic_data), u = u_val, optim.type = 2)
  stab_results$xi_censored[i] <- fit_cens$par[3]

  # 2. Threshold-Free EGPD-Gumbel
  stab_results$xi_egpd[i] <- res_xi
}

stab_results |>
  tidyr::pivot_longer(cols = c(xi_censored, xi_egpd), names_to = "model", values_to = "xi") |>
  ggplot(aes(x = prob, y = xi, color = model)) +
  geom_line() +
  geom_hline(yintercept = true_margin_par["xi"], linetype = "dashed", color = "red") +
  labs(
    title = "Threshold Stability Comparison",
    subtitle = paste("True xi =", true_margin_par["xi"], "| N =", n_sim),
    x = "Threshold Quantile (u)",
    y = "Estimated xi",
    color = "Model"
  ) +
  theme_minimal()


# ==============================================================================
# Monte Carlo Simulation: Sensitivity to Dependence Strength (Theta)
# ==============================================================================

# True marginal parameters
true_margin_par <- c(mu = 0, kappa = 6.0, sigma = 2.5, xi = 0.1)
theta_sequence <- c(1.5, 3.0, 4.5, 6.0) # From weak to ERA5-level dependence
n_sim <- 3500
mc_it <- 50 # Number of MC trials per theta (increase for final results)

# Storage for results
mc_results <- data.frame()

for (th in theta_sequence) {
  message(paste("Running simulations for Theta =", th))

  for (i in 1:mc_it) {
    # 1. Simulate data with current theta
    sim_data <- sim_model$simulate(
      n = n_sim,
      copula_param = th,
      margin_param = true_margin_par
    )$x

    # 2. Fit IID Model (The "Naïve" approach)
    fit_egpd_iid <- egpd::fitegpd(sim_data, type = 1)
    xi_egpd_iid <- fit_egpd_iid$estimate["xi"]

    u_fix <- quantile(sim_data, 0.90)
    fit_gpd_iid <- extRemes::fevd(sim_data, threshold = u_fix, type = "GP")
    xi_gpd_iid <- fit_gpd_iid$results$par["shape"]

    # 3. Fit Gumbel dependence +  EGPD
    synthetic_lag <- lag_df(sim_data)

    # Initialise theta
    tau <- cor(synthetic_lag$x, synthetic_lag$x_lag, method = "kendall")

    theta_tau <- 1 / (1 - tau)
    theta_tau

    # Initial theta_vec
    # [log_kappa, log_sigma, xi, log(theta_copula - 1)]
    start_params <- c(
      log(fit_egpd_iid$estimate["kappa"]),
      log(fit_egpd_iid$estimate["sigma"]),
      fit_egpd_iid$estimate["xi"],
      log(theta_tau - 1)
    )

    fit_res_egpd_gumbel <- fit_egpd_gumbel(sim_data, init_par = start_params)
    xi_egpd_gumbel <- fit_res_egpd_gumbel$par[3]

    # Fit Censored GPD + Gumbel (Winter2016 style)
    fit_gpd_gumbel_cens <- FitGpd(dat = synthetic_lag, u = u_fix, optim.type = 2)
    xi_gpd_gumbel_cens <- fit_gpd_gumbel_cens$par[3]

    # 4. Store
    mc_results <- rbind(mc_results, data.frame(
      theta = th,
      iteration = i,
      model = "IID_EGPD",
      xi_hat = xi_egpd_iid
    ))
    mc_results <- rbind(mc_results, data.frame(
      theta = th,
      iteration = i,
      model = "IID_GPD",
      xi_hat = xi_gpd_iid
    ))
     mc_results <- rbind(mc_results, data.frame(
      theta = th,
      iteration = i,
      model = "Censored_GPD_GUMBEL",
      xi_hat = xi_gpd_gumbel_cens
    ))
    mc_results <- rbind(mc_results, data.frame(
      theta = th,
      iteration = i,
      model = "EGPD_GUMBEL",
      xi_hat = xi_egpd_gumbel
    ))
  }
}

# save(mc_results, file = "sims/copula_markov/egpd_gumbel/mc_dependence_xi03.RData")
load("sims/copula_markov/egpd_gumbel/mc_dependence_xi01.RData")

# ==============================================================================
# VISUALIZATION OF THE "DEPENDENCE DEGRADATION"
# ==============================================================================
library(ggplot2)

ggplot(mc_results, aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = true_margin_par["xi"], linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True xi =", true_margin_par["xi"], "| N =", n_sim, "| MC Iterations =", mc_it),
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()

mc_results |>
  filter(model != "Censored_GPD_GUMBEL") |>
  ggplot(aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = true_margin_par["xi"], linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True xi =", true_margin_par["xi"], "| N =", n_sim, "| MC Iterations =", mc_it),
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()

mc_results |>
  filter(model == "Censored_GPD_GUMBEL" | model == "EGPD_GUMBEL") |>
  ggplot(aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = true_margin_par["xi"], linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True xi =", true_margin_par["xi"], "| N =", n_sim, "| MC Iterations =", mc_it),
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()


# Calculate median and MAD for each theta-model combination

summary_stats <- mc_results |>
  group_by(theta, model) |>
  summarise(
    median_xi = median(xi_hat),
    mad_xi = mad(xi_hat),
    .groups = "drop"
  )

ggplot(summary_stats, aes(x = as.factor(theta), y = median_xi, color = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = median_xi - mad_xi, ymax = median_xi + mad_xi),
    width = 0.2,
    linewidth = 0.6,
    position = position_dodge(width = 0.7)
  ) +
  geom_hline(yintercept = true_margin_par["xi"], linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True xi =", true_margin_par["xi"], "| N =", n_sim, "| MC Iterations =", mc_it, "| Error bars: MAD"),
    x = "Dependence Strength (Theta)",
    y = "Estimated xi (Median)",
    color = "Model"
  ) +
  theme_minimal()

