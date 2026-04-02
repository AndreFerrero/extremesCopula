# ==============================================================================
# SIMULATION STUDY: THE IMPACT OF DEPENDENCE AND THRESHOLD SELECTION
# ==============================================================================
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/margins/frechet.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")
source("winter2016/marg_dep_mods_funcs.R")

library(dplyr)
library(ggplot2)

lag_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) |>
    na.omit()

  return(df)
}

egpd_sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)
frechet_sim_model <- make_copula_markov_model(margin_frechet, copula_gumbel, stan_mod = NULL)

true_egpd <- c(mu = 0, kappa = 6, sigma = 1, xi = 0.1)
true_frechet <- c(scale = 9, shape = 2.5)


# 1. Setup True Parameters
mod <- "Frechet"

if (mod == "EGPD") {
  sim_model <- egpd_sim_model
  true_margin <- true_egpd
  true_tail_index <- true_margin["xi"]
} else if (mod == "Frechet") {
  sim_model <- frechet_sim_model
  true_margin <- true_frechet
  true_tail_index <- round(1 / true_margin["shape"], 3)
}

true_theta <- 3.0
n_sim <- 10000

data <- sim_model$simulate(
  n = n_sim,
  copula_param = true_theta,
  margin_param = true_margin,
  seed = 123+63
)$x

hist(data, breaks = 50, main = "Simulated Data", xlab = "Value")

# threshrange.plot(data, nint = 50)
# quantile(data, seq(0.9, 0.99, by = 0.01))

# ------------------------------------------------------------------------------
# POINT 1: PROVING THAT IGNORING DEPENDENCE BIASES XI
# ------------------------------------------------------------------------------

# Fit A: Independent (IID) EGPD Model (Ignoring Copula)
fit_iid_egpd <- egpd::fitegpd(data, type = 1)
xi_iid <- fit_iid_egpd$estimate["xi"]

# Fit B: Independent (IID) GPD Model (Standard Threshold Approach)
u_fix <- quantile(data, 0.98)
fit_iid_gpd <- extRemes::fevd(data, threshold = u_fix, type = "GP")
xi_gpd_iid <- fit_iid_gpd$results$par["shape"]

# Fit C: Your Proposition (Joint EGPD + Gumbel)
data_lag <- lag_df(data)

tau <- cor(data_lag$x, data_lag$x_lag, method = "kendall")

theta_tau <- 1 / (1 - tau)
theta_tau

# Initial values from marginal-only fit
init_marg <- egpd::fitegpd(data, method = "mle", model = 1)

# Initial theta_vec
# [log_kappa, log_sigma, xi, log(theta_copula - 1)]
start_params <- c(
  log(init_marg$estimate["kappa"]),
  log(init_marg$estimate["sigma"]),
  init_marg$estimate["xi"],
  log(theta_tau - 1)
)

fit_res_egpd_gumbel <- fit_egpd_gumbel(data, init_par = start_params)

res_kappa <- exp(fit_res_egpd_gumbel$par[1])
res_sigma <- exp(fit_res_egpd_gumbel$par[2])
res_xi <- fit_res_egpd_gumbel$par[3]
res_theta <- exp(fit_res_egpd_gumbel$par[4]) + 1

cat("--- XI Bias Analysis ---\n")
cat("True Tail Index: ", true_tail_index, "\n")
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
  u_val <- quantile(data, u_sequence[i])

  # 1. Censored Gumbel-GPD (Winter2016 style)
  fit_cens <- FitGpd(dat = lag_df(data), u = u_val, optim.type = 2)
  stab_results$xi_censored[i] <- fit_cens$par[3]

  # 2. Threshold-Free EGPD-Gumbel
  stab_results$xi_egpd[i] <- res_xi
}

stab_results |>
  tidyr::pivot_longer(cols = c(xi_censored, xi_egpd), names_to = "model", values_to = "xi") |>
  ggplot(aes(x = prob, y = xi, color = model)) +
  geom_line() +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = "Threshold Stability Comparison",
    subtitle = paste("True tail index =", true_tail_index, "| N =", n_sim),
    x = "Threshold Quantile (u)",
    y = "Estimated xi",
    color = "Model"
  ) +
  theme_minimal()


# ==============================================================================
# Monte Carlo Simulation: Sensitivity to Dependence Strength (Theta)
# ==============================================================================

theta_sequence <- c(1.5, 3.0, 4.5, 6.0)
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
      margin_param = true_margin
    )$x

    # 2. Fit IID Model (The "Naïve" approach)
    fit_egpd_iid <- egpd::fitegpd(sim_data, type = 1)
    xi_egpd_iid <- fit_egpd_iid$estimate["xi"]

    u_fix <- quantile(sim_data, 0.90)
    fit_gpd_iid <- extRemes::fevd(sim_data, threshold = u_fix, type = "GP")
    xi_gpd_iid <- fit_gpd_iid$results$par["shape"]

    # 3. Fit Gumbel dependence +  EGPD
    data_lag <- lag_df(sim_data)

    # Initialise theta
    tau <- cor(data_lag$x, data_lag$x_lag, method = "kendall")

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
    fit_gpd_gumbel_cens <- FitGpd(dat = data_lag, u = u_fix, optim.type = 2)
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

# save(mc_results, file = "sims/copula_markov/egpd_gumbel/res/frechet_dependence_alpha3.RData")
load("sims/copula_markov/egpd_gumbel/res/egpd_mc_dependence_xi01.RData")

# ==============================================================================
# VISUALIZATION OF THE "DEPENDENCE DEGRADATION"
# ==============================================================================
library(ggplot2)

ggplot(mc_results, aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True tail index =", true_tail_index, "| N =", n_sim, "| MC Iterations =", mc_it),
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()

mc_results |>
  filter(model != "Censored_GPD_GUMBEL") |>
  ggplot(aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True tail index =", true_tail_index, "| N =", n_sim, "| MC Iterations =", mc_it),
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()

mc_results |>
  filter(model == "Censored_GPD_GUMBEL" | model == "EGPD_GUMBEL") |>
  ggplot(aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True tail index =", true_tail_index, "| N =", n_sim, "| MC Iterations =", mc_it),
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
    mean_xi = mean(xi_hat),
    sd_xi = sd(xi_hat),
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
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Dependence Strength on Tail Index (xi) Estimation",
    subtitle = paste("True tail index =", true_tail_index, "| N =", n_sim, "| MC Iterations =", mc_it, "| Error bars: MAD"),
    x = "Dependence Strength (Theta)",
    y = "Estimated xi (Median)",
    color = "Model"
  ) +
  theme_minimal()


# ==============================================================================
# SIMPLE DIAGNOSTIC SIMULATION: EGPD + GUMBEL ONLY
# ==============================================================================

set.seed(123)

n_sim <- 2000
mc_it <- 100
theta_true <- 3.0 # fix dependence
results <- data.frame()

for (i in 1:mc_it) {
  # --- 1. Simulate data
  sim_data <- sim_model$simulate(
    n = n_sim,
    copula_param = theta_true,
    margin_param = true_margin
  )$x

  data_lag <- lag_df(sim_data)

  # --- 2. Initialisation
  fit_init <- egpd::fitegpd(sim_data, type = 1)

  tau <- cor(data_lag$x, data_lag$x_lag, method = "kendall")
  theta_tau <- 1 / (1 - tau)

  start_params <- c(
    log(fit_init$estimate["kappa"]),
    log(fit_init$estimate["sigma"]),
    fit_init$estimate["xi"],
    log(theta_tau - 1)
  )

  # --- 3. Fit model (with error handling)
  fit <- tryCatch(
    fit_egpd_gumbel(sim_data, init_par = start_params),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    results <- rbind(results, data.frame(
      iter = i,
      kappa = NA,
      sigma = NA,
      xi = NA,
      theta = NA,
      converged = 0
    ))
    next
  }

  # --- 4. Extract parameters
  kappa_hat <- exp(fit$par[1])
  sigma_hat <- exp(fit$par[2])
  xi_hat <- fit$par[3]
  theta_hat <- exp(fit$par[4]) + 1

  results <- rbind(results, data.frame(
    iter = i,
    kappa = kappa_hat,
    sigma = sigma_hat,
    xi = xi_hat,
    theta = theta_hat,
    converged = 1
  ))
}

library(ggplot2)

results |>
  tidyr::pivot_longer(cols = c(kappa, sigma, xi, theta)) |>
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free") +
  theme_minimal()

true_vals <- data.frame(
  name = c("kappa", "sigma", "xi", "theta"),
  truth = c(true_margin["kappa"], true_margin["sigma"], true_tail_index, theta_true)
)

results |>
  tidyr::pivot_longer(cols = c(kappa, sigma, xi, theta)) |>
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  geom_hline(data = true_vals, aes(yintercept = truth), color = "red", linetype = "dashed") +
  theme_minimal()


summary(results)

# Extreme values
results |>
  filter(abs(xi) > 1 | theta > 20)

pairs(results[, c("kappa", "sigma", "xi", "theta")])

