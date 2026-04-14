# ==============================================================================
# SIMULATION STUDY: THE IMPACT OF DEPENDENCE AND THRESHOLD SELECTION
# ==============================================================================

source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copulas/joe.r")
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

egpd_gumbel_sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)
egpd_joe_sim_model <- make_copula_markov_model(margin_egpd, copula_joe, stan_mod = NULL)

sim_model <- egpd_gumbel_sim_model

# ==============================================================================
# Monte Carlo Simulation: Sensitivity to Dependence Strength (Theta)
# ==============================================================================
library(parallel)
library(doParallel)

set.seed(123)

# Setup
theta_sequence <- c(1.5, 3.0, 4.5, 6.0)
n_sim <- 3500
mc_it <- 50
r_dc <- 1
threshold_probs <- c(0.9, 0.95, 0.96, 0.98)
true_margin_sim <- c(mu = 0, kappa = 6, sigma = 1, xi = 0.1)
true_tail_index <- true_margin_sim["xi"]

# Parallel setup
n_cores <- min(detectCores() - 1, length(theta_sequence))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary objects to cluster
clusterExport(cl, c(
  "n_sim", "mc_it", "r_dc", "threshold_probs", "true_margin_sim", "simulate_copula_markov",
  "sim_model", "lag_df", "fit_egpd_gumbel", "FitGpd", "sample_bisection",
  "margin_egpd", "copula_gumbel", "copula_joe"
))

# Parallel execution
mc_results <- foreach(
  th = theta_sequence,
  .combine = rbind,
  .packages = c("extRemes", "egpd")
) %dopar% {
  results <- vector("list", mc_it)

  cat(sprintf(
    "[%s] Starting theta = %.1f (%d iterations)\n",
    Sys.time(), th, mc_it
  ))

  for (i in 1:mc_it) {
    # Simulate
    sim_data <- sim_model$simulate(
      n = n_sim,
      copula_param = th,
      margin_param = true_margin_sim
    )

    # Store data in lag_df form
    data_lag <- lag_df(sim_data$x)

    # Fit threshold indep models
    fit_egpd_iid <- egpd::fitegpd(sim_data$x, type = 1)

    fit_egpd_dep <- fit_egpd_gumbel(sim_data$x, method = "Nelder-Mead")


    # Compute thresholds
    thresholds <- quantile(sim_data$x, probs = threshold_probs)

    thresh_results <- list()

    for (p in seq_along(threshold_probs)) {
      u_fix <- thresholds[p]

      fit_gpd_iid <- extRemes::fevd(sim_data$x, threshold = u_fix, type = "GP")

      sim_dcruns <- extRemes::decluster(sim_data$x, threshold = u_fix, r = r_dc)

      fit_gpd_dcruns <- extRemes::fevd(sim_dcruns, threshold = u_fix, type = "GP")

      fit_gpd_gumbel_cens <- FitGpd(dat = data_lag, u = u_fix, optim.type = 2)

      # --- store results for THIS threshold ---
      thresh_results[[p]] <- data.frame(
        theta = th,
        iteration = i,
        threshold = threshold_probs[p],
        model = c(
          "IID_EGPD", "IID_GPD", "GPD_DCRUNS",
          "Censored_GPD_GUMBEL", "EGPD_GUMBEL"
        ),
        xi_hat = c(
          fit_egpd_iid$estimate["xi"], # same across thresholds
          fit_gpd_iid$results$par["shape"],
          fit_gpd_dcruns$results$par["shape"],
          fit_gpd_gumbel_cens$par[3],
          fit_egpd_dep$estimate["xi"] # same across thresholds
        ),
        stringsAsFactors = FALSE
      )
    }

    # Combine all thresholds for this iteration
    results[[i]] <- do.call(rbind, thresh_results)
  }

  do.call(rbind, results)
}

stopCluster(cl)

save(mc_results, file = "sims/copula_markov/egpd_gumbel/res/egpd_xi01_complete_multiplethresh.RData")
load("sims/copula_markov/egpd_gumbel/res/egpd_xi01_complete_multiplethresh.RData")


# ==============================================================================
# VISUALIZATION OF THE "DEPENDENCE DEGRADATION"
# ==============================================================================
plot_title <- "Effect of Dependence Strength on Tail Index (xi) Estimation - Gumbel EGPD Simulation"
plot_subtitle <- paste(
  "True tail index =", true_tail_index,
  "| N =", n_sim,
  "| MC Iterations =", mc_it,
  "| Declustering run length =", r_dc
)

mc_results$model_type <- ifelse(
  mc_results$model %in% c("IID_GPD", "GPD_DCRUNS", "Censored_GPD_GUMBEL"),
  "Threshold-dependent",
  "Threshold-independent"
)

ggplot(mc_results, aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  facet_wrap(~threshold) +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()

mc_results |>
  filter(model != "Censored_GPD_GUMBEL") |>
  ggplot(aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  facet_wrap(~threshold) +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()

mc_results |>
  filter(model == "Censored_GPD_GUMBEL" | model == "EGPD_GUMBEL") |>
  ggplot(aes(x = as.factor(theta), y = xi_hat, fill = model)) +
  geom_boxplot() +
  facet_wrap(~threshold) +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "Dependence Strength (Theta)",
    y = "Estimated xi"
  ) +
  theme_minimal()


# Calculate median and MAD for each theta-model combination

summary_stats <- mc_results |>
  group_by(theta, model, threshold) |>
  summarise(
    median_xi = median(xi_hat),
    mad_xi = mad(xi_hat),
    mean_xi = mean(xi_hat),
    sd_xi = sd(xi_hat),
    .groups = "drop"
  )

ggplot(summary_stats, aes(x = as.factor(theta), y = mean_xi, color = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = mean_xi - sd_xi, ymax = mean_xi + sd_xi),
    width = 0.2,
    linewidth = 0.6,
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~threshold) +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "Dependence Strength (Theta)",
    y = "Estimated xi (Mean)",
    color = "Model"
  ) +
  theme_minimal()

ggplot(summary_stats, aes(x = as.factor(theta), y = median_xi, color = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = median_xi - mad_xi, ymax = median_xi + mad_xi),
    width = 0.2,
    linewidth = 0.6,
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~threshold) +
  geom_hline(yintercept = true_tail_index, linetype = "dashed", color = "red") +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "Dependence Strength (Theta)",
    y = "Estimated xi (Median)",
    color = "Model"
  ) +
  theme_minimal()


### Bigger simulation study for consistency

# --- Generator Library ---

# 1. Original Markov Gumbel (Your "Home Field" model)
gen_gumbel_markov <- function(n, dep_param, xi, kappa=6, sigma=1) {
  sim_data <- sim_model$simulate(
    n = n,
    copula_param = dep_param,
    margin_param = c(mu=0, kappa=kappa, sigma=sigma, xi=xi)
  )
  return(sim_data$x)
}

# 2. Raw AR(1) Pareto (Linear Misspecification)
gen_ar1_pareto <- function(n, dep_param, xi) {
  burn_in <- 500
  total_n <- n + burn_in
  innov <- (1 - runif(total_n))^(-xi)
  x <- numeric(total_n)
  x[1] <- innov[1]
  for(t in 2:total_n) x[t] <- dep_param * x[t-1] + innov[t]
  return(x[(burn_in + 1):total_n])
}

run_study <- function(generator_fn, 
                      dep_sequence, 
                      n_sequence, 
                      mc_it = 100, 
                      true_xi = 0.1,
                      threshold_probs = c(0.9, 0.95),
                      r_dc = 1,
                      extra_gen_args = list()) {
  
  n_cores <- min(parallel::detectCores() - 1, length(dep_sequence))
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
  
  results <- foreach(
    dep_val = dep_sequence,
    .combine = rbind,
    .packages = c("extRemes", "egpd", "tidyverse")
  ) %dopar% {
    
    # Helper to safely extract xi or return NA
    get_xi <- function(fit, type) {
      if (is.null(fit)) return(NA)
      if (type == "egpd") return(fit$estimate["xi"])
      if (type == "fevd") return(fit$results$par["shape"])
      if (type == "censored") return(fit$par[3])
      return(NA)
    }
    
    iter_results <- list()
    
    for (n_sim in n_sequence) {
      for (i in 1:mc_it) {
        gen_args <- c(list(n = n_sim, dep_param = dep_val, xi = true_xi), extra_gen_args)
        x_sim <- do.call(generator_fn, gen_args)
        
        data_lag <- lag_df(x_sim)
        thresholds <- quantile(x_sim, probs = threshold_probs)
        
        fit_egpd_iid <- tryCatch(egpd::fitegpd(x_sim, type = 1), error = function(e) NULL)
        fit_egpd_dep <- tryCatch(fit_egpd_gumbel(x_sim, method = "Nelder-Mead"), error = function(e) NULL)
        
        for (p in seq_along(threshold_probs)) {
          u_fix <- thresholds[p]
          
          fit_gpd_iid <- tryCatch(extRemes::fevd(x_sim, threshold = u_fix, type = "GP"), error = function(e) NULL)
          sim_dcruns <- tryCatch(extRemes::decluster(x_sim, threshold = u_fix, r = r_dc), error = function(e) NULL)
          fit_gpd_dcruns <- tryCatch(extRemes::fevd(sim_dcruns, threshold = u_fix, type = "GP"), error = function(e) NULL)
          fit_gpd_gumbel_cens <- tryCatch(FitGpd(dat = data_lag, u = u_fix, optim.type = 2), error = function(e) NULL)
          
          # Now xi_hat is guaranteed to have length 5
          xi_hat_vec <- c(
            get_xi(fit_egpd_iid, "egpd"),
            get_xi(fit_gpd_iid, "fevd"),
            get_xi(fit_gpd_dcruns, "fevd"),
            get_xi(fit_gpd_gumbel_cens, "censored"),
            get_xi(fit_egpd_dep, "egpd")
          )
          
          iter_results[[length(iter_results) + 1]] <- data.frame(
            dep_val = dep_val, 
            n = n_sim, 
            iteration = i, 
            threshold = threshold_probs[p],
            model = c("IID_EGPD", "IID_GPD", "GPD_DCRUNS", "Censored_GPD_GUMBEL", "EGPD_GUMBEL"),
            xi_hat = xi_hat_vec
          )
        }
      }
    }
    do.call(rbind, iter_results)
  }
  
  parallel::stopCluster(cl)
  return(results %>% mutate(bias = xi_hat - true_xi))
}

set.seed(123)

results_gumbel <- run_study(
  generator_fn = gen_gumbel_markov,
  dep_sequence = c(1.5, 3, 4.5, 6),
  n_sequence = c(500, 1000, 2000, 4000, 8000),
  mc_it = 200,
  extra_gen_args = list(kappa = 6)
)

save(results_gumbel, file = "sims/copula_markov/egpd_gumbel/res/egpd_gumbel_200mc_consistency_results_xi01.RData")


# Compute summary statistics (mean bias, sd, RMSE) for each combination
summary_results <- results_gumbel %>%
  group_by(dep_val, n, model, threshold) %>% # Use dep_val as per your run_study code
  summarise(
    # Calculations skip NAs
    mean_bias   = mean(bias, na.rm = TRUE),
    sd_bias     = sd(bias, na.rm = TRUE),
    median_bias = median(bias, na.rm = TRUE),
    mad_bias    = mad(bias, na.rm = TRUE),
    rmse        = sqrt(mean(bias^2, na.rm = TRUE)),
    
    # Track counts
    total_attempts = n(),                       # How many iterations we tried
    converged_obs  = sum(!is.na(xi_hat)),       # How many actually worked
    conv_rate      = sum(!is.na(xi_hat)) / n(), # Percentage success
    
    .groups = "drop"
  )

# Create faceted plot by theta and threshold
# You can choose to focus on one threshold or show all
selected_threshold <- 0.90  # Choose your preferred threshold

plot_data <- summary_results %>%
  filter(threshold == selected_threshold)

p <- ggplot(plot_data, aes(x = n, y = median_bias, color = model, linetype = model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ dep_val, nrow = 2, 
             labeller = labeller(theta = function(x) paste("θ =", x))) +
  scale_x_continuous(trans = "log10", breaks = n_sequence) +
  labs(
    title = sprintf("Bias in Tail Index Estimation (threshold = %.2f)", selected_threshold),
    subtitle = sprintf("True ξ = %.2f", true_tail_index),
    x = "Sample Size (n)",
    y = "Bias (ξ̂ - ξ)",
    color = "Method",
    linetype = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)


# If you want to verify n^(-α) convergence rate
p_loglog <- ggplot(plot_data, aes(x = n, y = abs(median_bias), color = model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  # Add reference lines for different convergence rates
  geom_abline(intercept = log10(0.5), slope = -0.5, 
              linetype = "dotted", color = "black", alpha = 0.5) +  # n^(-1/2)
  geom_abline(intercept = log10(1), slope = -1, 
              linetype = "dotted", color = "black", alpha = 0.5) +    # n^(-1)
  annotate("text", x = 1000, y = 0.5/sqrt(1000), label = "n^(-1/2)", 
           size = 3, hjust = -0.2) +
  annotate("text", x = 1000, y = 1/1000, label = "n^(-1)", 
           size = 3, hjust = -0.2) +
  facet_wrap(~ dep_val, nrow = 2) +
  scale_x_log10(breaks = n_sequence, labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = "Convergence Rate Analysis (log-log scale)",
    x = "Sample Size (n)",
    y = "|Bias|",
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_loglog)


