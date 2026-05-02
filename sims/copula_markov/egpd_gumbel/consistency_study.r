# ==============================================================================
# CLEAN SIMULATION STUDY: PARAMETER CONSISTENCY ACROSS MODELS
# ==============================================================================

# Load dependencies
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")
source("code/models/gpd_gumbel.r")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

lag_df <- function(x) {
  data.frame(x = x, x_lag = dplyr::lag(x)) |> na.omit()
}

# Extract parameters in standardized format
# ==============================================================================
# UPDATED HELPER FUNCTION
# ==============================================================================

# ==============================================================================
# UPDATED HELPER FUNCTION
# ==============================================================================

# Extract parameters in standardized format with threshold as a metadata column
extract_model_params <- function(fit, model_type, threshold = NA) {
  if (is.null(fit)) return(NULL)
  
  params <- switch(model_type,
    "IID_EGPD" = {
      list(kappa = fit$estimate["kappa"], sigma = fit$estimate["sigma"], xi = fit$estimate["xi"])
    },
    "EGPD_GUMBEL" = {
      list(kappa = fit$estimate["kappa"], sigma = fit$estimate["sigma"], xi = fit$estimate["xi"], theta = fit$estimate["theta"])
    },
    "IID_GPD" = {
      list(sigma_u = fit$results$par["scale"], xi = fit$results$par["shape"])
    },
    "GPD_DCRUNS" = {
      list(sigma_u = fit$results$par["scale"], xi = fit$results$par["shape"])
    },
    "Censored_GPD_GUMBEL" = {
      list(theta = 1/fit$par[1], sigma_u = fit$par[2], xi = fit$par[3], lambda_u = fit$par[4])
    }
  )
  
  data.frame(
    model = model_type,
    parameter = names(params),
    estimate = as.numeric(unlist(params)),
    threshold = threshold, # This is now a dedicated column
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# UPDATED SIMULATION FUNCTION
# ==============================================================================

run_consistency_study <- function(
  dep_sequence,
  n_sequence,
  mc_iterations = 5,
  true_params = list(kappa = 6, sigma = 1, xi = 0.1),
  threshold_probs = c(0.90, 0.95),
  r_dc = 1,
  seed = 123
) {
  
  set.seed(seed)
  sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)
  
  all_results <- list()
  
  for (theta_true in dep_sequence) {
    cat(sprintf("\n=== Processing theta = %.1f ===\n", theta_true))
    
    for (n_sim in n_sequence) {
      cat(sprintf("  n = %d\n", n_sim))
      
      for (iter in 1:mc_iterations) {
        cat(sprintf("    Iteration %d/%d...", iter, mc_iterations))
        
        tryCatch({
          # Generate data
          sim_data <- sim_model$simulate(
            n = n_sim,
            copula_param = theta_true,
            margin_param = c(mu = 0, kappa = true_params$kappa, sigma = true_params$sigma, xi = true_params$xi)
          )
          
          x_sim <- sim_data$x
          data_lag <- lag_df(x_sim)
          threshold_vals <- quantile(x_sim, probs = threshold_probs)
          
          # 1. Fit threshold-independent models (Threshold = NA or "Global")
          fit_egpd_iid <- tryCatch(egpd::fitegpd(x_sim, type = 1), error = function(e) NULL)
          fit_egpd_gumbel <- tryCatch(fit_egpd_gumbel_copula(x_sim), error = function(e) NULL)
          
          params_iid <- extract_model_params(fit_egpd_iid, "IID_EGPD", threshold = NA)
          params_gumbel <- extract_model_params(fit_egpd_gumbel, "EGPD_GUMBEL", threshold = NA)
          
          # 2. Fit threshold-dependent models
          for (p_idx in seq_along(threshold_probs)) {
            u_fix <- threshold_vals[p_idx]
            thresh_prob <- threshold_probs[p_idx]
            
            fit_gpd_iid <- tryCatch(extRemes::fevd(x_sim, threshold = u_fix, type = "GP"), error = function(e) NULL)
            
            sim_dcruns <- tryCatch(extRemes::decluster(x_sim, threshold = u_fix, r = r_dc), error = function(e) NULL)
            fit_gpd_dcruns <- if (!is.null(sim_dcruns)) {
              tryCatch(extRemes::fevd(sim_dcruns, threshold = u_fix, type = "GP"), error = function(e) NULL)
            } else NULL
            
            fit_gpd_cens <- tryCatch(FitGpd(dat = data_lag, u = u_fix, optim.type = 2), error = function(e) NULL)
            
            # Combine all for this threshold
            # Note: We include params_iid and params_gumbel in each threshold block 
            # so you can compare them side-by-side in plots/tables.
            threshold_res <- bind_rows(
              params_iid,
              params_gumbel,
              extract_model_params(fit_gpd_iid, "IID_GPD", threshold = thresh_prob),
              extract_model_params(fit_gpd_dcruns, "GPD_DCRUNS", threshold = thresh_prob),
              extract_model_params(fit_gpd_cens, "Censored_GPD_GUMBEL", threshold = thresh_prob)
            )
            
            if (nrow(threshold_res) > 0) {
              threshold_res$theta_true <- theta_true
              threshold_res$n <- n_sim
              threshold_res$iteration <- iter
              all_results[[length(all_results) + 1]] <- threshold_res
            }
          }
          cat(" done\n")
        }, error = function(e) cat(sprintf(" ERROR: %s\n", e$message)))
      }
    }
  }
  
  results_df <- bind_rows(all_results)
  
  # Final Processing
  results_df <- results_df %>%
    mutate(
      true_value = case_when(
        parameter == "kappa" ~ true_params$kappa,
        parameter == "sigma" ~ true_params$sigma,
        parameter == "xi"    ~ true_params$xi,
        parameter == "theta" ~ theta_true,
        TRUE ~ NA_real_
      ),
      bias = estimate - true_value
    )
  
  return(results_df)
}

# ==============================================================================
# RUN
# ==============================================================================

cat("Starting run...\n")

results <- run_consistency_study(
  dep_sequence = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
  n_sequence = c(500, 1000, 2000, 4000, 8000),
  mc_iterations = 100,
  true_params = list(kappa = 6, sigma = 1, xi = 0.1),
  threshold_probs = c(0.90, 0.95),
  seed = 123
)

# Compare xi estimates across different thresholds and models
results %>%
  filter(parameter == "xi") %>%
  ggplot(aes(x = factor(threshold), y = estimate, fill = model)) +
  geom_boxplot() +
  facet_wrap(~n) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed")