==============================================================================
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

# 1. Original Markov Gumbel (Your "Home Field" model)
gen_gumbel_markov <- function(n, dep_param, xi, kappa=6, sigma=1) {
  sim_data <- sim_model$simulate(
    n = n,
    copula_param = dep_param,
    margin_param = c(mu=0, kappa=kappa, sigma=sigma, xi=xi)
  )
  return(sim_data$x)
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