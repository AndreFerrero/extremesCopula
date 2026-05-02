library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

script_dir <- here("sims", "copula_markov", "egpd_gumbel")

# ==============================================================================
# HELPER
# ==============================================================================

lag_df <- function(x) {
  data.frame(x = x, x_lag = dplyr::lag(x)) |> na.omit()
}

extract_model_params <- function(fit, model_type, threshold = NA) {
  if (is.null(fit)) {
    return(NULL)
  }

  params <- switch(model_type,
    "IID_EGPD" = list(
      kappa = fit$estimate["kappa"],
      sigma = fit$estimate["sigma"],
      xi = fit$estimate["xi"]
    ),
    "EGPD_GUMBEL" = list(
      kappa = fit$estimate["kappa"],
      sigma = fit$estimate["sigma"],
      xi = fit$estimate["xi"],
      theta = fit$estimate["theta"]
    ),
    "EGPD_JOE" = list(
      kappa = fit$estimate["kappa"],
      sigma = fit$estimate["sigma"],
      xi = fit$estimate["xi"],
      theta = fit$estimate["theta"]
    ),
    "IID_GPD" = list(
      sigma_u = fit$results$par["scale"],
      xi = fit$results$par["shape"]
    ),
    "GPD_Declustering" = list(
      sigma_u = fit$results$par["scale"],
      xi = fit$results$par["shape"]
    ),
    "Censored_GPD_GUMBEL" = list(
      theta = 1 / fit$par[1],
      sigma_u = fit$par[2],
      xi = fit$par[3],
      lambda_u = fit$par[4]
    ),
    "Hill" = list(xi = fit),
    "Hill_BC" = list(xi = fit)
  )

  data.frame(
    model = model_type,
    parameter = names(params),
    estimate = as.numeric(unlist(params)),
    threshold = threshold,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# PARALLEL SIMULATION
# ==============================================================================

run_consistency_study_parallel <- function(
  dep_sequence,
  n_sequence,
  mc_iterations = 5,
  true_params = list(kappa = 6, sigma = 1, xi = 0.1),
  threshold_probs = c(0.90, 0.95),
  seed = 123,
  system_cores = NULL
) {
  param_grid <- expand.grid(
    theta_true = dep_sequence,
    n_sim = n_sequence,
    KEEP.OUT.ATTRS = FALSE
  )

  if (is.null(system_cores)) {
    system_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  }

  n_jobs <- nrow(param_grid)
  print(sprintf("Requested %d cores, %d jobs to run.", system_cores, n_jobs))

  n_cores <- min(system_cores, n_jobs)

  cat(sprintf("Using %d cores for %d jobs\n", n_cores, n_jobs))

  cat(sprintf("MC iterations per job: %d\n", mc_iterations))

  print("starting parallel execution...")

  cl <- makeCluster(n_cores)

  # --- Load EVERYTHING needed on workers ---
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(extRemes)
    library(egpd)
    library(here)
    code_mod_dir <- here("code", "models")

    source(here(code_mod_dir, "copula_markov/copula_markov_model.r"))
    source(here(code_mod_dir, "margins/egpd.r"))
    source(here(code_mod_dir, "copulas/gumbel.r"))
    source(here(code_mod_dir, "copula_markov/egpd_joe.r"))
    source(here(code_mod_dir, "copula_markov/egpd_gumbel.r"))
    source(here(code_mod_dir, "gpd_gumbel.r"))
    source(here(code_mod_dir, "hill_estimator.r"))
  })

  # --- Export only lightweight objects ---
  clusterExport(cl, varlist = c(
    "lag_df",
    "extract_model_params",
    "true_params",
    "threshold_probs",
    "mc_iterations",
    "seed",
    "param_grid"
  ), envir = environment())

  results_list <- parLapplyLB(cl, seq_len(n_jobs), function(job_idx) {
    theta_true <- param_grid$theta_true[job_idx]
    n_sim <- param_grid$n_sim[job_idx]

    set.seed(seed + job_idx)

    sim_model <- make_copula_markov_model(
      margin_egpd,
      copula_gumbel,
      stan_mod = NULL
    )

    job_results <- list()

    for (iter in 1:mc_iterations) {
      tryCatch(
        {
          sim_data <- sim_model$simulate(
            n = n_sim,
            copula_param = theta_true,
            margin_param = c(
              mu = 0,
              kappa = true_params$kappa,
              sigma = true_params$sigma,
              xi = true_params$xi
            )
          )

          x_sim <- sim_data$x
          data_lag <- lag_df(x_sim)
          threshold_vals <- quantile(x_sim, probs = threshold_probs)

          # --- threshold-independent ---
          fit_egpd_iid <- tryCatch(
            egpd::fitegpd(x_sim, type = 1),
            error = function(e) NULL
          )

          fit_egpd_gumbel <- tryCatch(
            fit_egpd_gumbel_copula(x_sim, optim.method = "L-BFGS-B"),
            error = function(e) NULL
          )

          fit_egpd_joe <- tryCatch(
            fit_egpd_joe_copula(x_sim, optim.method = "L-BFGS-B"),
            error = function(e) NULL
          )

          params_iid <- extract_model_params(
            fit_egpd_iid, "IID_EGPD",
            threshold = NA
          )

          params_gumbel <- extract_model_params(
            fit_egpd_gumbel, "EGPD_GUMBEL",
            threshold = NA
          )

          params_joe <- extract_model_params(
            fit_egpd_joe, "EGPD_JOE",
            threshold = NA
          )

          # --- threshold-dependent ---
          for (p_idx in seq_along(threshold_probs)) {
            u_fix <- threshold_vals[p_idx]
            thresh_prob <- threshold_probs[p_idx]
            k_val <- floor(n_sim * (1 - thresh_prob))

            fit_gpd_iid <- tryCatch(
              extRemes::fevd(x_sim, threshold = u_fix, type = "GP"),
              error = function(e) NULL
            )

            sim_dc <- tryCatch(
              extRemes::decluster(
                x_sim,
                threshold = u_fix,
                method = "intervals"
              ),
              error = function(e) NULL
            )

            fit_gpd_dc <- if (!is.null(sim_dc)) {
              tryCatch(
                extRemes::fevd(sim_dc, threshold = u_fix, type = "GP"),
                error = function(e) NULL
              )
            } else {
              NULL
            }

            fit_gpd_cens <- tryCatch(
              FitGpd(dat = data_lag, u = u_fix, optim.type = 2),
              error = function(e) NULL
            )

            hill_estimates <- tryCatch(
              hill_bc_hat(x_sim, k_target = k_val),
              error = function(e) NULL
            )

            threshold_res <- dplyr::bind_rows(
              params_iid,
              params_gumbel,
              params_joe,
              extract_model_params(
                fit_gpd_iid,
                "IID_GPD",
                threshold = thresh_prob
              ),
              extract_model_params(
                fit_gpd_dc,
                "GPD_Declustering",
                threshold = thresh_prob
              ),
              extract_model_params(
                fit_gpd_cens,
                "Censored_GPD_GUMBEL",
                threshold = thresh_prob
              ),
              extract_model_params(
                hill_estimates$hill,
                "Hill",
                threshold = thresh_prob
              ),
              extract_model_params(
                hill_estimates$hill_bc,
                "Hill_BC",
                threshold = thresh_prob
              )
            )

            if (nrow(threshold_res) > 0) {
              threshold_res$theta_true <- theta_true
              threshold_res$n <- n_sim
              threshold_res$iteration <- iter
              job_results[[length(job_results) + 1]] <- threshold_res
            }
          }
        },
        error = function(e) {
          message(sprintf(
            "Error (theta=%.2f, n=%d, iter=%d): %s",
            theta_true, n_sim, iter, e$message
          ))
        }
      )
    }

    dplyr::bind_rows(job_results)
  })

  stopCluster(cl)

  # --- Combine results ---
  results_df <- dplyr::bind_rows(results_list)

  if (nrow(results_df) == 0) {
    print(results_df)
    stop("No results returned from workers.")
  }

  # --- Post-processing ---
  results_df <- results_df %>%
    mutate(
      true_value = case_when(
        parameter == "kappa" ~ true_params$kappa,
        parameter == "sigma" ~ true_params$sigma,
        parameter == "xi" ~ true_params$xi,
        parameter == "theta" ~ theta_true,
        TRUE ~ NA_real_
      ),
      bias = estimate - true_value
    )

  return(results_df)
}

results <- run_consistency_study_parallel(
  dep_sequence = c(1.5, 3, 4.5, 6),
  n_sequence = c(250, 500, 1000, 2000, 4000, 8000, 10000),
  mc_iterations = 200,
  true_params = list(kappa = 6, sigma = 1, xi = 0.1),
  threshold_probs = c(0.90, 0.95),
  seed = 123
)

save(results, file = here(script_dir, "res/consistency_lbfgsb_gumbeldata_joe_hill_kappa6_sigma1_xi01.RData"))
print("Simulation completed and results saved.")
