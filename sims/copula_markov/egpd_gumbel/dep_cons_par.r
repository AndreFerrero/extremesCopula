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
  if (is.null(fit) || inherits(fit, "error")) {
    return(NA)
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
    "EGPD_GUMBEL_IFM" = list(
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
# MAIN FUNCTION
# ==============================================================================

run_consistency_study_parallel <- function(
  dep_sequence,
  n_sequence,
  mc_iterations = 200,
  true_params = list(kappa = 6, sigma = 1, xi = 0.1),
  threshold_probs = c(0.90, 0.95),
  seed = 123,
  system_cores = NULL
) {
  dir.create("logs", showWarnings = FALSE)

  param_grid <- expand.grid(
    theta_true = dep_sequence,
    n_sim = n_sequence,
    KEEP.OUT.ATTRS = FALSE
  )

  if (is.null(system_cores)) {
    system_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  }

  n_jobs <- nrow(param_grid)
  n_cores <- min(system_cores, n_jobs)

  cat(sprintf("Using %d cores for %d jobs\n", n_cores, n_jobs))
  cat(sprintf("MC iterations per job: %d\n", mc_iterations))

  cl <- makeCluster(n_cores)

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

  # --------------------------------------------------------------------------
  # LOG FILES
  # --------------------------------------------------------------------------
  progress_file <- "logs/progress_global.log"
  audit_file <- "logs/audit_state.log"
  state_dir <- "logs/state"

  dir.create(state_dir, showWarnings = FALSE)

  progress_step <- 10 # THROTTLE

  clusterExport(cl, varlist = c(
    "lag_df",
    "extract_model_params",
    "true_params",
    "threshold_probs",
    "mc_iterations",
    "seed",
    "param_grid",
    "progress_file",
    "audit_file",
    "state_dir",
    "progress_step"
  ), envir = environment())

  # --------------------------------------------------------------------------
  # PARALLEL EXECUTION
  # --------------------------------------------------------------------------
  results_list <- parLapplyLB(cl, seq_len(n_jobs), function(job_idx) {
    theta_true <- param_grid$theta_true[job_idx]
    n_sim <- param_grid$n_sim[job_idx]

    job_id <- sprintf("theta_%g_n_%d", theta_true, n_sim)
    state_file <- file.path(state_dir, paste0(job_id, ".rds"))


    sim_model <- make_copula_markov_model(
      margin_egpd,
      copula_gumbel,
      stan_mod = NULL
    )

    # ----------------------------------------------------------------------
    # LOAD STATE (CRITICAL FIX)
    # ----------------------------------------------------------------------
    if (file.exists(state_file)) {
      state <- readRDS(state_file)
      completed_iters <- state$completed
    } else {
      completed_iters <- integer(0)
    }

    job_results <- vector("list", mc_iterations)

    for (iter in 1:mc_iterations) {
      set.seed(seed + job_idx * 10000 + iter)

      iter_result <- NULL
      error_msg <- NULL

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

          params_iid <- extract_model_params(
            tryCatch(egpd::fitegpd(x_sim, type = 1), error = function(e) e),
            "IID_EGPD"
          )

          params_gumbel <- extract_model_params(
            tryCatch(fit_egpd_gumbel_copula(x_sim), error = function(e) e),
            "EGPD_GUMBEL"
          )

          params_gumbel_ifm <- extract_model_params(
            tryCatch(fit_egpd_gumbel_copula(x_sim, method = "ifm"), error = function(e) e),
            "EGPD_GUMBEL_IFM"
          )
          params_joe <- extract_model_params(
            tryCatch(fit_egpd_joe_copula(x_sim), error = function(e) e),
            "EGPD_JOE"
          )


          k_val <- floor(n_sim^0.7)

          hill_fit <- tryCatch(
              hill_bc_hat(x_sim, k_val),
              error = function(e) e
            )

          hill <- extract_model_params(
            hill_fit["Hill"],
            "Hill"
          )

          hill_bc <- extract_model_params(
            hill_fit["Hill_BC"],
            "Hill_BC"
          )

          iter_rows <- list()

          # ----------------------------------------------------------
          # 1. NON-THRESHOLD MODELS (run ONCE)
          # ----------------------------------------------------------
          base_res <- dplyr::bind_rows(
            params_iid,
            params_gumbel,
            params_gumbel_ifm,
            params_joe,
            hill,
            hill_bc
          )

          iter_rows[[1]] <- base_res

          for (p_idx in seq_along(threshold_probs)) {
            # ------------------------------------------------------------------
            # PRECOMPUTE ALL INPUTS
            # ------------------------------------------------------------------
            u_fix <- threshold_vals[p_idx]
            thresh_prob <- threshold_probs[p_idx]

            # ------------------------------------------------------------------
            # FIT MODELS
            # ------------------------------------------------------------------
            fit_iid_gpd <- tryCatch(
              extRemes::fevd(x_sim, threshold = u_fix, type = "GP"),
              error = function(e) e
            )

            fit_declust <- tryCatch(
              extRemes::fevd(
                extRemes::decluster(x_sim, threshold = u_fix),
                threshold = u_fix,
                type = "GP"
              ),
              error = function(e) e
            )

            fit_censored <- tryCatch(
              FitGpd(data_lag, u_fix),
              error = function(e) e
            )

            # ------------------------------------------------------------------
            # EXTRACT PARAMETERS
            # ------------------------------------------------------------------

            res_thresh <- dplyr::bind_rows(
              extract_model_params(fit_iid_gpd, "IID_GPD", thresh_prob),
              extract_model_params(fit_declust, "GPD_Declustering", thresh_prob),
              extract_model_params(fit_censored, "Censored_GPD_GUMBEL", thresh_prob)
            )

            # ------------------------------------------------------------------
            # STORE IF VALID
            # ------------------------------------------------------------------
            if (!is.null(res_thresh) && nrow(res_thresh) > 0) {
              iter_rows[[length(iter_rows) + 1]] <- res_thresh
            }
          }

          if (length(iter_rows) > 0) {
            iter_result <- bind_rows(iter_rows)
          }
        },
        error = function(e) {
          error_msg <<- e$message
        }
      )

      if (is.null(iter_result)) {
        iter_result <- data.frame(
          model = NA,
          parameter = NA,
          estimate = NA,
          threshold = NA
        )
      }

      iter_result$theta_true <- theta_true
      iter_result$n <- n_sim
      iter_result$iteration <- iter
      iter_result$error <- error_msg

      job_results[[iter]] <- iter_result

      # ----------------------------------------------------------------------
      # UPDATE STATE FILE (CRITICAL FIX)
      # ----------------------------------------------------------------------
      completed_iters <- unique(c(completed_iters, iter))

      saveRDS(
        list(
          completed = completed_iters,
          last_update = Sys.time()
        ),
        state_file
      )

      # ----------------------------------------------------------------------
      # THROTTLED LOGGING
      # ----------------------------------------------------------------------
      if (iter %% progress_step == 0 || iter == mc_iterations) {
        cat(
          sprintf(
            "theta=%.2f n=%d done=%d/%d missing=%d\n",
            theta_true,
            n_sim,
            length(completed_iters),
            mc_iterations,
            mc_iterations - length(completed_iters)
          ),
          file = audit_file,
          append = TRUE
        )
      }

      cat(
        sprintf(
          "DONE %.2f %d %d\n",
          theta_true, n_sim, iter
        ),
        file = progress_file,
        append = TRUE
      )
    }

    bind_rows(job_results)
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

now <- Sys.time()
cat("Launching at ", now, "\n")

results <- run_consistency_study_parallel(
  dep_sequence = c(1, 2, 4, 6),
  n_sequence = c(250, 500, 1000, 2000, 4000, 8000),
  mc_iterations = 300,
  true_params = list(kappa = 2, sigma = 1, xi = 0.1),
  threshold_probs = c(0.90, 0.95),
  seed = 123
)

after <- Sys.time()
cat("Simulation completed at ", after, "\n")

diff <- as.numeric(difftime(after, now, units = "hours"))
cat("Elapsed time", diff, "\n")

save(results, file = here(script_dir, "res/consistency_mc300_u9095_neldermead_gumbeldata_ifm_joe_hill_kn07_kappa2_sigma1_xi01.RData"))

print("Results saved")
