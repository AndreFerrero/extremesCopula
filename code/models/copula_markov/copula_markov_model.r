make_copula_markov_model <- function(
  margin,
  copula,
  stan_file
) {
  # ----------------------------------------
  # Compile (auto-cached by rstan)
  # ----------------------------------------

  message("Compiling (or reusing cached) Stan model...")

  stan_model_obj <- rstan::stan_model(stan_file)

  # ----------------------------------------
  # Build object
  # ----------------------------------------

  list(
    margin = margin,
    copula = copula,
    stan_model = stan_model_obj,
    stan_file = stan_file,

    # -----------------------
    # SIMULATION
    # -----------------------
    simulate = function(n,
                        copula_param,
                        margin_param,
                        seed = NULL) {
      if (!is.null(seed)) set.seed(seed)

      simulate_copula_markov(
        n = n,
        copula = copula,
        copula_param = copula_param,
        margin = margin,
        margin_param = margin_param
      )
    },

    # -----------------------
    # FIT
    # -----------------------
    fit = function(x,
                   iter = 2000,
                   chains = 4,
                   seed = 42,
                   prior_check = 0,
                   run_ppc = 1,
                   I = 25,
                   compute_loo = TRUE,
                   adapt_delta = NULL) {
      options(mc.cores = chains)
      # control list
      control_list <- if (!is.null(adapt_delta)) list(adapt_delta = adapt_delta) else list()

      # Stan data
      stan_data <- list(
        T = length(x),
        x = x,
        prior_check = prior_check,
        run_ppc = run_ppc
      )

      if (!is.null(I)) stan_data$I <- I

      message("Sampling...")

      fit_obj <- rstan::sampling(
        object = stan_model_obj,
        data = stan_data,
        iter = iter,
        chains = chains,
        seed = seed,
        control = control_list
      )

      # Extract all posterior draws
      posterior_list <- rstan::extract(fit_obj, permuted = TRUE)

      # Remove generated quantities we don't want
      exclude <- c("x_rep", "log_lik", "lp__")
      param_names <- setdiff(names(posterior_list), exclude)

      # Convert to dataframe
      draws_df <- as.data.frame(posterior_list[param_names])

      margin_draws <- draws_df[, margin$name_param, drop = FALSE]
      copula_draws <- draws_df[, copula$name_param, drop = FALSE]

      # Build output
      out <- list(
        fit = fit_obj,
        x = x,
        draws = draws_df,
        margin_draws = margin_draws,
        copula_draws = copula_draws
      )

      # PPC
      if (run_ppc == 1) {
        out$ppc <- rstan::extract(fit_obj, "x_rep")$x_rep
      }

      # LOO
      if (compute_loo) {
        log_lik <- rstan::extract(fit_obj, "log_lik")$log_lik
        out$loo <- loo::loo(log_lik)
      }

      return(out)
    },

    # Conditional consecutive probability of extreme
    lambda_z_posterior = function(z_level, fit) {
      lambda_z <- function(z_level,
                           margin_par,
                           copula_par) {
        # 1. Transform physical z to uniform scale
        u <- margin$cdf(z_level, margin_par)
        u <- pmax(pmin(u, 1 - 1e-10), 1e-10)

        # 2. Copula diagonal C(u,u)
        c_uu <- copula$cdf_diag(u, copula_par)

        # 3. Conditional exceedance probability
        lambda <- (1 - 2 * u + c_uu) / (1 - u)

        return(lambda)
      }

      margin_draws <- fit$margin_draws
      copula_draws <- fit$copula_draws

      n_draws <- nrow(margin_draws)
      out <- numeric(n_draws)

      for (i in seq_len(n_draws)) {
        p_margin <- as.numeric(margin_draws[i, ])
        names(p_margin) <- colnames(margin_draws)

        p_copula <- as.numeric(copula_draws[i, ])
        names(p_copula) <- colnames(copula_draws)
        out[i] <- lambda_z(z_level, p_margin, p_copula)
      }

      out
    }
  )
}
