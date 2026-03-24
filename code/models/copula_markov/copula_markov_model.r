sample_bisection <- function(w, v_prev, fun, copula_param, iterations) {
  low <- 1e-10
  high <- 1 - 1e-10

  for (i in 1:iterations) {
    mid <- (low + high) / 2
    if (fun(mid, v_prev, copula_param) < w) {
      low <- mid
    } else {
      high <- mid
    }
  }
  return((low + high) / 2)
}

simulate_copula_markov <- function(
  n,
  copula, copula_param,
  margin, margin_param,
  bisec_it = 20
) {
  u <- numeric(n)
  u[1] <- runif(1, 1e-5, 1 - 1e-5)

  # Check if the copula object has a pre-defined analytical inverse
  has_h_inv <- !is.null(copula$h_inv)

  for (t in 2:n) {
    v_prev <- u[t - 1]
    w_target <- runif(1)

    if (has_h_inv) {
      # Analytical Inverse
      u[t] <- copula$h_inv(w_target, v_prev, copula_param)
    } else {
      # Bisection Method
      u[t] <- sample_bisection(
        w_target, v_prev,
        copula$h_dist, copula_param,
        bisec_it
      )
    }
  }

  x <- margin$quantile(u, margin_param)

  return(list(x = x, u = u, method_used = ifelse(has_h_inv, "analytical", "bisection")))
}


make_copula_markov_model <- function(
  margin,
  copula,
  stan_mod
) {
  # ----------------------------------------
  # Build object
  # ----------------------------------------

  list(
    margin = margin,
    copula = copula,

    # -----------------------
    # SIMULATION
    # -----------------------
    simulate = function(n,
                        copula_param,
                        margin_param,
                        seed = NULL) {
      if (!is.null(seed)) set.seed(seed)

      out <- simulate_copula_markov(
        n = n,
        copula = copula,
        copula_param = copula_param,
        margin = margin,
        margin_param = margin_param
      )

      return(out)
    },

    # -----------------------
    # FIT
    # -----------------------
    fit = function(x,
                   iter = 2000,
                   chains = 4,
                   seed = 42,
                   prior_check = 0,
                   run_ppc = 0,
                   run_ei_mcmc = NULL,
                   I = NULL,
                   init_par = NULL,
                   adapt_delta = NULL,
                   extract_df = FALSE) {
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

      # Some models don't need I
      if (!is.null(I)) stan_data$I <- I

      if (is.null(init_par)) {
        init_par <- "random"
      }

      if (!is.null(run_ei_mcmc)) stan_data$run_ei_mcmc <- run_ei_mcmc

      message("Sampling...")

      fit_obj <- rstan::sampling(
        object = stan_mod,
        data = stan_data,
        iter = iter,
        chains = chains,
        seed = seed,
        init = init_par,
        control = control_list
      )

      out <- list(
        fit = fit_obj,
        x = x
      )

      if (extract_df) {
        # Extract all posterior draws
        posterior_list <- rstan::extract(fit_obj, permuted = TRUE)

        # Remove generated quantities we don't want
        exclude <- ifelse(run_ppc == 0,
          NULL,
          c("x_rep", "lp__")
        )

        param_names <- setdiff(names(posterior_list), exclude)

        # Convert to dataframe
        draws_df <- as.data.frame(posterior_list[param_names])

        margin_draws <- draws_df[, margin$name_param, drop = FALSE]
        copula_draws <- draws_df[, copula$name_param, drop = FALSE]

        out$margin_draws <- margin_draws
        out$copula_draws <- copula_draws
      }

      # PPC
      if (run_ppc == 1) {
        out$ppc <- rstan::extract(fit_obj, "x_rep")$x_rep
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
