# libs/models/builders/semibsl_logposterior.R

source("libs/packages.R")

#' Build a SemiBSL log-posterior function
#'
#' @param simulator Function: simulates n_obs samples given param vector
#' @param sum_stats Function: computes summary statistics from simulated data
#' @param n_sim Number of simulations per likelihood evaluation
#' @param log_prior Function: returns log prior for param
#' @param inverse_transform function to transform parameters in the original space
#' @param log_jacobian jacobian function of the transformation
#'
#' @return Function(param) -> log-posterior
build_bsl_logposterior <- function(
  copula, margin, param_map,
  data,
  simulator,
  sum_stats,
  n_sim,
  synthetic_loglik,
  inverse_transform = NULL,
  log_jacobian = NULL,
  margin_prior = NULL,
  copula_prior = NULL
) {

  # Precompute observed summary statistics
  y_obs <- sum_stats(data)

  function(param_init) {

    # Apply inverse transform if given
    # param is always the parameter in the original space
    # no transformation leaves param = param_init
    param <- if (!is.null(inverse_transform)) inverse_transform(param_init) else param_init

    param_m <- param[param_map$margin]
    param_c <- param[param_map$copula]

    # --- Prior ---
    # do.call is to insert the optional prior parameters to overide the defaults
    logprior_m <- do.call(margin$log_prior, c(list(param_m), margin_prior))
    logprior_c <- do.call(copula$log_prior, c(list(param_c), copula_prior))
    logprior <- logprior_m + logprior_c
    
    if (!is.finite(logprior)) return(-Inf)

    # Simulate n_sim datasets and compute summary statistics
    sim_stats <- matrix(NA, nrow = n_sim, ncol = length(y_obs))
    for (i in seq_len(n_sim)) {
      x_sim <- simulator(param, length(data))
      if (any(is.na(x_sim))) return(-Inf)
      sim_stats[i, ] <- sum_stats(x_sim)
    }
    
    loglik <- synthetic_loglik(sim_stats, y_obs)

    # Jacobian adjustment if in transformed space
    logjac <- if (!is.null(log_jacobian)) log_jacobian(param_init) else 0

    logprior + loglik + logjac
  }
}
