# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Load all source files
load_all <- function() {
  source("R/copula_generators.R")
  source("R/gev_functions.R")
  source("R/asymptotic_model.R")
  source("R/estimation_dmle.R")
  source("R/estimation_mle.R")
  source("R/simulation.R")
  source("R/diagnostics.R")
}

#' Print estimation results
print_results <- function(fit, theta_true = NULL) {

  cat("\n====================================\n")
  cat("ESTIMATION RESULTS\n")
  cat("====================================\n")

  if (!is.null(theta_true)) {
    cat("theta_true =", theta_true, "\n")
  }

  cat("theta_hat  =", fit$theta, "\n")
  cat("mu_hat     =", fit$mu, "\n")
  cat("sigma_hat  =", fit$sigma, "\n")
  cat("xi_hat     =", fit$xi, "\n")

  if (!is.null(fit$convergence)) {
    cat("convergence:", fit$convergence, "\n")
  }
}
