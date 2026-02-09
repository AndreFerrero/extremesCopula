synthetic_gaussian_bsl <- function(sim_stats, y_obs) {

    mu_hat <- colMeans(sim_stats)
    Sigma_hat <- cov(sim_stats)

    if (any(!is.finite(Sigma_hat))) return(-Inf)

    mvtnorm::dmvnorm(y_obs, mu_hat, Sigma_hat, log = TRUE)
  }