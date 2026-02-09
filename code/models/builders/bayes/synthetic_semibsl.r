synthetic_semibsl <- function(sim_stats, y_obs) {

    n_sim <- nrow(sim_stats)
    d <- ncol(sim_stats)

    loglik_marg <- 0
    eta_obs <- numeric(d)
    eta_sim <- matrix(NA, n_sim, d)

    for (j in 1:d) {

      sd_sim <- sd(sim_stats[, j])
      iqr_sim <- IQR(sim_stats[, j])
      h <- 0.9 * min(sd_sim, iqr_sim / 1.34) * n_sim^(-0.2)
      if (h == 0) h <- 1e-6

      z <- (y_obs[j] - sim_stats[, j]) / h
      log_vals <- dnorm(z, log = TRUE)

      loglik_marg <- loglik_marg +
        max(log_vals) +
        log(sum(exp(log_vals - max(log_vals)))) -
        log(n_sim) - log(h)

      u_obs <- mean(pnorm(z))
      u_obs <- min(max(u_obs, 1e-6), 1 - 1e-6)
      eta_obs[j] <- qnorm(u_obs)

      u_sim <- rank(sim_stats[, j]) / (n_sim + 1)
      eta_sim[, j] <- qnorm(pmin(pmax(u_sim, 1e-6), 1 - 1e-6))
    }

    R_hat <- cor(eta_sim)
    R_hat <- R_hat + diag(1e-6, d)

    R_inv <- solve(R_hat)
    log_det_R <- determinant(R_hat, logarithm = TRUE)$modulus
    quad <- t(eta_obs) %*% (R_inv - diag(d)) %*% eta_obs

    loglik_marg - 0.5 * log_det_R - 0.5 * quad
  }
