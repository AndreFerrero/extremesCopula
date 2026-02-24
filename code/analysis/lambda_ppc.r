lambda_z_empirical <- function(x, prob) {
  if (length(x) < 2) {
    stop("x must have length >= 2")
  }

  z_level <- quantile(x, prob)

  x_t <- x[-length(x)]
  x_t1 <- x[-1]

  exceed_t <- x_t > z_level
  exceed_t1 <- x_t1 > z_level

  den <- sum(exceed_t)

  if (den == 0) {
    return(NA_real_)
  }

  num <- sum(exceed_t & exceed_t1)

  num / den
}

ppc_lambda <- function(x_obs, x_ppc, prob) {
  # Arguments:
  # x_obs: observed series
  # x_ppc: T x N_draws matrix, e.g., fit$ppc
  # z_level: threshold to compute lambda(z)

  n_draws <- ncol(x_ppc)
  lambda_ppc <- numeric(n_draws)

  for (i in seq_len(n_draws)) {
    lambda_ppc[i] <- lambda_z_empirical(x_ppc[i, ], prob)
  }

  # Observed lambda
  lambda_obs <- lambda_z_empirical(x_obs, prob)

  out <- list(
    ppc_stat = lambda_ppc,
    obs_stat = lambda_obs
  )

  return(out)
}
