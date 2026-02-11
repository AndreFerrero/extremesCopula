# Markov Simulation Engine using Bisection for Copula Inversion

simulate_copula_markov <- function(n, copula, theta, margin, margin_param, burn = 100, bisec_it = 20) {
  total_n <- n + burn
  U <- numeric(total_n)
  U[1] <- runif(1, 1e-5, 1 - 1e-5)

  # Bisection algo: inverting conditional distribution function
  for (t in 2:total_n) {
    v_prev <- U[t - 1]
    w_target <- runif(1)
    low <- 1e-10
    high <- 1 - 1e-10
    for (i in 1:bisec_it) {
      mid <- (low + high) / 2
      if (copula$h_dist(mid, v_prev, theta) < w_target) low <- mid else high <- mid
    }
    U[t] <- (low + high) / 2
  }
  U_final <- U[(burn + 1):total_n]
  X <- margin$quantile(U_final, margin_param)
  return(list(X = X, U = U_final, n_failures = 0))
}