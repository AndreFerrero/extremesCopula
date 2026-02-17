# Optimized Markov Simulation Engine
simulate_copula_markov <- function(n, copula, theta, margin, margin_param, burn = 100, bisec_it = 20) {
  total_n <- n + burn
  U <- numeric(total_n)
  U[1] <- runif(1, 1e-5, 1 - 1e-5)

  # Check if the copula object has a pre-defined analytical inverse
  has_h_inv <- !is.null(copula$h_inv)

  for (t in 2:total_n) {
    v_prev <- U[t - 1]
    w_target <- runif(1)
    
    if (has_h_inv) {
      # FAST PATH: Use analytical inversion (e.g., Gaussian, Clayton)
      U[t] <- copula$h_inv(w_target, v_prev, theta)
    } else {
      # ROBUST PATH: Fallback to Bisection (e.g., Gumbel, Joe)
      low <- 1e-10
      high <- 1 - 1e-10
      for (i in 1:bisec_it) {
        mid <- (low + high) / 2
        if (copula$h_dist(mid, v_prev, theta) < w_target) low <- mid else high <- mid
      }
      U[t] <- (low + high) / 2
    }
  }
  
  U_final <- U[(burn + 1):total_n]
  X <- margin$quantile(U_final, margin_param)
  
  return(list(X = X, U = U_final, method_used = ifelse(has_h_inv, "analytical", "bisection")))
}