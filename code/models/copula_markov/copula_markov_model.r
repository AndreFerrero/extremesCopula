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

copula_step <- function(w, v_prev, copula, copula_param, bisec_it = 20) {

  has_h_inv <- !is.null(copula$h_inv)

  if (has_h_inv) {
    return(copula$h_inv(w, v_prev, copula_param))
  }

  sample_bisection(
    w, v_prev,
    copula$h_dist,
    copula_param,
    bisec_it
  )
}

simulate_copula_markov <- function(
  n,
  copula, copula_param,
  margin, margin_param,
  bisec_it = 20
) {

  u <- numeric(n)
  u[1] <- runif(1, 1e-5, 1 - 1e-5)

  for (t in 2:n) {
    v_prev <- u[t - 1]
    w <- runif(1)

    u[t] <- copula_step(
      w, v_prev,
      copula, copula_param,
      bisec_it
    )
  }

  x <- margin$quantile(u, margin_param)

  list(x = x, u = u,
       method_used = ifelse(is.null(copula$h_inv), "bisection", "analytical"))
}

make_copula_markov_model <- function(
  margin,
  copula,
  stan_mod = NULL
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
    }
  )
}
