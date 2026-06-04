source("code/models/copula_markov/copula_markov_model.r")

simulate_mixture_copula_markov <- function(
  n,
  copula1, copula2,
  pi,
  param1, param2,
  margin, margin_param,
  bisec_it = 20
) {

  u <- numeric(n)
  u[1] <- runif(1, 1e-5, 1 - 1e-5)

  for (t in 2:n) {

    v_prev <- u[t - 1]
    w <- runif(1)

    # latent regime
    z <- rbinom(1, 1, pi)

    if (z == 1) {
      u[t] <- copula_step(
        w, v_prev,
        copula1, param1,
        bisec_it
      )
    } else {
      u[t] <- copula_step(
        w, v_prev,
        copula2, param2,
        bisec_it
      )
    }
  }

  x <- margin$quantile(u, margin_param)

  list(x = x, u = u)
}

mix_data <- simulate_mixture_copula_markov(
  n = 1000,
  copula1 = copula_gumbel, copula2 = copula_gaussian,
  pi = 0.1,
  param1 = 2, param2 = 0.4,
  margin = margin_egpd, margin_param = c(mu = 0, sigma = 1, xi = 0.1, kappa = 2)
)

plot(mix_data$x, type = "l")
