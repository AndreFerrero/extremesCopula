n <- 10000

joe_data <- joe_model$simulate(
  n = n,
  copula_param = 2,
  margin_param = c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2)
)

lag_plot(joe_data$u)
lag_plot(joe_data$x)
cor(joe_data$u[-n], joe_data$u[-1], method = "kendall")

gumbel_data <- gumbel_model$simulate(
  n = n,
  copula_param = 2,
  margin_param = c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2)
)

lag_plot(gumbel_data$u)

cor(gumbel_data$u[-n], gumbel_data$u[-1], method = "kendall")
