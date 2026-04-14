source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/bernstein_egpd_gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")

sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

true_margin <- c(mu = 0, kappa = 6, sigma = 1, xi = 0.3)

true_theta <- 3.0
n_sim <- 5000

x <- sim_model$simulate(
  n = n_sim,
  copula_param = true_theta,
  margin_param = true_margin,
  seed = 123
)$x

m_degree <- 5

get_init_egpd_gumbel(x)
get_init_bern_egpd_gumbel(x, m_degree)

fit_res_egpd_gumbel <- fit_egpd_gumbel(x)
fit_bern_egpd_gumbel <- fit_egpd_bernstein_gumbel(x, m = m_degree)

res_kappa <- fit_res_egpd_gumbel$estimate["kappa"]
res_sigma <- fit_res_egpd_gumbel$estimate["sigma"]
res_xi <- fit_res_egpd_gumbel$estimate["xi"]
res_theta <- fit_res_egpd_gumbel$estimate["theta"]
res_kappa_bern <- fit_bern_egpd_gumbel$estimate["kappa"]
res_sigma_bern <- fit_bern_egpd_gumbel$estimate["sigma"]
res_xi_bern <- fit_bern_egpd_gumbel$estimate["xi"]
res_theta_bern <- fit_bern_egpd_gumbel$estimate["theta"]
res_weights <- fit_bern_egpd_gumbel$estimate[paste0("w", 1:m_degree)]

paste("Estimated kappa:", round(res_kappa, 3))
paste("Estimated sigma:", round(res_sigma, 3))
paste("Estimated xi:", round(res_xi, 3))
paste("Estimated theta:", round(res_theta, 3))

paste("Estimated kappa (Bernstein):", round(res_kappa_bern, 3))
paste("Estimated sigma (Bernstein):", round(res_sigma_bern, 3))
paste("Estimated xi (Bernstein):", round(res_xi_bern, 3))
paste("Estimated theta (Bernstein):", round(res_theta_bern, 3))

paste("Estimated Bernstein weights:", round(res_weights, 3))

# Generate x sequence for the smooth red line
x_seq <- seq(min(x), max(x), length.out = 500)
y_fit <- egpd:::degpd_density(
  x_seq,
  sigma = res_sigma, xi = res_xi, kappa = res_kappa, type = 1
)
y_bern_fit <- egpd:::.bernstein_full_density(
  x_seq, res_sigma_bern, res_xi_bern, res_kappa_bern, res_weights, m_degree
)

# Plot
hist(x,
  breaks = 40, prob = TRUE, col = "lightblue", border = "darkgray",
  main = "Histogram and Fitted Density", xlab = "x"
)
lines(x_seq, y_fit, col = "blue", lwd = 2)
lines(x_seq, y_bern_fit, col = "red", lwd = 2)

legend("topright", legend = c("EGPD-Gumbel Fit", "Bernstein EGPD-Gumbel Fit"),
       col = c("blue", "red"), lwd = 2)
