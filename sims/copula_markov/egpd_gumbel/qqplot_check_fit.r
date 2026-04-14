source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")

sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

true_margin <- c(mu = 0, kappa = 6, sigma = 1, xi = 0.1)

true_theta <- 3.0
n_sim <- 1000

data <- sim_model$simulate(
  n = n_sim,
  copula_param = true_theta,
  margin_param = true_margin,
  seed = 123
)$x

fit_res_egpd_gumbel <- fit_egpd_gumbel(data)

# 1. Calculate the PIT values using the estimated parameters
pits <- get_pit_values(
  x = data, 
  theta_vec = fit_res_egpd_gumbel$par, 
  h_dist_fn = copula_gumbel$h_dist
)

# 2. Generate the plots
plot_diag_plots(pits)

# 3. (Optional) Statistical Check
# If the model is good, these p-values should be > 0.05
ks_marginal <- ks.test(pits$u, "punif")
ks_conditional <- ks.test(pits$w, "punif")

message("KS Test p-value (Marginal): ", round(ks_marginal$p.value, 4))
message("KS Test p-value (Conditional): ", round(ks_conditional$p.value, 4))