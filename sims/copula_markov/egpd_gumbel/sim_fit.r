source("code/models/copula_markov/load_models.r")

# ---- Simulate ----
gumbel_sim <- gumbel_model$simulate(
  n = 1000,
  copula_param = 2,
  margin_param = c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2),
  seed = 46
)

gumbel_egpd2_sim <- gumbel_egpd2_model$simulate(
  n = 1000,
  copula_param = 2,
  margin_param = c(mu = 0, p = 0.5, kappa1 = 2, kappa2 = 2, sigma = 2, xi = 0.2),
  seed = 46
)

gumbel_egpd4_sim <- gumbel_egpd4_model$simulate(
  n = 1000,
  copula_param = 2,
  margin_param = c(mu = 0, kappa = 2, delta = 2, sigma = 2, xi = 0.2),
  seed = 46
)

gaussian_sim <- gaussian_model$simulate(
  n = 1000,
  copula_param = 0.7,
  margin_param = c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2),
  seed = 46
)

joe_sim <- joe_model$simulate(
  n = 1000,
  copula_param = 2,
  margin_param = c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2),
  seed = 46
)

# Independent Gumbel

gumbel_indep_sim <- gumbel_model$simulate(
  n = 5000,
  copula_param = 1,
  margin_param = c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.3),
  seed = 46
)

indep_gumbel_fit <- gumbel_model$fit(gumbel_indep_sim$x)

# ---- Model fit ----

gumbel_egpd2_fit <- gumbel_egpd2_model$fit(gumbel_egpd2_sim$x,
  I = 32,
  run_ppc = 1,
  iter = 2000,
)

gumbel_egpd4_fit <- gumbel_egpd4_model$fit(gumbel_egpd4_sim$x,
  I = 32,
  run_ppc = 1,
  iter = 1000,
)


### Gumbel data fit
gumbel_fit <- gumbel_model$fit(gumbel_sim$x,
  I = 32,
  run_ppc = 1,
  run_ei_mcmc = 1
)

gumbeljoe_fit <- joe_model$fit(gumbel_sim$x,
  run_ppc = 1,
  iter = 1000,
  I = 25,
  run_ei_mcmc = 1
)

### Gaussian data Fit
ggumbel_fit <- gumbel_model$fit(gaussian_sim$x,
  I = 32,
  run_ppc = 1
)

gaussian_fit <- gaussian_model$fit(gaussian_sim$x,
  run_ppc = 1
)

gjoe_fit <- joe_model$fit(gaussian_sim$x,
  I = 32,
  run_ppc = 1
)

### Joe data fit

joe_fit <- joe_model$fit(joe_sim$x,
  run_ppc = 1,
  iter = 1000,
  I = 25,
  run_ei_mcmc = 1
)

jgumbel_fit <- gumbel_model$fit(joe_sim$x,
  iter = 1000,
  run_ppc = 1,
  I = 25,
  run_ei_mcmc = 1
)
