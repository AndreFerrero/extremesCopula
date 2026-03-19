# ---- Load packages ----
source("code/packages.r")

# ---- Load model components ----
source("code/models/margins/egp.r")
source("code/models/copulas/gumbel.r")
source("code/models/copulas/joe.r")
source("code/models/copulas/gaussian.r")
source("code/models/copula_markov/copula_markov_model.R")

rstan_options(auto_write = TRUE)
gaussian_stan <- rstan::stan_model("code/stan/gaussian_egpd.stan")
gumbel_stan <- rstan::stan_model("code/stan/gumbel_egpd.stan")
joe_stan <- rstan::stan_model("code/stan/joe_egpd.stan")

gumbel_model <- make_copula_markov_model(
  margin = margin_egp,
  copula = copula_gumbel,
  stan_mod = gumbel_stan
)

gaussian_model <- make_copula_markov_model(
  margin = margin_egp,
  copula = copula_gaussian,
  stan_mod = gaussian_stan
)

joe_model <- make_copula_markov_model(
  margin = margin_egp,
  copula = copula_joe,
  stan_mod = joe_stan
)

# ---- Simulate ----
gumbel_sim <- gumbel_model$simulate(
  n = 1000,
  copula_param = 2,
  margin_param = c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2),
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
### Gumbel data fit
gumbel_fit <- gumbel_model$fit(gumbel_sim$x,
  I = 32,
  run_ppc = 1,
  ei_mcmc = 5000
)

gumbeljoe_fit <- joe_model$fit(gumbel_sim$x,
  run_ppc = 1,
  iter = 1000,
  I = 25,
  ei_mcmc = 5000
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
  ei_mcmc = 5000
)

jgumbel_fit <- gumbel_model$fit(joe_sim$x,
  iter = 1000,
  run_ppc = 1,
  I = 25,
  ei_mcmc = 5000
)
