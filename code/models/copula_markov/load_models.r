# ---- Load packages ----
source("code/packages.r")

# ---- Load model components ----
source("code/models/margins/egp.r")
source("code/models/copulas/gumbel.r")
source("code/models/copulas/joe.r")
source("code/models/copulas/gaussian.r")
source("code/models/copula_markov/copula_markov_model.R")

rstan_options(auto_write = TRUE)
# gaussian_stan <- rstan::stan_model("code/stan/gaussian_egpd.stan")
gumbel_stan <- rstan::stan_model("code/stan/gumbel_egpd.stan")
# joe_stan <- rstan::stan_model("code/stan/joe_egpd.stan")

gumbel_model <- make_copula_markov_model(
  margin = margin_egp,
  copula = copula_gumbel,
  stan_mod = gumbel_stan
)

# gaussian_model <- make_copula_markov_model(
#   margin = margin_egp,
#   copula = copula_gaussian,
#   stan_mod = gaussian_stan
# )

# joe_model <- make_copula_markov_model(
#   margin = margin_egp,
#   copula = copula_joe,
#   stan_mod = joe_stan
# )

