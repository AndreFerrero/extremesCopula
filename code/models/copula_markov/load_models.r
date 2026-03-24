# ---- Load packages ----
source("code/packages.r")

# ---- Load model components ----
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copulas/joe.r")
source("code/models/copulas/gaussian.r")
source("code/models/copula_markov/copula_markov_model.R")

rstan_options(auto_write = TRUE)
# gaussian_stan <- rstan::stan_model("code/stan/gaussian_egpd.stan")

gumbel_stan <- rstan::stan_model("code/stan/gumbel_egpd.stan")
# mxi_gumbel_stan <- rstan::stan_model("code/stan/mxi_gumbel_egpd.stan")

gumbel_egpd2_stan <- rstan::stan_model("code/stan/gumbel_egpd2.stan")
gumbel_egpd4_stan <- rstan::stan_model("code/stan/gumbel_egpd4.stan")

joe_stan <- rstan::stan_model("code/stan/joe_egpd.stan")

gumbel_model <- make_copula_markov_model(
  margin = margin_egpd,
  copula = copula_gumbel,
  stan_mod = gumbel_stan
)

# mxi_gumbel_model <- make_copula_markov_model(
#   margin = margin_egp,
#   copula = copula_gumbel,
#   stan_mod = mxi_gumbel_stan
# )

gumbel_egpd2_model <- make_copula_markov_model(
  margin = margin_egpd2,
  copula = copula_gumbel,
  stan_mod = gumbel_egpd2_stan
)

gumbel_egpd4_model <- make_copula_markov_model(
  margin = margin_egpd4,
  copula = copula_gumbel,
  stan_mod = gumbel_egpd4_stan
)

# gaussian_model <- make_copula_markov_model(
#   margin = margin_egpd,
#   copula = copula_gaussian,
#   stan_mod = gaussian_stan
# )

joe_model <- make_copula_markov_model(
  margin = margin_egpd,
  copula = copula_joe,
  stan_mod = joe_stan
)

