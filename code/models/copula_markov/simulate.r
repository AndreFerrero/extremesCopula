source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copulas/joe.r")


egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)
egpd_joe_model <- make_copula_markov_model(margin_egpd, copula_joe, stan_mod = NULL)


margin_param <- c(mu = 0, kappa = 6, sigma = 1, xi = 0.3)
copula_param <- 1.1
n <- 100000

egpd_gumbel_data <- egpd_gumbel_model$simulate(
    n = n,
    margin_param = margin_param,
    copula_param = copula_param,
    seed = NULL
)
