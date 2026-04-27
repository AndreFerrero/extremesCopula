source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copulas/joe.r")


egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)
egpd_joe_model <- make_copula_markov_model(margin_egpd, copula_joe, stan_mod = NULL)


margin_param <- c(mu = 0, kappa = 3, sigma = 1, xi = 0.1)
copula_param <- 2
n <- 1000

egpd_gumbel_data <- egpd_gumbel_model$simulate(
    n = n,
    margin_param = margin_param,
    copula_param = copula_param,
    seed = 46
)

egpd_data <- egpd::regpd(
    n,
    kappa = margin_param["kappa"],
    sigma = margin_param["sigma"],
    xi = margin_param["xi"]
)

# egpd_joe_data <- egpd_joe_model$simulate(
#     n = n,
#     margin_param = margin_param,
#     copula_param = copula_param,
#     seed = NULL
# )