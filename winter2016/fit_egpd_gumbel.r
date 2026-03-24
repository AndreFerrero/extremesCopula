source("code/models/copula_markov/load_models.r")
source("winter2016/load_data.r")

# Fit the gumbel copula Markov model to the data
fit_gumbel <- gumbel_model$fit(
    x = data,
    iter = 1000,
    chains = 4,
    seed = 46,
    I = 25,
    run_ppc = 1,
)

print(fit_gumbel$fit, pars = c("mu", "kappa", "sigma", "xi", "theta"))

save(fit_gumbel, file = "winter2016/winter2016fit_gumbel_egpd_mcmc1000_seed46_ei5000_chains4.rdata")

mcmc_trace(fit_gumbel$fit, pars = c("mu", "kappa", "sigma", "xi", "theta"))
mcmc_acf(fit_gumbel$fit, pars = c("mu", "kappa", "sigma", "xi", "theta"))

plot(fit_gumbel$margin_draws[, 4], type = "l")
acf(fit_gumbel$margin_draws[, 2])

ppc_stat(data, fit_gumbel$ppc, stat = "max") +
  ggtitle("Posterior Predictive Check: Maximum Values")

q90 <- function(x) quantile(x, probs = 0.90)
ppc_stat(data, fit_gumbel$ppc, stat = "q90")

