source("code/models/copula_markov/load_models.r")
source("winter2016/load_data.r")

init_fun <- function(x, chain_id, seed = 46) {
  set.seed(seed + chain_id) # Ensure different initial values for each chain
  repeat {
    sigma <- sd(x) * runif(1, 0.9, 1.1)
    xi <- runif(1, 0.1, 0.2)
    upper_bound <- sigma / xi
    if (max(x) < upper_bound) break
  }
  list(sigma = sigma, xi = xi, kappa = runif(1, 1.8, 2.2), thetam1 = runif(1, 0.9, 1.1))
}

n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) init_fun(data, chain_id = id))

# Fit the gumbel copula Markov model to the data

fit_gumbel <- gumbel_model$fit(
	x = data,
	iter = 1000,
	chains = n_chains,
	seed = 46,
	I = 25,
	init_par = init_ll,
	run_ppc = 1,
	extract_df = TRUE
)

fit_gumbel_egpd_noshift <- gumbel_egpd_noshift_model$fit(
	x = data - min(data) + 1e-8,
	iter = 1000,
	chains = n_chains,
	seed = 46,
	I = 25,
	init_par = init_ll,
	run_ppc = 1
)

mxi_fit_gumbel <- mxi_gumbel_model$fit(
  x = data,
  iter = 1000,
  chains = n_chains,
  seed = 46,
  I = 25,
  init_par = init_ll,
  run_ppc = 1
)

fit_gumbel_egpd4_noshift <- gumbel_egpd4_noshift_model$fit(
  x = data - min(data) + 1e-8,
  iter = 1000,
  chains = n_chains,
  seed = 46,
  I = 25,
  run_ppc = 1,
)

save(fit_gumbel_egpd4_noshift, file = "winter2016/winter2016fit_gumbel_egpd4_noshift_mcmc1000_seed46_chains4.rdata")

load("winter2016/winter2016fit_gumbel_egpd_noshift_mcmc1000_seed46_chains4.rdata")


FIT_RES <- fit_gumbel_egpd_noshift
params <- c("kappa", "sigma", "xi", "theta")


print(FIT_RES$fit, pars = params)

mcmc_trace(FIT_RES$fit, pars = params)
mcmc_acf(FIT_RES$fit, pars = params)

pairs(FIT_RES$fit, pars = params)

ppc_dens_overlay(data - min(data), FIT_RES$ppc[1:500, ]) +
	ggtitle("Posterior Predictive Check: Density Overlay")

ppc_stat(data - min(data), FIT_RES$ppc, stat = "max") +
  ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(data - min(data), FIT_RES$ppc, stat = "min") +
  ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(data - min(data), FIT_RES$ppc, stat = "median") +
  ggtitle("Posterior Predictive Check: Median Values")

ppc_stat(data - min(data), FIT_RES$ppc, stat = "sd") +
  ggtitle("Posterior Predictive Check: SD Values")

q90 <- function(x) quantile(x, probs = 0.90)
ppc_stat(data - min(data), FIT_RES$ppc, stat = "q90")

q10 <- function(x) quantile(x, probs = 0.10)
ppc_stat(data - min(data), FIT_RES$ppc, stat = "q10")
