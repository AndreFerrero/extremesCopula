library(bayesplot, help, pos = 2, lib.loc = NULL)
library(ggplot2, help, pos = 2, lib.loc = NULL)

#######
# Bayesian GUMBEL + EGPD Power
#######

init_fun <- function(x, chain_id, seed = 46) {
  set.seed(seed + chain_id) # Ensure different initial values for each chain
  repeat {
    sigma <- sd(x) * runif(1, 1, 5)
    xi <- runif(1, 0.05, 0.2)
    upper_bound <- sigma / xi
    if (max(x) < upper_bound) break
  }
  list(sigma = sigma, xi = xi, kappa = runif(1, 1.8, 5), thetam1 = runif(1, 0.9, 7))
}

n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) init_fun(winter_hourly_gust, chain_id = id))

egpd_gumbel_bayes_fit <- gumbel_egpd_noshift_model$fit(
  winter_hourly_gust,
  iter = 2000,
  run_ppc = 1,
  adapt_delta = 0.9,
  I = 16,
  init_par = init_ll
)

# save(egpd_gumbel_bayes_fit, file = "ERA5/egpd_gumbel_bayesfit.RData")
load("ERA5/egpd_gumbel_bayesfit.RData")

params <- c("kappa", "sigma", "xi", "theta")

print(egpd_gumbel_bayes_fit$fit, pars = params)

mcmc_trace(
  egpd_gumbel_bayes_fit$fit,
  pars = params
)

mcmc_acf(egpd_gumbel_bayes_fit$fit, pars = params)

pairs(egpd_gumbel_bayes_fit$fit, pars = params)

draws <- rstan::extract(egpd_gumbel_bayes_fit$fit, pars = params)

cor(cbind(draws$kappa, draws$sigma, draws$xi, draws$theta))


ppc_dens_overlay(winter_hourly_gust, egpd_gumbel_bayes_fit$ppc[1:100, ]) +
  ggtitle("Posterior Predictive Check: Density Overlay")

ppc_stat(winter_hourly_gust, egpd_gumbel_bayes_fit$ppc, stat = "max") +
  ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(winter_hourly_gust, egpd_gumbel_bayes_fit$ppc, stat = "min") +
  ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(winter_hourly_gust, egpd_gumbel_bayes_fit$ppc, stat = "median") +
  ggtitle("Posterior Predictive Check: Median Values")

ppc_stat(winter_hourly_gust, egpd_gumbel_bayes_fit$ppc, stat = "sd") +
  ggtitle("Posterior Predictive Check: SD Values")

q90 <- function(x) quantile(x, probs = 0.90)
ppc_stat(winter_hourly_gust, egpd_gumbel_bayes_fit$ppc, stat = "q90")

q10 <- function(x) quantile(x, probs = 0.10)
ppc_stat(winter_hourly_gust, egpd_gumbel_bayes_fit$ppc, stat = "q10")

## Bayesian p-values
winter_max <- winter_hourly_gust |> max()
winter_q90 <- quantile(winter_hourly_gust, probs = 0.90)

max_ppc <- apply(egpd_gumbel_bayes_fit$ppc, 1, max)
q90_ppc <- apply(egpd_gumbel_bayes_fit$ppc, 1, q90)

sum(max_ppc > winter_max) / length(max_ppc)
sum(q90_ppc > winter_q90) / length(q90_ppc)


# cluster of extremes posterior predictive distribution

#---------------------------------------
# Inputs:
#---------------------------------------
ppc <- egpd_gumbel_bayes_fit$ppc # posterior predictive: rows=MCMC draws, cols=time points
u <- 0.95
block_size <- 24 # e.g., yearly block for block maxima (adjust as needed)
n_draws <- nrow(ppc)
n_time <- ncol(ppc)

#---------------------------------------
# Storage for summaries
#---------------------------------------
cluster_sizes_list <- vector("list", n_draws)
prob_run3 <- numeric(n_draws)
prob_run8 <- numeric(n_draws)
joint_prob <- numeric(n_draws)
block_maxima <- vector("list", n_draws)

#---------------------------------------
# Loop over posterior predictive draws
#---------------------------------------

for (i in 1:n_draws) {
  sim <- ppc[i, ]

  #--- 1. Identify exceedances
  exceed <- sim > quantile(sim, probs = u)
  r <- rle(exceed)

  # Cluster sizes (only for TRUE clusters)
  cluster_sizes <- r$lengths[r$values]
  cluster_sizes_list[[i]] <- cluster_sizes

  # Probability of at least 3 consecutive exceedances
  prob_run3[i] <- as.numeric(any(cluster_sizes >= 3))
  prob_run8[i] <- as.numeric(any(cluster_sizes >= 8))

  # Joint exceedance probability P(X_t>u, X_{t+1}>u)
  joint_prob[i] <- mean(sim[-n_time] > quantile(sim, probs = u) & sim[-1] > quantile(sim, probs = u))

  #--- 2. Block maxima
  # Split into blocks
  n_blocks <- floor(n_time / block_size)
  block_max <- sapply(1:n_blocks, function(b) {
    start <- (b - 1) * block_size + 1
    end <- b * block_size
    max(sim[start:end])
  })
  block_maxima[[i]] <- block_max
}

#---------------------------------------
# Posterior summaries
#---------------------------------------

# Cluster sizes: mean, median, quantiles
all_clusters <- unlist(cluster_sizes_list)
cluster_summary <- quantile(all_clusters, probs = c(0.025, 0.5, 0.975))

# Probability of runs >=3
prob_run3_summary <- quantile(prob_run3, probs = c(0.025, 0.5, 0.975))

# Probability of runs >=3
prob_run8_summary <- quantile(prob_run8, probs = c(0.025, 0.5, 0.975))

# Joint exceedance probability
joint_prob_summary <- quantile(joint_prob, probs = c(0.025, 0.5, 0.975))

# Block maxima summary (e.g., return levels)
block_max_all <- unlist(block_maxima)
block_max_summary <- quantile(block_max_all, probs = c(0.025, 0.5, 0.975))

#---------------------------------------
# Output summaries
#---------------------------------------
list(
  cluster_summary = cluster_summary,
  prob_run3_summary = prob_run3_summary,
  prob_run8_summary = prob_run8_summary,
  joint_prob_summary = joint_prob_summary,
  block_max_summary = block_max_summary
)
