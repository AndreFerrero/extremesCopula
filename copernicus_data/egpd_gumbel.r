source("copernicus_data/load_data.r")

source("winter2016/marg_dep_mods_funcs.R")

source("code/models/copula_markov/load_models.r")

library(bayesplot)

winter_hourly_gust <- data$fg10[data$season == "Winter"]

lag_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) |>
    na.omit()

  return(df)
}

winter_lag <- lag_df(winter_hourly_gust)

tau <- cor(winter_lag$x, winter_lag$x_lag, method = "kendall")

1 / (1 - tau)

###
# GUMBEL + GPD - Censored likelihood - Winter2016 implementation
###
u.thresh <- quantile(winter_hourly_gust, probs = 0.90) # Modelling threshold u on uniform scale

cens_gpd_gumbel <- FitGpd(dat = winter_lag, u = u.thresh)

gamma <- cens_gpd_gumbel$par[1]
sig.u <- cens_gpd_gumbel$par[2]
xi <- cens_gpd_gumbel$par[3]
lambda.u <- cens_gpd_gumbel$par[4]

gamma
1 / gamma
sig.u
xi
lambda.u

# Grid of thresholds
thresholds <- quantile(winter_hourly_gust, probs = seq(0.9, 0.97, by = 0.01))

n_rep <- 50

# Storage arrays
results <- list()

for (j in seq_along(thresholds)) {
  
  u <- thresholds[j]
  
  gamma_vec  <- numeric(n_rep)
  sigma_vec  <- numeric(n_rep)
  xi_vec     <- numeric(n_rep)
  lambda_vec <- numeric(n_rep)
  
  for (i in 1:n_rep) {
    
    fit <- FitGpd(dat = winter_lag, u = u)
    
    gamma_vec[i]  <- fit$par[1]
    sigma_vec[i]  <- fit$par[2]
    xi_vec[i]     <- fit$par[3]
    lambda_vec[i] <- fit$par[4]
  }
  
  results[[j]] <- data.frame(
    threshold = u,
    gamma_mean = mean(gamma_vec, na.rm = TRUE),
    gamma_sd   = sd(gamma_vec, na.rm = TRUE),
    
    sigma_mean = mean(sigma_vec, na.rm = TRUE),
    sigma_sd   = sd(sigma_vec, na.rm = TRUE),
    
    xi_mean    = mean(xi_vec, na.rm = TRUE),
    xi_sd      = sd(xi_vec, na.rm = TRUE),
    
    lambda_mean = mean(lambda_vec, na.rm = TRUE),
    lambda_sd   = sd(lambda_vec, na.rm = TRUE)
  )
}

# Combine results
res_df <- do.call(rbind, results)

plot(res_df$threshold, res_df$xi_mean, type = "b", pch = 16,
     ylim = range(res_df$xi_mean + 2*res_df$xi_sd,
                  res_df$xi_mean - 2*res_df$xi_sd),
     xlab = "Threshold", ylab = "xi",
     main = "Threshold Stability (xi) with stochastic fits")

arrows(res_df$threshold,
       res_df$xi_mean - 2*res_df$xi_sd,
       res_df$threshold,
       res_df$xi_mean + 2*res_df$xi_sd,
       angle = 90, code = 3, length = 0.05)

plot(res_df$threshold, res_df$sigma_mean, type = "b", pch = 16,
     ylim = range(res_df$sigma_mean + 2*res_df$sigma_sd,
                  res_df$sigma_mean - 2*res_df$sigma_sd),
     xlab = "Threshold", ylab = "sigma",
     main = "Threshold Stability (sigma)")

arrows(res_df$threshold,
       res_df$sigma_mean - 2*res_df$sigma_sd,
       res_df$threshold,
       res_df$sigma_mean + 2*res_df$sigma_sd,
       angle = 90, code = 3, length = 0.05)

plot(res_df$threshold, res_df$gamma_mean, type = "b", pch = 16,
     ylim = range(res_df$gamma_mean + 2*res_df$gamma_sd,
                  res_df$gamma_mean - 2*res_df$gamma_sd),
     xlab = "Threshold", ylab = "gamma",
     main = "Threshold Stability (dependence parameter)")

arrows(res_df$threshold,
       res_df$gamma_mean - 2*res_df$gamma_sd,
       res_df$threshold,
       res_df$gamma_mean + 2*res_df$gamma_sd,
       angle = 90, code = 3, length = 0.05)

plot(res_df$threshold, res_df$lambda_mean, type = "b", pch = 16,
     ylim = range(res_df$lambda_mean + 2*res_df$lambda_sd,
                  res_df$lambda_mean - 2*res_df$lambda_sd),
     xlab = "Threshold", ylab = "lambda",
     main = "Threshold Stability (exceedance rate)")

arrows(res_df$threshold,
       res_df$lambda_mean - 2*res_df$lambda_sd,
       res_df$threshold,
       res_df$lambda_mean + 2*res_df$lambda_sd,
       angle = 90, code = 3, length = 0.05)

###
# FULL bayesian model GUMBEL + EGPD
###

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

egpd_gumbel_fit <- gumbel_egpd_noshift_model$fit(
  winter_hourly_gust,
  iter = 2000,
  run_ppc = 1,
  adapt_delta = 0.9,
  I = 16,
	init_par = init_ll
)

params <- c("kappa", "sigma", "xi", "theta")

print(egpd_gumbel_fit$fit, pars = params)

mcmc_trace(
  egpd_gumbel_fit$fit,
  pars = params
)

mcmc_acf(egpd_gumbel_fit$fit, pars = params)

pairs(egpd_gumbel_fit$fit, pars = params)

ppc_dens_overlay(winter_hourly_gust, egpd_gumbel_fit$ppc[1:1000, ]) +
	ggtitle("Posterior Predictive Check: Density Overlay")

ppc_stat(winter_hourly_gust, egpd_gumbel_fit$ppc, stat = "max") +
  ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(winter_hourly_gust, egpd_gumbel_fit$ppc, stat = "min") +
  ggtitle("Posterior Predictive Check: Maximum Values")

ppc_stat(winter_hourly_gust, egpd_gumbel_fit$ppc, stat = "median") +
  ggtitle("Posterior Predictive Check: Median Values")

ppc_stat(winter_hourly_gust, egpd_gumbel_fit$ppc, stat = "sd") +
  ggtitle("Posterior Predictive Check: SD Values")

q90 <- function(x) quantile(x, probs = 0.90)
ppc_stat(winter_hourly_gust, egpd_gumbel_fit$ppc, stat = "q90")

q10 <- function(x) quantile(x, probs = 0.10)
ppc_stat(winter_hourly_gust, egpd_gumbel_fit$ppc, stat = "q10")
