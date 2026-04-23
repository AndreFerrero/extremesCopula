library(egpd)
library(evd)
library(ggplot2)
library(tidyr)
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")


egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

# --- 2. FITTING COMPONENTS ---

# A. Proposed Model: DS-EGPD (Filter then Fit EGPD)
fit_ds_egpd <- function(x) {
  nll_scale <- function(pars, x) {
    w <- pars[1]; a <- pars[2]; b <- pars[3]
    if(abs(b) >= 1) return(1e10)
    sig <- numeric(length(x))
    sig[1] <- mean(x) # Initial guess
    for(t in 2:length(x)) sig[t] <- exp(w + a*log(x[t-1]) + b*log(sig[t-1]))
    z <- x / sig
    -sum(dlnorm(z, 0, 1, log = TRUE))
  }
  
  opt_scale <- optim(c(0.1, 0.2, 0.5), nll_scale, x = x)
  w_h <- opt_scale$par[1]; a_h <- opt_scale$par[2]; b_h <- opt_scale$par[3]
  
  sig_h <- numeric(length(x))
  sig_h[1] <- mean(x)
  for(t in 2:length(x)) sig_h[t] <- exp(w_h + a_h*log(x[t-1]) + b_h*log(sig_h[t-1]))
  z_h <- x / sig_h
  
  fit_z <- tryCatch(egpd::fitegpd(z_h, type = 1, family = "egpd"), error = function(e) NULL)
  if(is.null(fit_z)) return(NA)
  return(fit_z$estimate["xi"])
}

# B. Filtered GPD (95% Threshold)
fit_filtered_gpd <- function(x, threshold_q = 0.95) {
  # Reuse filter logic
  nll_scale <- function(pars, x) {
    w <- pars[1]; a <- pars[2]; b <- pars[3]
    sig <- numeric(length(x))
    sig[1] <- mean(x)
    for(t in 2:length(x)) sig[t] <- exp(w + a*log(x[t-1]) + b*log(sig[t-1]))
    -sum(dlnorm(x/sig, 0, 1, log = TRUE))
  }
  opt <- optim(c(0.1, 0.2, 0.5), nll_scale, x = x)
  sig_h <- numeric(length(x)); sig_h[1] <- mean(x)
  for(t in 2:length(x)) sig_h[t] <- exp(opt$par[1] + opt$par[2]*log(x[t-1]) + opt$par[3]*log(sig_h[t-1]))
  z_h <- x / sig_h
  
  u <- quantile(z_h, threshold_q)
  fit <- tryCatch(extRemes::fevd(x, threshold = u, type = "GP"), error = function(e) NULL)
  if(is.null(fit)) return(NA)
  return(fit$results$par["shape"])
}

# --- 3. THE SIMULATION STUDY ---

set.seed(42)
n_sim <- 50
sample_size <- 2000
true_xi <- 0.3
true_kappa <- 2.0
true_sigma <- 1.0
true_theta <- 5.0 # Strong Gumbel Dependence

results <- data.frame()

for(i in 1:n_sim) {
  # Generate Data from the COPULA MARKOV model
  x_obs <- egpd_gumbel_model$simulate(
    n = sample_size,
    margin_param = c(mu = 0, sigma=true_sigma, xi=true_xi, kappa=true_kappa),
    copula_param = true_theta
  )$x
  
  # 1. Proposed: DS-EGPD (Wrong filter, but cleaning)
  xi_ds_egpd <- fit_ds_egpd(x_obs)
  
  # 2. IID EGPD (Ignores Gumbel Clustering)
  xi_iid_egpd <- egpd::fitegpd(x_obs, type = 1, family = "egpd")$estimate["xi"]
  
  # 3. Raw GPD (95% threshold)
  xi_raw_gpd <- extRemes::fevd(x_obs, threshold = quantile(x_obs, 0.95), type = "GP")$results$par["shape"]
  
  # 4. Filtered GPD (DS-GPD)
  xi_ds_gpd <- fit_filtered_gpd(x_obs)
  
  results <- rbind(results, data.frame(
    sim = i,
    DS_EGPD = xi_ds_egpd,
    IID_EGPD = xi_iid_egpd,
    GPD = xi_raw_gpd,
    DS_GPD = xi_ds_gpd
  ))
  if(i %% 20 == 0) message(paste("Iteration", i))
}

# --- 4. VISUALIZATION ---
res_long <- pivot_longer(results, cols = -sim, names_to = "Model", values_to = "xi_est")

ggplot(res_long, aes(x = Model, y = xi_est, fill = Model)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = true_xi, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Misspecification Robustness: Gumbel-Markov Data fitted with Scale Filter",
       subtitle = paste0("True xi = ", true_xi, " | Gumbel theta = ", true_theta, " | Filter is log-linear"),
       y = "Estimated xi") +
  theme_minimal()