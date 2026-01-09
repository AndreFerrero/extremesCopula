# ==============================================================
# Maxima and High-Quantile Analysis Script
# ==============================================================
# This script simulates data from a copula + lognormal margin model
# and studies the distribution of maxima and high quantiles
# for different copulas and parameters, compared to independence.
# ==============================================================

library(ggplot2)
library(dplyr)

# ------------------------------
# 1. Source your libs
# ------------------------------
source("libs/models/copulas/gumbel.R")     
source("libs/models/copulas/joe.r")    
source("libs/models/margins/lognormal.R")  
source("libs/models/builders/simulator.R") 

param_map <- list(
  margin = c("mu", "sigma"),
  copula = "theta"
)

# ------------------------------
# 2. Simulator wrapper
# ------------------------------
build_and_simulate <- function(copula, margin, param, n) {
  simulator <- build_simulator(copula, margin, param_map)
  X <- simulator(param, n)
  if (all(is.na(X))) stop("Simulator returned NA: check parameters")
  return(X)
}

# ------------------------------
# 3. Function to compute maxima and high quantiles
# ------------------------------
analyze_maxima_quantiles <- function(X, quantiles = c(0.90, 0.95, 0.99), top_k = 5) {
  
  # maxima
  max_X <- max(X)
  
  # top-k
  top_vals <- sort(X, decreasing = TRUE)[1:top_k]
  
  # high quantiles
  q_vals <- quantile(X, probs = quantiles)
  
  # return as list
  list(
    max = max_X,
    mean_top_k = mean(top_vals),
    top_k = top_vals,
    high_quantiles = q_vals
  )
}

# ------------------------------
# 4. Parameters for simulation
# ------------------------------
copula_list <- list(
  gumbel = copula_gumbel,
  joe    = copula_joe
)

theta_vec <- c(1, 1.5, 2, 3)
n_reps <- 200
n_obs <- 50000

param_base <- c(mu = 0, sigma = 1, theta = 1)

# ------------------------------
# 5. Run simulations and store maxima/quantiles
# ------------------------------
results <- list()

set.seed(123)

for(cop_name in names(copula_list)) {
  cat("Simulating for copula:", cop_name, "...\n")
  cop <- copula_list[[cop_name]]
  
  for(th in theta_vec) {
    replicate_stats <- replicate(n_reps, {
      param_base["theta"] <- th
      X <- build_and_simulate(cop, margin_lognormal, param_base, n_obs)
      analyze_maxima_quantiles(X)
    }, simplify = FALSE)
    
    # convert to dataframe
    df <- do.call(rbind, lapply(replicate_stats, function(res) {
      data.frame(
        max = res$max,
        mean_top_k = res$mean_top_k,
        q90 = res$high_quantiles[1],
        q95 = res$high_quantiles[2],
        q99 = res$high_quantiles[3]
      )
    }))
    df$theta <- th
    df$copula <- cop_name
    results[[paste0(cop_name, "_", th)]] <- df
  }
}

# Combine all results
results_df <- do.call(rbind, results)

# ------------------------------
# 6. Summarize
# ------------------------------
summary_df <- results_df %>%
  group_by(copula, theta) %>%
  summarize(
    max_med = median(max),
    max_low = quantile(max, 0.1),
    max_high = quantile(max, 0.9),
    mean_top_k_med = median(mean_top_k),
    q95_med = median(q95),
    q99_med = median(q99)
  )

print(summary_df)

# ------------------------------
# 7. Plot max and high quantiles vs copula parameter
# ------------------------------
# Maxima
ggplot(summary_df, aes(x = theta, y = max_med, color = copula)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = max_low, ymax = max_high), width = 0.1) +
  labs(title = "Median maximum vs copula parameter", x = "Copula parameter theta", y = "Maximum") +
  theme_minimal()

# Mean top-k
ggplot(summary_df, aes(x = theta, y = mean_top_k_med, color = copula)) +
  geom_line() + geom_point() +
  labs(title = "Mean top-k values vs copula parameter", x = "Copula parameter theta", y = "Mean top-k") +
  theme_minimal()

# q99 quantile
ggplot(summary_df, aes(x = theta, y = q99_med, color = copula)) +
  geom_line() + geom_point() +
  labs(title = "99% quantile vs copula parameter", x = "Copula parameter theta", y = "Quantile q99") +
  theme_minimal()
