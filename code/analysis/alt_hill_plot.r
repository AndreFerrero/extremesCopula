n_sequence = c(250, 500, 1000, 2000, 6000, 8000)
beta = 0.7

kn_beta <- sapply(n_sequence, FUN = function(n) floor(n^beta))
sapply(n_sequence, FUN = function(n) n^beta/n * 100)



# 1. Generate one large "ideal" sample from your simulation
library(ggplot2)
library(tidyr)
library(dplyr)

calc_alt_hill_data <- function(x) {
  n <- length(x)
  order_stats <- sort(x, decreasing = TRUE)
  
  # --- Step A: Estimate Rho Star once ---
  search_upper <- floor(min(n - 1, (2 * n) / log(log(n))))
  search_range <- seq(search_upper, 2, by = -1) # Start at 2 to avoid threshold error
  
  get_moments <- function(k) {
    # Ensure threshold exists at k+1
    log_excess <- log(order_stats[1:k]) - log(order_stats[k+1])
    c(M1 = mean(log_excess), M2 = mean(log_excess^2), 
      M3 = mean(log_excess^3), M4 = mean(log_excess^4))
  }
  
  rho_star <- NA
  for (k_val in search_range) {
    m <- get_moments(k_val)
    num <- (m["M4"] - 24 * m["M1"]^4) * (m["M2"] - 2 * m["M1"]^2)
    den <- m["M3"] - 6 * m["M1"]^3
    Sk2 <- (3 / 4) * (num / den^2)
    if (!is.na(Sk2) && Sk2 > 2/3 && Sk2 < 3/4) {
      rho_star <- (-4 + 6 * Sk2 + sqrt(3 * Sk2 - 2)) / (4 * Sk2 - 3)
      break
    }
  }
  
  if (is.na(rho_star)) {
    warning("Rho_star could not be estimated. Only standard Hill will be plotted.")
    rho_star <- 0 # Fallback to no correction
  } else {
    message(paste("Using Fixed Rho:", round(rho_star, 4)), " k =", k_val)
  }

  # --- Step B: Calculate Estimates for ALL k ---
  # We use a loop that explicitly creates a dataframe to ensure column names are kept
  results_list <- list()
  k_vec <- 5:(n-2) # Start at 5 to allow stable moment calculation
  
  term_rho <- rho_star / (1 - rho_star)
  
  for (i in seq_along(k_vec)) {
    k <- k_vec[i]
    log_excess <- log(order_stats[1:k]) - log(order_stats[k+1])
    m1 <- mean(log_excess)
    m2 <- mean(log_excess^2)
    
    # Correction logic: handle potential division by zero
    bias_corr <- if (m1 > 0 && term_rho != 0) {
      (m2 - 2 * m1^2) / (2 * m1 * term_rho)
    } else {
      0
    }
    
    results_list[[i]] <- data.frame(
      k = k,
      beta = log(k)/log(n),
      Hill = m1,
      Hill_BC = m1 - bias_corr
    )
  }
  
  # Combine and pivot
  dplyr::bind_rows(results_list) %>%
    tidyr::pivot_longer(
      cols = c("Hill", "Hill_BC"), 
      names_to = "Estimator", 
      values_to = "Estimate"
    )
}

margin_param <- c(mu = 0, kappa = 3, sigma = 1, xi = 0.1)
copula_param <- 2
n <- 8000

egpd_gumbel_data <- egpd_gumbel_model$simulate(
    n = n,
    margin_param = margin_param,
    copula_param = copula_param
)

# Run and Plot
plot_df <- calc_alt_hill_data(egpd_gumbel_data$x)

ggplot(plot_df, aes(x = beta, y = Estimate, color = Estimator)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 0.1, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("Hill" = "#377eb8", "Hill_BC" = "#e41a1c")) +
  labs(title = "Alternative Hill Plot (Resnick & Stărică Transformation)",
       subtitle = "Comparing Standard Hill vs. Moment-Based Bias Correction",
       x = expression(beta == log(k)/log(n)),
       y = expression(hat(xi))) +
  coord_cartesian(ylim = c(-0.1, 0.5)) + 
  theme_minimal()

ggplot(plot_df, aes(x = k, y = Estimate, color = Estimator)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 0.1, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("Hill" = "#377eb8", "Hill_BC" = "#e41a1c")) +
  labs(title = "Alternative Hill Plot (Resnick & Stărică Transformation)",
       subtitle = "Comparing Standard Hill vs. Moment-Based Bias Correction",
       x = expression(k),
       y = expression(hat(xi))) +
  coord_cartesian(ylim = c(-0.1, 0.5)) + 
  theme_minimal()
