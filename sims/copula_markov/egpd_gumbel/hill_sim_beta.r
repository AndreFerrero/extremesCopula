library(ggplot2)
library(dplyr)
library(tidyr)
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")

egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

# 1. Main Simulation Parameters
n_sim <- 2000
mc_iterations <- 100 # Adjust to 200 or 500 for smoother lines if time permits
kappa_values <- c(1, 3, 6)
true_xi <- 0.1
beta_sequence <- seq(0.4, 0.85, by = 0.02) # The range of threshold depths

# 2. Main Simulation Loop (Sequential)
all_results <- list()
counter <- 1

for (k_val in kappa_values) {
  cat(sprintf("\n--- Starting Simulation for Kappa = %d ---\n", k_val))
  
  # Matrix to store MC iterations for averaging
  # Rows = Beta levels, Cols = MC iterations
  hill_mc_matrix <- matrix(NA, nrow = length(beta_sequence), ncol = mc_iterations)
  bc_mc_matrix <- matrix(NA, nrow = length(beta_sequence), ncol = mc_iterations)
  
  for (m in 1:mc_iterations) {
    # Progress indicator
    if (m %% 10 == 0) cat(sprintf("  Iteration %d/%d\n", m, mc_iterations))
    
    # a. Generate Data
    set.seed(123 + k_val * 1000 + m) # Consistent seeds
    x <- egpd_gumbel_model$simulate(
        n = n_sim,
        margin_param = c(mu = 0, kappa = k_val, sigma = 1, xi = true_xi),
        copula_param = 2
    )$x
    order_stats <- sort(x, decreasing = TRUE)
    
    # b. Estimate Rho Star ONCE per sample (de Haan 2016 logic)
    search_upper <- floor(min(n_sim - 1, (2 * n_sim) / log(log(n_sim))))
    search_range <- seq(search_upper, 20, length.out = 30)
    rho_star <- -1 # Default fallback
    
    for (kv in search_range) {
      top_k <- order_stats[1:kv]
      thresh <- order_stats[kv+1]
      lx <- log(top_k) - log(thresh)
      m_vec <- c(M1=mean(lx), M2=mean(lx^2), M3=mean(lx^3), M4=mean(lx^4))
      
      num <- (m_vec[4] - 24 * m_vec[1]^4) * (m_vec[2] - 2 * m_vec[1]^2)
      den <- m_vec[3] - 6 * m_vec[1]^3
      Sk2 <- (3 / 4) * (num / den^2)
      
      if (!is.na(Sk2) && Sk2 > 2/3 && Sk2 < 3/4) {
        rho_star <- (-4 + 6 * Sk2 + sqrt(3 * Sk2 - 2)) / (4 * Sk2 - 3)
        break
      }
    }
    
    # c. Calculate Hill and Hill_BC for every Beta level
    term_rho <- rho_star / (1 - rho_star)
    
    for (b_idx in seq_along(beta_sequence)) {
      beta <- beta_sequence[b_idx]
      k <- floor(n_sim^beta)
      
      lx_target <- log(order_stats[1:k]) - log(order_stats[k+1])
      m1 <- mean(lx_target)
      m2 <- mean(lx_target^2)
      
      # Correction logic
      bias_corr <- if (m1 > 0 && term_rho != 0) {
        (m2 - 2 * m1^2) / (2 * m1 * term_rho)
      } else { 0 }
      
      hill_mc_matrix[b_idx, m] <- m1
      bc_mc_matrix[b_idx, m] <- m1 - bias_corr
    }
  }
  
  # d. Compute Averages for this Kappa
  results_kappa <- data.frame(
    beta = beta_sequence,
    kappa = as.factor(k_val),
    Mean_Hill = rowMeans(hill_mc_matrix, na.rm = TRUE),
    Mean_Hill_BC = rowMeans(bc_mc_matrix, na.rm = TRUE)
  )
  
  all_results[[counter]] <- results_kappa
  counter <- counter + 1
}

# 3. Combine and Reshape for Plotting
final_df <- dplyr::bind_rows(all_results) %>%
  tidyr::pivot_longer(cols = starts_with("Mean"), 
                      names_to = "Estimator", 
                      values_to = "Estimate")


# 4. Plot results
ggplot(final_df, aes(x = beta, y = Estimate, color = Estimator)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = true_xi, color = "black", linetype = "dashed") +
  facet_wrap(~kappa, labeller = label_both) +
  scale_color_manual(values = c("Mean_Hill" = "#377eb8", "Mean_Hill_BC" = "#e41a1c")) +
  labs(title = "Average Bias Trajectories: Hill vs. Hill_BC",
       subtitle = sprintf("MC Iterations: %d | n = %d | xi = %.1f", mc_iterations, n_sim, true_xi),
       x = expression(beta == log(k)/log(n)),
       y = expression(E[hat(xi)])) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_minimal()
