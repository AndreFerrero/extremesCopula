# Asymptotically Unbiased Estimator for EVI (de Haan et al., 2016)
# Grounded Version for Dependent Data (e.g., Copula Markov Processes)

#' Calculate Hill Estimator and Moments
calc_moments <- function(data, k, alpha) {
  n <- length(data)
  sorted_data <- sort(data)
  top_k <- sorted_data[(n-k+1):n]
  threshold <- sorted_data[n-k]
  log_excess <- log(top_k) - log(threshold)
  return(mean(log_excess^alpha))
}

#' Estimate the second-order parameter rho (alpha = 2 case)
#' Now "quiet" (no printing) to allow for efficient looping
estimate_rho <- function(data, k) {
  M1 <- calc_moments(data, k, 1)
  M2 <- calc_moments(data, k, 2)
  M3 <- calc_moments(data, k, 3)
  M4 <- calc_moments(data, k, 4)
  
  num <- (M4 - 24*M1^4) * (M2 - 2*M1^2)
  den <- M3 - 6*M1^3
  Sk2 <- (3/4) * (num / den^2)
  
  # Existence condition from Section 6.2
  if (is.na(Sk2) || Sk2 < 2/3 || Sk2 > 3/4) return(NA)
  
  rho_hat <- (-4 + 6*Sk2 + sqrt(3*Sk2 - 2)) / (4*Sk2 - 3)
  return(rho_hat)
}

#' Plot EVI Stability using the Supremum Rho (de Haan et al., 2016, Section 6.2)
#' @param data Numeric vector of observations
#' @param k_range Vector of kn values to plot (the lower-order sequence)
plot_unbiased_evi_supremum <- function(data, k_range) {
  n <- length(data)
  m <- length(data[data > 0])
  
  # --- 1. Find k_rho as the Supremum (Page 335, Section 6.2) ---
  # "The parameter rho is estimated as rho_hat_krho with 
  # krho = sup{ k : k <= min(m-1, 2m/loglogm) and rho_hat exists }"
  
  search_upper_bound <- floor(min(m - 1, (2 * m) / log(log(m))))
  # Search from the upper bound downwards to find the SUPREMUM
  search_range <- search_upper_bound:10 
  
  best_krho <- NA
  rho_star <- NA
  
  for (k_val in search_range) {
    r <- estimate_rho(data, k_val)
    if (!is.na(r)) {
      rho_star <- r
      best_krho <- k_val
      break # Stop at the first (highest) valid k
    }
  }
  
  if (is.na(rho_star)) {
    stop("Could not find any k where rho_hat exists.")
  }
  
  cat("Supremum k_rho found at:", best_krho, "\n")
  cat("Fixed Rho estimate for correction:", round(rho_star, 4), "\n")
  
  # --- 2. Calculate EVI across k_range using the FIXED rho_star ---
  results <- data.frame(
    k = k_range,
    hill = sapply(k_range, function(k) calc_moments(data, k, 1)),
    unbiased = NA
  )
  
  # Constant term for bias correction
  term_rho <- rho_star / (1 - rho_star)
  
  for (i in 1:nrow(results)) {
    k_val <- results$k[i]
    # We only apply correction where kn < krho (as per theory on page 326)
    # though in practice we can plot further to see where it breaks.
    
    gamma_h <- results$hill[i]
    m2 <- calc_moments(data, k_val, 2)
    
    # Equation 4.2 using the FIXED rho_star
    bias_corr <- (m2 - 2 * gamma_h^2) / (2 * gamma_h * term_rho)
    results$unbiased[i] <- gamma_h - bias_corr
  }
  
  # --- 3. Plotting ---
  ymin <- min(c(results$hill, results$unbiased), na.rm = TRUE)
  ymax <- max(c(results$hill, results$unbiased), na.rm = TRUE)
  
  # Visual limit to keep plot readable
  y_lims <- c(max(-0.5, ymin), min(1.2, ymax))
  
  plot(results$k, results$hill, type = "l", col = "red", lwd = 2, lty = 2,
       ylim = y_lims,
       xlab = "kn (Number of top observations)", 
       ylab = expression(hat(gamma)),
       main = "EVI Stability")
  
  lines(results$k, results$unbiased, col = "blue", lwd = 2)
  
  # Mark the k_rho used for the correction
  abline(v = best_krho, col = "darkgreen", lty = 3)
  text(best_krho, y_lims[1] + 0.1, "k_rho", col = "darkgreen", pos = 2)
  
  grid()
  legend("topright", 
         legend = c("Hill (Biased)", paste0("Unbiased (fixed rho=", round(rho_star,2), ")")),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n")
  
  return(invisible(results))
}
