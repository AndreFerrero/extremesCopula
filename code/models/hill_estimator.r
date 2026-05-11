# --- Hill Estimator and Bias-Corrected Hill ---
#' @param data A positive numeric vector of data points.
#' @param k_target The number of top order statistics to use

hill_bc_hat <- function(data, k_target) {
  # data is a positive vector

  n <- length(data)
  order_stats <- sort(data, decreasing = TRUE)

  calc_all_moments <- function(sorted_data, k) {
    top_k <- sorted_data[1:k]
    threshold <- sorted_data[k+1]
    log_excess <- log(top_k) - log(threshold)

    c(
      M1 = mean(log_excess),
      M2 = mean(log_excess^2),
      M3 = mean(log_excess^3),
      M4 = mean(log_excess^4)
    )
  }

  estimate_rho <- function(moments) {
    num <- (moments["M4"] - 24 * moments["M1"]^4) * (moments["M2"] - 2 * moments["M1"]^2)
    den <- moments["M3"] - 6 * moments["M1"]^3
    Sk2 <- (3 / 4) * (num / den^2)

    if (is.na(Sk2) || Sk2 < 2 / 3 || Sk2 > 3 / 4) {
      return(NA)
    }

    (-4 + 6 * Sk2 + sqrt(3 * Sk2 - 2)) / (4 * Sk2 - 3)
  }

  # 1. Estimate Rho using your Supremum logic
  search_upper <- floor(min(n - 1, (2 * n) / log(log(n))))
  search_range <- seq(search_upper, 1, by = -1)

  rho_star <- NA
  for (k_val in search_range) {
    M_k <- calc_all_moments(order_stats, k_val)

    r <- estimate_rho(M_k) # Uses your logic
    if (!is.na(r)) {
      rho_star <- r
      break
    }
  }

  if (is.na(rho_star)) {
    return(NA)
  }

  print(paste("Estimated Rho:", round(rho_star, 4), "at k =", k_val)
  )
  M <- calc_all_moments(order_stats, k_target) # Moments at k_target
  
  # 2. Apply correction at the specific k_target
  gamma_h <- M["M1"] # Standard Hill estimator at k_target

  # Moment 2 calculation
  log_excess <- log(order_stats[1:k_target]) - log(order_stats[k_target + 1])
  m2 <- mean(log_excess^2)

  term_rho <- rho_star / (1 - rho_star)
  bias_corr <- (m2 - 2 * gamma_h^2) / (2 * gamma_h * term_rho)

  xi_hat <- c(gamma_h, gamma_h - bias_corr)
  names(xi_hat) <- c("Hill", "Hill_BC")
  return(xi_hat)
}