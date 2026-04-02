source("ERA5/load_data.R")

# =========================
# 3. COPULA FUNCTION (ECDF)
# =========================
get_lag_copula <- function(x, lag = 1) {
  # Use base R for speed and to avoid NA issues with dplyr::lag inside loops
  n_total <- length(x)
  
  # Align vectors based on lag
  x_t     <- x[(lag + 1):n_total]
  x_prev  <- x[1:(n_total - lag)]
  
  n <- length(x_t)
  
  # Probability Integral Transform (PIT) using ECDF ranks
  u_t    <- rank(x_t) / (n + 1)
  u_prev <- rank(x_prev) / (n + 1)
  
  return(data.frame(U_t = u_t, U_lag = u_prev))
}

# =========================
# 4. GENERIC PLOTTING FUNS
# =========================

plot_hist_by_season <- function(df, var, title_prefix) {
  par(mfrow = c(2,2))
  
  for (s in c("Winter","Spring","Summer","Autumn")) {
    hist(df[[var]][df$season == s],
         breaks = 20,
         main = paste(title_prefix, "-", s),
         xlab = var,
         col = "lightblue")
  }
  
  par(mfrow = c(1,1))
}

plot_acf_by_season <- function(df, var, title_prefix) {
  par(mfrow = c(2,2))
  
  for (s in c("Winter","Spring","Summer","Autumn")) {
    acf(df[[var]][df$season == s],
        main = paste("ACF:", title_prefix, "-", s))
  }
  
  par(mfrow = c(1,1))
}

plot_pacf_by_season <- function(df, var, title_prefix) {
  par(mfrow = c(2,2))
  
  for (s in c("Winter","Spring","Summer","Autumn")) {
    pacf(df[[var]][df$season == s],
        main = paste("PACF:", title_prefix, "-", s))
  }
  
  par(mfrow = c(1,1))
}

plot_copula_by_season <- function(df, var, title_prefix, lag = 1, u = 0.9, l = 0.05) {
  seasons <- c("Winter","Spring","Summer","Autumn")
  
  for (s in seasons) {
    
    x <- df[[var]][df$season == s]
    cop <- get_lag_copula(x, lag = lag)
    
    par(mfrow = c(1,3))
    
    # Full copula
    plot(cop$U_lag, cop$U_t,
         pch = 16, cex = 0.4,
         main = paste(title_prefix, "-", s),
         xlab = paste0("U[t-", lag, "]"),
         ylab = paste0("U[t]"))
    
    # Upper tail
    cop_u <- cop[cop$U_t > u & cop$U_lag > u, ]
    
    plot(cop_u$U_lag, cop_u$U_t,
         pch = 16, cex = 0.4,
         main = paste("Upper Tail -", s),
         xlab = paste0("U[t-", lag, "]"),
         ylab = paste0("U[t]"))
    
    # Lower tail
    cop_l <- cop[cop$U_t < l & cop$U_lag < l, ]
    
    plot(cop_l$U_lag, cop_l$U_t,
         pch = 16, cex = 0.4,
         main = paste("Lower Tail -", s),
         xlab = paste0("U[t-", lag, "]"),
         ylab = paste0("U[t]"))
  }
  
  par(mfrow = c(1,1))
}


# =========================
# 5. RUN ANALYSIS
# =========================

# ---- HOURLY ----
plot_hist_by_season(data, "fg10", "Hourly Wind Gust")
plot_acf_by_season(data, "fg10", "Hourly Wind Gust")
plot_pacf_by_season(data, "fg10", "Hourly Wind Gust")
plot_copula_by_season(data, "fg10", lag = 1, "Hourly Copula")
plot_copula_by_season(data, "fg10", lag = 2, "Hourly Lag-2 Analysis")

# ---- DAILY MAX ----
plot_hist_by_season(daily_max, "fg10", "Daily Max Wind Gust")
plot_acf_by_season(daily_max, "fg10", "Daily Max Wind Gust")
plot_pacf_by_season(daily_max, "fg10", "Daily Max Wind Gust")
plot_copula_by_season(daily_max, "fg10", "Daily Max Copula")
