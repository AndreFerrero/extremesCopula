source("ERA5/load_data.R")

# =========================
# 3. COPULA FUNCTION (ECDF)
# =========================
make_copula_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) %>%
    na.omit()

  n <- nrow(df)

  df$U_t   <- rank(df$x) / (n + 1)
  df$U_t_1 <- rank(df$x_lag) / (n + 1)

  return(df)
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

plot_copula_by_season <- function(df, var, title_prefix, u = 0.9, l = 0.05) {
  seasons <- c("Winter","Spring","Summer","Autumn")
  
  for (s in seasons) {
    
    x <- df[[var]][df$season == s]
    cop <- make_copula_df(x)
    
    par(mfrow = c(1,3))
    
    # Full copula
    plot(cop$U_t_1, cop$U_t,
         pch = 16, cex = 0.4,
         main = paste(title_prefix, "-", s),
         xlab = expression(U[t-1]),
         ylab = expression(U[t]))
    
    # Upper tail
    cop_u <- cop[cop$U_t > u & cop$U_t_1 > u, ]
    
    plot(cop_u$U_t_1, cop_u$U_t,
         pch = 16, cex = 0.4,
         main = paste("Upper Tail -", s),
         xlab = expression(U[t-1]),
         ylab = expression(U[t]))
    
    # Lower tail
    cop_l <- cop[cop$U_t < l & cop$U_t_1 < l, ]
    
    plot(cop_l$U_t_1, cop_l$U_t,
         pch = 16, cex = 0.4,
         main = paste("Lower Tail -", s),
         xlab = expression(U[t-1]),
         ylab = expression(U[t]))
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
plot_copula_by_season(data, "fg10", "Hourly Copula")

# ---- DAILY MAX ----
plot_hist_by_season(daily_max, "fg10", "Daily Max Wind Gust")
plot_acf_by_season(daily_max, "fg10", "Daily Max Wind Gust")
plot_pacf_by_season(daily_max, "fg10", "Daily Max Wind Gust")
plot_copula_by_season(daily_max, "fg10", "Daily Max Copula")
