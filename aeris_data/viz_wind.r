# =========================
# LOAD DATA
# =========================

source("code/handy_funs.r")
load("aeris_data/mtp_aereoport/wind_hourly_data.RData")

library(dplyr)

# =========================
# HELPER FUNCTIONS
# =========================

make_season_list <- function(winter, spring, summer, autumn) {
  list(
    winter = winter,
    spring = spring,
    summer = summer,
    autumn = autumn
  )
}

plot_hist_seasons <- function(w, sp, su, a, main_title = "") {
  par(mfrow = c(2, 2))
  
  hist(w, breaks = 50, main = paste(main_title, "Winter"), xlab = "Wind Speed")
  hist(sp, breaks = 50, main = paste(main_title, "Spring"), xlab = "Wind Speed")
  hist(su, breaks = 50, main = paste(main_title, "Summer"), xlab = "Wind Speed")
  hist(a, breaks = 50, main = paste(main_title, "Autumn"), xlab = "Wind Speed")
  
  par(mfrow = c(1, 1))
}

plot_lags <- function(w, sp, su, a, title_prefix = "") {
  par(mfrow = c(2, 2))
  lag_plot(w, main = paste(title_prefix, "Winter"))
  lag_plot(sp, main = paste(title_prefix, "Spring"))
  lag_plot(su, main = paste(title_prefix, "Summer"))
  lag_plot(a, main = paste(title_prefix, "Autumn"))
  par(mfrow = c(1, 1))
}

plot_acf <- function(w, sp, su, a, title_prefix = "") {
  par(mfrow = c(2, 2))
  acf(w, main = paste(title_prefix, "Winter"))
  acf(sp, main = paste(title_prefix, "Spring"))
  acf(su, main = paste(title_prefix, "Summer"))
  acf(a, main = paste(title_prefix, "Autumn"))
  par(mfrow = c(1, 1))
}

make_copula_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag)
  df <- df[!is.na(df$x) & !is.na(df$x_lag), ]

  n <- nrow(df)

  df$U_t   <- rank(df$x) / (n + 1)
  df$U_t_1 <- rank(df$x_lag) / (n + 1)

  return(df)
}

plot_copula_base <- function(copula_list, main_title = "") {
  par(mfrow = c(2, 2))
  
  for (season in names(copula_list)) {
    df <- copula_list[[season]]
    
    plot(
      df$U_t_1, df$U_t,
      pch = 16,
      cex = 0.5,
      main = paste(main_title, "-", season),
      xlab = expression(U[t-1]),
      ylab = expression(U[t])
    )
  }
  
  par(mfrow = c(1, 1))
}

plot_upper_tail <- function(copula_list, threshold = 0.9) {
  par(mfrow = c(2, 2))
  
  for (season in names(copula_list)) {
    df <- copula_list[[season]]
    
    df_tail <- df[df$U_t > threshold & df$U_t_1 > threshold, ]
    
    plot(
      df_tail$U_t_1, df_tail$U_t,
      pch = 16,
      cex = 0.5,
      xlim = c(threshold, 1),
      ylim = c(threshold, 1),
      main = paste("Upper Tail -", season),
      xlab = expression(U[t-1]),
      ylab = expression(U[t])
    )
  }
  
  par(mfrow = c(1, 1))
}

# =========================
# TIME SERIES
# =========================

plot(
  wind_full,
  type = "l",
  main = "Wind Speed",
  xlab = "Time Index",
  ylab = "Wind Speed"
)

# =========================
# HISTOGRAMS (ALL DATA)
# =========================

plot_hist_seasons(
  wind_winter, wind_spring, wind_summer, wind_autumn,
  "Histogram -"
)

# =========================
# HISTOGRAMS (POSITIVE ONLY)
# =========================

plot_hist_seasons(
  wind_pos_winter, wind_pos_spring, wind_pos_summer, wind_pos_autumn,
  "Histogram (Positive) -"
)

# =========================
# LAG PLOTS
# =========================

plot_lags(
  wind_winter, wind_spring, wind_summer, wind_autumn,
  "Lag -"
)

plot_lags(
  wind_pos_winter, wind_pos_spring, wind_pos_summer, wind_pos_autumn,
  "Lag (Positive) -"
)

# =========================
# ACF PLOTS
# =========================

plot_acf(
  wind_winter, wind_spring, wind_summer, wind_autumn,
  "ACF -"
)

plot_acf(
  wind_pos_winter, wind_pos_spring, wind_pos_summer, wind_pos_autumn,
  "ACF (Positive) -"
)

# =========================
# COPULA ANALYSIS
# =========================

copula_pos_list <- list(
  winter = make_copula_df(wind_pos_winter),
  spring = make_copula_df(wind_pos_spring),
  summer = make_copula_df(wind_pos_summer),
  autumn = make_copula_df(wind_pos_autumn)
)

plot_copula_base(
  copula_pos_list,
  "Lagged Copula of Positive Wind Speed"
)

# =========================
# UPPER TAIL ANALYSIS
# =========================

plot_upper_tail(
  copula_pos_list,
  threshold = 0.90
)
