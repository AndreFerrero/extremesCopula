source("ERA5/load_data.r")
source("megpd/megpd_fit.r")

winter_hourly_gust <- data$fg10[data$season == "Winter"]

lag_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) |>
    na.omit()

  return(df)
}

winter_lag <- lag_df(winter_hourly_gust)
winter_R <- winter_lag$x + winter_lag$x_lag
winter_V <- log(winter_lag$x/winter_lag$x_lag)

# megpd_markov_fit_m5 <- fit_megpd_bivariate(winter_lag[,1], winter_lag[,2], m = 5)
megpd_markov_fit_m10 <- fit_megpd_bivariate(winter_lag[,1], winter_lag[,2], m = 10)

# summary(megpd_markov_fit_m5$radial_fit)
summary(megpd_markov_fit_m10$radial_fit)

plot(megpd_markov_fit_m10$angular_fit)

plot_delta(megpd_markov_fit_m10$angular_fit, winter_R)

plot(log(winter_lag$x), log(winter_lag$x_lag),
     pch = 16, col = rgb(0,0,0,0.3))

plot(winter_R, winter_V,
     pch = 16, col = rgb(0,0,0,0.3),
     main = "Angular spread vs radial component",
     xlab = "R", ylab = "V")

fitted_chi <- plot_fitted_chi(fit_megpd = megpd_markov_fit_m10, n_mc = 50)

emp_chi <- texmex::chi(winter_lag)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 2], col = "blue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 1], col = "lightblue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 3], col = "lightblue", lwd = 2)
