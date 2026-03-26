# =========================
# LOAD DATA
# =========================

source("code/handy_funs.r")
load("aeris_data/mtp_aereoport/wind_hourly_data.RData")

library(ggplot2)
library(dplyr)

# =========================
# HELPER FUNCTIONS
# =========================

make_season_df <- function(winter, spring, summer, autumn) {
  data.frame(
    value = c(winter, spring, summer, autumn),
    season = rep(
      c("winter", "spring", "summer", "autumn"),
      times = c(length(winter), length(spring), length(summer), length(autumn))
    )
  )
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

  df <- data.frame(x = x, x_lag = x_lag) %>%
    filter(!is.na(x) & !is.na(x_lag))

  df$U_t   <- rank(df$x) / (nrow(df) + 1)
  df$U_t_1 <- rank(df$x_lag) / (nrow(df) + 1)

  return(df)
}

make_copula_plot <- function(data_list, title) {
  df <- bind_rows(data_list, .id = "season")

  ggplot(df, aes(x = U_t_1, y = U_t)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~season) +
    labs(
      title = title,
      x = expression(U[t - 1]),
      y = expression(U[t])
    ) +
    theme_bw()
}

# =========================
# TIME SERIES
# =========================

plot(wind_full,
     type = "l",
     main = "Wind Speed",
     xlab = "Time Index",
     ylab = "Wind Speed")

# =========================
# HISTOGRAMS (ALL DATA)
# =========================

df_hist <- make_season_df(
  wind_winter, wind_spring, wind_summer, wind_autumn
)

ggplot(df_hist, aes(x = value)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~season, scales = "free_y") +
  labs(
    title = "Histogram of Wind Speed by Season",
    x = "Wind Speed",
    y = "Count"
  ) +
  theme_bw()

# =========================
# LAG PLOTS
# =========================

plot_lags(
  wind_winter, wind_spring, wind_summer, wind_autumn,
  "Lag -"
)

plot_lags(
  wind_pos_winter, wind_pos_spring, wind_pos_summer, wind_pos_autumn,
  "Lag -"
)

# =========================
# ACF PLOTS
# =========================

plot_acf(
  wind_winter, wind_spring, wind_summer, wind_autumn,
  "ACF -"
)

plot_acf(
  wind_pos_winter, wind_spring, wind_summer, wind_autumn,
  "ACF -"
)

# =========================
# COPULA ANALYSIS
# =========================

copula_list <- list(
  winter = make_copula_df(wind_winter),
  spring = make_copula_df(wind_spring),
  summer = make_copula_df(wind_summer),
  autumn = make_copula_df(wind_autumn)
)

df_copula <- bind_rows(copula_list, .id = "season")

make_copula_plot(
  copula_list,
  "Lagged Copula of Wind Speed"
)

# =========================
# UPPER TAIL ANALYSIS
# =========================

df_zoom <- df_copula %>%
  filter(U_t > 0.9 & U_t_1 > 0.9)

ggplot(df_zoom, aes(x = U_t_1, y = U_t)) +
  geom_point(alpha = 0.2, size = 0.5) +
  facet_wrap(~season) +
  coord_cartesian(xlim = c(0.9, 1), ylim = c(0.9, 1)) +
  labs(
    title = "Upper Tail Dependence (Wind Speed)",
    x = expression(U[t - 1]),
    y = expression(U[t])
  ) +
  theme_bw()
