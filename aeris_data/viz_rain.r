# =========================
# LOAD DATA
# =========================
source("code/handy_funs.r")
load("aeris_data/mtp_aereoport/rain_hourly_data.RData")

library(ggplot2)
library(dplyr)

# =========================
# Time series
# =========================

plot(rain_full, type = "l", main = "Rainfall (mm)", xlab = "Time Index", ylab = "Precipitation (mm)")

df_hist <- data.frame(
  value = c(rain_winter, rain_spring, rain_summer, rain_autumn),
  season = rep(c("winter", "spring", "summer", "autumn"),
    times = c(
      length(rain_winter),
      length(rain_spring),
      length(rain_summer),
      length(rain_autumn)
    )
  )
)

ggplot(df_hist, aes(x = value)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~season, scales = "free_y") +
  labs(
    title = "Histogram of Rainfall by Season",
    x = "Precipitation (mm)",
    y = "Count"
  ) +
  theme_minimal()

par(mfrow = c(2, 2))
lag_plot(rain_winter)
lag_plot(rain_spring)
lag_plot(rain_summer)
lag_plot(rain_autumn)
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
acf(rain_winter, main = "ACF of Rainfall (Winter)")
acf(rain_spring, main = "ACF of Rainfall (Spring)")
acf(rain_summer, main = "ACF of Rainfall (Summer)")
acf(rain_autumn, main = "ACF of Rainfall (Autumn)")
par(mfrow = c(1, 1))

# =========================
# HISTOGRAMS (POSITIVE RAIN)
# =========================

df_hist_pos <- data.frame(
  value = c(rain_pos_winter, rain_pos_spring, rain_pos_summer, rain_pos_autumn),
  season = rep(c("winter", "spring", "summer", "autumn"),
    times = c(
      length(rain_pos_winter),
      length(rain_pos_spring),
      length(rain_pos_summer),
      length(rain_pos_autumn)
    )
  )
)

ggplot(df_hist_pos, aes(x = value)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~season, scales = "free_y") +
  labs(
    title = "Histogram of Positive Rainfall by Season",
    x = "Precipitation (mm)",
    y = "Count"
  ) +
  theme_minimal()

par(mfrow = c(2, 2))
lag_plot(rain_pos_winter)
lag_plot(rain_pos_spring)
lag_plot(rain_pos_summer)
lag_plot(rain_pos_autumn)
par(mfrow = c(1, 1))

# =========================
# FUNCTION: COPULA DATA
# =========================

make_copula_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) %>%
    filter(!is.na(x) & !is.na(x_lag))

  df$U_t <- rank(df$x) / (nrow(df) + 1)
  df$U_t_1 <- rank(df$x_lag) / (nrow(df) + 1)

  return(df)
}

# =========================
# COPULA PER SEASON
# =========================
copula_list <- list(
  winter = make_copula_df(rain_winter),
  spring = make_copula_df(rain_spring),
  summer = make_copula_df(rain_summer),
  autumn = make_copula_df(rain_autumn)
)

df_copula <- bind_rows(copula_list, .id = "season")

ggplot(df_copula, aes(x = U_t_1, y = U_t)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~season) +
  labs(
    title = "Lagged Copula of Rainfall Intensity by Season",
    x = expression(U[t - 1]),
    y = expression(U[t])
  ) +
  theme_bw()

copula_pos_list <- list(
  winter = make_copula_df(rain_pos_winter),
  spring = make_copula_df(rain_pos_spring),
  summer = make_copula_df(rain_pos_summer),
  autumn = make_copula_df(rain_pos_autumn)
)

df_copula_pos <- bind_rows(copula_list, .id = "season")

# 1. Filter the data first to only include the "Heavy Rain" region (U > 0.8)
df_zoom <- df_copula_pos[df_copula_pos$U_t_1 > 0.98 & df_copula_pos$U_t > 0.98, ]

# 2. Plot the subsetted data
ggplot(df_zoom, aes(x = U_t_1, y = U_t)) +
  # Add jitter slightly if points are still sitting on top of each other
  # because of quantization, otherwise use geom_point
  geom_point(alpha = 0.2, size = 0.5) +
  facet_wrap(~season) +
  # Use xlim/ylim here since the data is already filtered
  scale_x_continuous(limits = c(0.98, 1.0)) +
  scale_y_continuous(limits = c(0.98, 1.0)) +
  labs(
    title = "Zoomed Upper Tail (U > 0.98)",
    subtitle = "Positive Rainfall Intensity Transitions",
    x = expression(U[t - 1]),
    y = expression(U[t])
  ) +
  theme_bw()


# =========================
# ACF OF POSITIVE RAINFALL
# =========================
par(mfrow = c(2, 2))
acf(rain_pos_winter, main = "ACF of Positive Rainfall (Winter)")
acf(rain_pos_spring, main = "ACF of Positive Rainfall (Spring)")
acf(rain_pos_summer, main = "ACF of Positive Rainfall (Summer)")
acf(rain_pos_autumn, main = "ACF of Positive Rainfall (Autumn)")
par(mfrow = c(1, 1))
