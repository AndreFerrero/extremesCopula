source("ERA5/load_data.R")

winter_hourly_gust <- data$fg10[data$season == "Winter"]

threshold <- quantile(winter_hourly_gust, 0.9, na.rm = TRUE)

library(gganimate)

df <- data.frame(
  time = 1:length(winter_hourly_gust),
  gust = winter_hourly_gust,
  exceed = winter_hourly_gust > threshold
)

df$gust_above <- ifelse(df$gust > threshold, df$gust, NA)

p <- ggplot(df, aes(time, gust)) +
  geom_line(color = "grey80") +  # full series
  geom_line(aes(y = gust_above), color = "red", linewidth = 1) +  # ONLY above threshold
  geom_hline(yintercept = threshold, color = "red", linewidth = 1) +
  transition_reveal(time) +
  labs(title = "Time: {frame_along}") +
  theme_minimal()

animate(p, nframes = 200, fps = 20)
anim_save("ERA5/winter_gust_animation.gif")
