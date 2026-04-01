source("code/handy_funs.r")
library(dplyr)

# =========================
# 1. LOAD & PREP DATA
# =========================
data <- read.csv("copernicus_data/cop_data.csv")

data <- data %>%
  filter(!is.na(fg10)) %>%
  mutate(valid_time = as.POSIXct(valid_time, tz = "UTC"))

# Keep last 10 years
cutoff_date <- as.POSIXct("2025-01-01 00:00:00", tz = "UTC")

data <- data %>%
  filter(valid_time >= cutoff_date)

# Add time features
data <- data %>%
  mutate(
    date = as.Date(valid_time),
    month = as.numeric(format(valid_time, "%m")),
    season = case_when(
      month %in% c(12,1,2) ~ "Winter",
      month %in% 3:5 ~ "Spring",
      month %in% 6:8 ~ "Summer",
      TRUE ~ "Autumn"
    )
  )

# =========================
# 2. DAILY MAX
# =========================
daily_max <- data %>%
  group_by(date) %>%
  summarise(fg10 = max(fg10), .groups = "drop") %>%
  mutate(
    month = as.numeric(format(date, "%m")),
    season = case_when(
      month %in% c(12,1,2) ~ "Winter",
      month %in% 3:5 ~ "Spring",
      month %in% 6:8 ~ "Summer",
      TRUE ~ "Autumn"
    )
  )

