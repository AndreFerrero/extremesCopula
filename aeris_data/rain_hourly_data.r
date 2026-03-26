# =========================
# LOAD + PREPARE DATA
# =========================

library(ncdf4)
library(dplyr)
library(lubridate)

# Folder
nc_folder <- "aeris_data/mtp_aereoport/"
files <- paste0(nc_folder, "34154001_MONTPELLIER-AEROPORT_MTO_1H_", 2006:2025, ".nc")

# Store data
precip_list <- list()

for (f in files) {
  nc_data <- nc_open(f)

  time <- ncvar_get(nc_data, "time")
  precip <- ncvar_get(nc_data, "cumul_precip")

  precip[is.nan(precip)] <- NA

  time_posix <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")

  df_year <- data.frame(
    time = time_posix,
    precipitation_mm = precip
  )

  precip_list[[f]] <- df_year

  nc_close(nc_data)
}

# Combine all
df <- bind_rows(precip_list) %>%
  filter(!is.na(precipitation_mm)) %>%
  arrange(time)

# Add season
df <- df %>%
  mutate(
    month = month(time),
    season = case_when(
      month %in% c(12,1,2)  ~ "winter",
      month %in% c(3,4,5)   ~ "spring",
      month %in% c(6,7,8)   ~ "summer",
      month %in% c(9,10,11) ~ "autumn"
    )
  )

# =========================
# OBJECTS
# =========================

# 1. Full vector (with zeros)
rain_full <- df$precipitation_mm

# 2. 4 seasonal vectors (with zeros)
rain_winter <- df %>% filter(season == "winter") %>% pull(precipitation_mm)
rain_spring <- df %>% filter(season == "spring") %>% pull(precipitation_mm)
rain_summer <- df %>% filter(season == "summer") %>% pull(precipitation_mm)
rain_autumn <- df %>% filter(season == "autumn") %>% pull(precipitation_mm)

# 3. Full vector (positive only)
xl <- 0
rain_pos_full <- df %>%
  filter(precipitation_mm > xl) %>%
  pull(precipitation_mm)

# 4. 4 seasonal vectors (positive only)
rain_pos_winter <- df %>% filter(season == "winter", precipitation_mm > xl) %>% pull(precipitation_mm)
rain_pos_spring <- df %>% filter(season == "spring", precipitation_mm > xl) %>% pull(precipitation_mm)
rain_pos_summer <- df %>% filter(season == "summer", precipitation_mm > xl) %>% pull(precipitation_mm)
rain_pos_autumn <- df %>% filter(season == "autumn", precipitation_mm > xl) %>% pull(precipitation_mm)

# Keep time if needed later
time_full <- df$time

# =========================
# SAVE EVERYTHING
# =========================

save(
  rain_full,
  rain_winter, rain_spring, rain_summer, rain_autumn,
  rain_pos_full,
  rain_pos_winter, rain_pos_spring, rain_pos_summer, rain_pos_autumn,
  time_full,
  file = "aeris_data/mtp_aereoport/rain_hourly_data.RData"
)