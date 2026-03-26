# =========================
# LOAD + PREPARE DATA (WIND)
# =========================

library(ncdf4)
library(dplyr)
library(lubridate)

# Folder
nc_folder <- "aeris_data/mtp_aereoport/"
files <- paste0(nc_folder, "34154001_MONTPELLIER-AEROPORT_MTO_1H_", 2009:2025, ".nc")

# Store data
wind_list <- list()

for (f in files) {
    nc_data <- nc_open(f)

    time <- ncvar_get(nc_data, "time")
    wind <- ncvar_get(nc_data, "ws")

    wind[is.nan(wind)] <- NA

    time_posix <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")

    df_year <- data.frame(
        time = time_posix,
        wind_speed = wind
    )

    wind_list[[f]] <- df_year

    nc_close(nc_data)
}

# Combine all
df <- bind_rows(wind_list) %>%
    filter(!is.na(wind_speed)) %>%
    arrange(time)

# Add season
df <- df %>%
    mutate(
        month = month(time),
        season = case_when(
            month %in% c(12, 1, 2) ~ "winter",
            month %in% c(3, 4, 5) ~ "spring",
            month %in% c(6, 7, 8) ~ "summer",
            month %in% c(9, 10, 11) ~ "autumn"
        )
    )

df_12h <- df %>%
    mutate(
        # Floor time to nearest 6-hour block
        time_12h = floor_date(time, unit = "12 hours")
    ) %>%
    group_by(time_12h) %>%
    summarise(
        wind_max_12h = max(wind_speed, na.rm = TRUE)
    ) %>%
    ungroup()

# Add season to 12h maxima
df_12h <- df_12h %>%
    mutate(
        month = month(time_12h),
        season = case_when(
            month %in% c(12, 1, 2) ~ "winter",
            month %in% c(3, 4, 5) ~ "spring",
            month %in% c(6, 7, 8) ~ "summer",
            month %in% c(9, 10, 11) ~ "autumn"
        )
    )

df_daily <- df %>%
    mutate(
        date = as.Date(time) # floor to day
    ) %>%
    group_by(date) %>%
    summarise(
        wind_max_daily = max(wind_speed, na.rm = TRUE)
    ) %>%
    ungroup()

# Add season to daily maxima
df_daily <- df_daily %>%
    mutate(
        month = month(date),
        season = case_when(
            month %in% c(12, 1, 2) ~ "winter",
            month %in% c(3, 4, 5) ~ "spring",
            month %in% c(6, 7, 8) ~ "summer",
            month %in% c(9, 10, 11) ~ "autumn"
        )
    )

# =========================
# OBJECTS
# =========================

# 1. Full vector (wind, includes zeros)
wind_full <- df$wind_speed

# 2. 4 seasonal vectors
wind_winter <- df %>%
    filter(season == "winter") %>%
    pull(wind_speed)
wind_spring <- df %>%
    filter(season == "spring") %>%
    pull(wind_speed)
wind_summer <- df %>%
    filter(season == "summer") %>%
    pull(wind_speed)
wind_autumn <- df %>%
    filter(season == "autumn") %>%
    pull(wind_speed)

# Keep time
time_full <- df$time

xl <- 0
wind_pos_full <- df %>%
    filter(wind_speed > xl) %>%
    pull(wind_speed)

# 4. 4 seasonal vectors (positive only)
wind_pos_winter <- df %>%
    filter(season == "winter", wind_speed > xl) %>%
    pull(wind_speed)
wind_pos_spring <- df %>%
    filter(season == "spring", wind_speed > xl) %>%
    pull(wind_speed)
wind_pos_summer <- df %>%
    filter(season == "summer", wind_speed > xl) %>%
    pull(wind_speed)
wind_pos_autumn <- df %>%
    filter(season == "autumn", wind_speed > xl) %>%
    pull(wind_speed)

# =========================
# SAVE EVERYTHING
# =========================

save(
    wind_full,
    wind_winter, wind_spring, wind_summer, wind_autumn,
    time_full,
    wind_pos_full,
    wind_pos_winter, wind_pos_spring, wind_pos_summer, wind_pos_autumn,
    file = "aeris_data/mtp_aereoport/wind_hourly_data.RData"
)
