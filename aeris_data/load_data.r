# =========================
# LOAD + PREPARE DATA (WIND)
# =========================

library(ncdf4)
library(dplyr)
library(lubridate)

nc_folder <- "aeris_data/mtp_aereoport/data"
files <- paste0(nc_folder, "34154001_MONTPELLIER-AEROPORT_MTO_1H_", 2009:2025, ".nc")

# ------------------------------------------------------------------
# 1. Read raw hourly data
# ------------------------------------------------------------------

wind_list <- lapply(files, function(f) {
  nc_data <- nc_open(f)
  time  <- ncvar_get(nc_data, "time")
  wind  <- ncvar_get(nc_data, "ws")
  nc_close(nc_data)

  wind[is.nan(wind)] <- NA
  data.frame(
    time       = as.POSIXct(time, origin = "1970-01-01", tz = "UTC"),
    wind_speed = wind
  )
})

df <- bind_rows(wind_list) %>%
  filter(!is.na(wind_speed)) %>%
  arrange(time) %>%
  mutate(
    month  = month(time),
    season = case_when(
      month %in% c(12, 1, 2) ~ "winter",
      month %in% c(3, 4, 5)  ~ "spring",
      month %in% c(6, 7, 8)  ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn"
    )
  )

# ------------------------------------------------------------------
# 2. Helper: split a numeric column into full + seasonal vectors,
#    each paired with its own time index so time-series plots always
#    have matching lengths regardless of aggregation level.
#
#    time_col : name of the timestamp column in df (default "time")
#    col      : name of the value column
#    threshold: if provided, keep only rows strictly > threshold
# ------------------------------------------------------------------

split_by_season <- function(df, col, time_col = NULL, threshold = NULL) {
  d <- df
  if (!is.null(threshold)) d <- filter(d, .data[[col]] > threshold)

  # Detect the time column automatically if not supplied
  if (is.null(time_col)) {
    candidates <- c("time", "date", "time_12h")
    time_col   <- candidates[candidates %in% names(d)][1]
  }

  make_entry <- function(rows) {
    list(
      values = rows[[col]],
      time   = rows[[time_col]]
    )
  }

  list(
    full   = make_entry(d),
    winter = make_entry(d[d$season == "winter", ]),
    spring = make_entry(d[d$season == "spring", ]),
    summer = make_entry(d[d$season == "summer", ]),
    autumn = make_entry(d[d$season == "autumn", ])
  )
}

# ------------------------------------------------------------------
# 3. Build aggregated data frames
# ------------------------------------------------------------------

df_daily <- df %>%
  mutate(date = as.Date(time)) %>%
  group_by(date) %>%
  summarise(wind_max_daily = max(wind_speed, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    month  = month(date),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "winter",
      month %in% c(3, 4, 5)   ~ "spring",
      month %in% c(6, 7, 8)   ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn"
    )
  )

df_12h <- df %>%
  mutate(time_12h = floor_date(time, unit = "12 hours")) %>%
  group_by(time_12h) %>%
  summarise(wind_max_12h = max(wind_speed, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    month  = month(time_12h),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "winter",
      month %in% c(3, 4, 5)   ~ "spring",
      month %in% c(6, 7, 8)   ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn"
    )
  )

# ------------------------------------------------------------------
# 4. Assemble the central wind_data list
#
#    Structure:
#      wind_data[[dataset_key]]$full
#      wind_data[[dataset_key]]$winter  /  $spring  /  $summer  /  $autumn
#
#    To add a new dataset: append one more entry below.
# ------------------------------------------------------------------

wind_data <- list(
  hourly     = split_by_season(df,       "wind_speed"),
  hourly_pos = split_by_season(df,       "wind_speed", threshold = 0),
  daily_max  = split_by_season(df_daily, "wind_max_daily"),
  daily_max_pos = split_by_season(df_daily, "wind_max_daily", threshold = 0),
  max_12h    = split_by_season(df_12h,   "wind_max_12h"),
  max_12h_pos = split_by_season(df_12h,  "wind_max_12h", threshold = 0)
)

# ------------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------------

save(wind_data, file = "aeris_data/mtp_aereoport/wind_data.RData")
message("Saved wind_data with keys: ", paste(names(wind_data), collapse = ", "))