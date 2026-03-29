# =========================
# EDA VISUALISATION — WIND
# =========================
#
# HOW TO USE
# ----------
# 1. Set `config$data_key` to any key in wind_data (see load_wind_data.R).
# 2. Toggle the plots you want in `config$plots`.
# 3. Set `config$save_plots = TRUE` to write PDFs, or FALSE to screen only.
# 4. Source the file. Everything else is automatic.
#
# Available keys (extend in load_wind_data.R):
#   "hourly"        — raw hourly wind speed (all values incl. zeros)
#   "hourly_pos"    — hourly, positive only (> 0)
#   "daily_max"     — daily maxima
#   "daily_max_pos" — daily maxima, positive only
#   "max_12h"       — 12-hour block maxima
#   "max_12h_pos"   — 12-hour block maxima, positive only
# ==========================================================================

source("code/handy_funs.r")
load("aeris_data/mtp_aereoport/wind_data.RData")

library(dplyr)

# ==========================================================================
# CONFIG — only change things here
# ==========================================================================

config <- list(
  data_key = "daily_max_pos",          # which dataset to visualise
  seasons  = c("winter", "spring", "summer", "autumn"),  # subset or reorder
  plots    = c(                    # toggle any subset
    "time_series",
    "histogram",
    "acf",
    "lag",
    "copula",
    "copula_upper_tail"
  ),
  copula_tail_threshold = 0.90,    # quantile threshold for upper-tail plot
  hist_breaks = 50,
  xlab = "Wind speed (m/s)",       # axis label used across all plots

  # -------------------------------------------------------------------
  # Saving
  # Set save_plots = TRUE to write one PDF per plot type.
  #
  # save_dir can be:
  #   - a single string: used for all data_keys
  #       save_dir = "outputs/eda"
  #
  #   - a named list to route each key to its own folder;
  #     use ".default" as a fallback for any unlisted key:
  #       save_dir = list(
  #         daily_max     = "outputs/daily_max",
  #         hourly        = "outputs/hourly",
  #         .default      = "outputs/other"
  #       )
  #
  # Directories are created automatically if they do not exist.
  # Output files are named: <data_key>_<plot_type>.pdf
  # -------------------------------------------------------------------
  save_plots = TRUE,
  save_dir   = "aeris_data/mtp_aereoport/mtp_figures/max_daily"
)

# ==========================================================================
# INTERNALS — no need to touch below this line
# ==========================================================================

# Pull the chosen dataset (a list with $full + four seasonal entries,
# each containing $values and $time)
d <- wind_data[[ config$data_key ]]

if (is.null(d)) {
  stop(
    "data_key '", config$data_key, "' not found in wind_data.\n",
    "Available keys: ", paste(names(wind_data), collapse = ", ")
  )
}

seasons <- config$seasons

# --------------------------------------------------------------------------
# Resolve output directory for the current data_key
# --------------------------------------------------------------------------

resolve_save_dir <- function() {
  sd <- config$save_dir
  if (is.character(sd)) return(sd)                       # single path
  if (is.list(sd)) {
    key <- config$data_key
    if (!is.null(sd[[ key ]])) return(sd[[ key ]])       # exact match
    if (!is.null(sd$.default))  return(sd$.default)      # fallback
    stop("save_dir has no entry for '", key, "' and no '.default'.")
  }
  stop("config$save_dir must be a string or a named list.")
}

# Open a PDF device (creating the directory if needed), run expr, close it.
# Falls back to screen only if save_plots = FALSE.
with_pdf <- function(plot_type, expr) {
  if (!isTRUE(config$save_plots)) {
    force(expr)
    return(invisible(NULL))
  }
  dir  <- resolve_save_dir()
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(config$data_key, "_", plot_type, ".pdf"))
  pdf(path)
  on.exit({ dev.off(); message("Saved: ", path) })
  force(expr)
}

# --------------------------------------------------------------------------
# Helper: 2x2 seasonal panel
# Each plot function receives a list entry with $values and $time
# --------------------------------------------------------------------------

with_season_panels <- function(seasons, fn, ...) {
  par(mfrow = c(2, 2))
  for (s in seasons) fn(d[[s]], season = s, ...)
  par(mfrow = c(1, 1))
}

# --------------------------------------------------------------------------
# Helper: compute empirical copula data frame (lag-1)
# --------------------------------------------------------------------------

make_copula_df <- function(x) {
  x_lag <- dplyr::lag(x)
  df    <- data.frame(x = x, x_lag = x_lag)
  df    <- df[!is.na(df$x) & !is.na(df$x_lag), ]
  n     <- nrow(df)
  df$U_t   <- rank(df$x)     / (n + 1)
  df$U_t_1 <- rank(df$x_lag) / (n + 1)
  df
}

# --------------------------------------------------------------------------
# Individual plot functions
# Each seasonal function receives a list `entry` with $values and $time
# --------------------------------------------------------------------------

plot_time_series <- function(entry, season, ...) {
  plot(entry$time, entry$values,
       type = "l",
       main = paste("Time series -", config$data_key, "/", season),
       xlab = "Time", ylab = config$xlab)
}

plot_histogram <- function(entry, season, ...) {
  hist(entry$values,
       breaks = config$hist_breaks,
       main   = paste("Histogram -", config$data_key, "/", season),
       xlab   = config$xlab)
}

plot_acf_panel <- function(entry, season, ...) {
  acf(entry$values, main = paste("ACF -", config$data_key, "/", season))
}

plot_lag_panel <- function(entry, season, ...) {
  lag_plot(entry$values, main = paste("Lag -", config$data_key, "/", season))
}

plot_copula_panel <- function(entry, season, ...) {
  df <- make_copula_df(entry$values)
  plot(df$U_t_1, df$U_t,
       pch  = 16, cex = 0.5,
       main = paste("Copula -", config$data_key, "/", season),
       xlab = expression(U[t-1]),
       ylab = expression(U[t]))
}

plot_copula_tail_panel <- function(entry, season, ...) {
  thr <- config$copula_tail_threshold
  df  <- make_copula_df(entry$values)
  df  <- df[df$U_t > thr & df$U_t_1 > thr, ]
  plot(df$U_t_1, df$U_t,
       pch   = 16, cex = 0.5,
       xlim  = c(thr, 1), ylim = c(thr, 1),
       main  = paste0("Upper tail (>", thr, ") - ", season),
       xlab  = expression(U[t-1]),
       ylab  = expression(U[t]))
}

# --------------------------------------------------------------------------
# Dispatch — each plot type wrapped in with_pdf()
# --------------------------------------------------------------------------

for (p in config$plots) {
  switch(p,

    time_series = with_pdf("time_series", {
      with_season_panels(seasons, plot_time_series)
    }),

    histogram = with_pdf("histogram", {
      with_season_panels(seasons, plot_histogram)
    }),

    acf = with_pdf("acf", {
      with_season_panels(seasons, plot_acf_panel)
    }),

    lag = with_pdf("lag", {
      with_season_panels(seasons, plot_lag_panel)
    }),

    copula = with_pdf("copula", {
      with_season_panels(seasons, plot_copula_panel)
    }),

    copula_upper_tail = with_pdf("copula_upper_tail", {
      with_season_panels(seasons, plot_copula_tail_panel)
    }),

    warning("Unknown plot type '", p, "' — skipping.")
  )
}