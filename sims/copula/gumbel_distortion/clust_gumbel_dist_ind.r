library(here)
library(tidyverse)
library(copula)
library(evd)
library(future)
library(future.apply)

set.seed(123)

# ============================
# Paths
# ============================
work_dir <- "/home/ferreroa/work_rainstsimu/extrCopula"
dist_dir <- here(work_dir, "sims", "copula", "gumbel_distortion")
plots_dir <- here(dist_dir, "plots")
res_dir   <- here(dist_dir, "res")
common_dir <- here(work_dir, "sims", "common")

# ============================
# Helper functions
# ============================
source(here(common_dir, "handy_funs.r"))  # provides rCopFrechet etc.

# ============================
# Simulation settings
# ============================
n <- 100000   # size of each sample
B <- 10000     # number of maxima per repetition
R <- 100     # number of repetitions
alpha <- 2
scenarios <- c("iid", "theta_1.5", "theta_2.5")

# ============================
# Parallel plan (over repetitions)
# ============================
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
plan(multisession, workers = n_cores)  # PSOCK, robust on Linux

# ============================
# Container for results
# ============================
all_results_by_scenario <- list()

# ============================
# Simulation loop
# ============================
for (scen in scenarios) {
  cat("Starting scenario:", scen, "...\n")

  # Parallelise only over repetitions R
  scenario_results <- future_lapply(seq_len(R), function(r) {

    # Generate B maxima sequentially
    max_vec <- numeric(B)

    if (scen == "iid") {
      for (i in seq_len(B)) {
        X <- rfrechet(n, shape = alpha)
        max_vec[i] <- max(X)
      }
    } else {
      tval <- as.numeric(sub("theta_", "", scen))
      for (i in seq_len(B)) {
        Gcop <- gumbelCopula(param = tval, dim = n)
        X <- rCopFrechet(alpha, Gcop)
        max_vec[i] <- max(X)
      }
    }

    # Fit GEV safely
    fit <- tryCatch(fgev(max_vec, std.err = FALSE)$estimate,
                    error = function(e) c(NA, NA, NA))

    data.frame(
      repetition = r,
      location = fit[1],
      scale    = fit[2],
      shape    = fit[3]
    )

  }, future.seed = TRUE)  # parallel-safe RNG

  all_results_by_scenario[[scen]] <- bind_rows(scenario_results)
}

# ============================
# Combine and save
# ============================
all_results_long <- bind_rows(lapply(names(all_results_by_scenario), function(scen) {
  df <- all_results_by_scenario[[scen]]
  df$scenario <- scen
  df
}))

save(all_results_long, file = here(res_dir, "clust_gumbel_wholemax_res.Rdata"))

# ============================
# Shutdown workers
# ============================
plan(sequential)

# ============================================================
# Summaries and quick plots
# ============================================================
# Convert to long format for plotting distributions of parameter estimates across R repetitions
all_long <- all_results_long %>%
  pivot_longer(cols = c(location, scale, shape), names_to = "parameter", values_to = "estimate")

# ============================================================
# Remove extreme outliers before plotting
# ============================================================
clean_long <- all_long

clean_long$keep <- TRUE   # logical flag

params <- unique(clean_long$parameter)

for (p in params) {
  idx <- clean_long$parameter == p
  est <- clean_long$estimate[idx]

  # compute quantile limits for *this* parameter
  q_low  <- quantile(est, 0.005, na.rm = TRUE)
  q_high <- quantile(est, 0.995, na.rm = TRUE)

  # mark as FALSE anything outside the limits
  clean_long$keep[idx] <- est >= q_low & est <= q_high
}

# keep filtered rows only
clean_long <- clean_long[clean_long$keep, ]

# remove the temporary flag
clean_long$keep <- NULL


# Example boxplot across scenarios (one point per repetition)
my_cols <- c("iid" = "#66c2a5", "theta_1.5" = "#fc8d62", "theta_2.5" = "#8da0cb")

# Create the plot object
p <- ggplot(clean_long, aes(x = scenario, y = estimate, fill = scenario)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~parameter, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = my_cols) +
  labs(
    title = "Distribution of GEV parameter estimates over R repetitions (whole-sequence maxima)",
    subtitle = paste("Each repetition: generate", B, "maxima (each max from n =", n, "observations). R =", R),
    x = NULL, y = "Estimated parameter"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

save(p, file = here(plots_dir, "sim_plot.Rdata"))