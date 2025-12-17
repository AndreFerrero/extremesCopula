# =============================================================================
# Compare computational cost: BSL vs SemiBSL inner likelihood
# =============================================================================

source("libs/packages.R")
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")
source("libs/models/builders/simulator.r")
source("libs/models/builders/synthetic_gaussian_bsl.r")
source("libs/models/builders/synthetic_semibsl.r")

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
set.seed(123)

param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

true_param <- c(mu = 0, sigma = 1, theta = 2)
n_obs <- 1000
sum_stats <- function(x) {
  m <- median(x)
  md <- mad(x)
  if (md == 0) md <- 1e-6
  smax <- (max(x) - m) / md
  c(m, md, smax)
}

# Observed summaries
y_obs <- sum_stats(simulator(true_param, n_obs))

# -----------------------------------------------------------------------------
# Benchmark across n_sim
# -----------------------------------------------------------------------------

n_sim_grid <- c(100, 200, 500, 1000, 2000)
n_rep <- 50

results <- data.frame()

for (n_sim in n_sim_grid) {

  times_bsl <- numeric(n_rep)
  times_semibsl <- numeric(n_rep)

  for (r in 1:n_rep) {

    sim_stats <- matrix(NA, n_sim, length(y_obs))
    for (i in 1:n_sim) {
      sim_stats[i, ] <- sum_stats(simulator(true_param, n_obs))
    }

    times_bsl[r] <- system.time(
      synthetic_gaussian_bsl(sim_stats, y_obs)
    )["elapsed"]

    times_semibsl[r] <- system.time(
      synthetic_semibsl(sim_stats, y_obs)
    )["elapsed"]
  }

  results <- rbind(
    results,
    data.frame(
      n_sim = n_sim,
      method = "BSL",
      mean_time = mean(times_bsl)
    ),
    data.frame(
      n_sim = n_sim,
      method = "SemiBSL",
      mean_time = mean(times_semibsl)
    )
  )
}

print(results)

# =============================================================================
# Compare computational cost: BSL vs SemiBSL inner likelihood
# =============================================================================

source("libs/packages.R")
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")
source("libs/models/builders/simulator.r")
source("libs/models/builders/synthetic_gaussian_bsl.r")
source("libs/models/builders/synthetic_semibsl.r")

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
set.seed(123)

param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

true_param <- c(mu = 0, sigma = 1, theta = 2)
n_obs <- 1000
sum_stats <- function(x) {
  m <- median(x)
  md <- mad(x)
  if (md == 0) md <- 1e-6
  smax <- (max(x) - m) / md
  c(m, md, smax)
}

# Observed summaries
y_obs <- sum_stats(simulator(true_param, n_obs))

# -----------------------------------------------------------------------------
# Benchmark across n_sim
# -----------------------------------------------------------------------------

n_sim_grid <- c(100, 200, 500, 1000, 2000)
n_rep <- 50

results <- data.frame()

for (n_sim in n_sim_grid) {

  times_bsl <- numeric(n_rep)
  times_semibsl <- numeric(n_rep)

  for (r in 1:n_rep) {

    sim_stats <- matrix(NA, n_sim, length(y_obs))
    for (i in 1:n_sim) {
      sim_stats[i, ] <- sum_stats(simulator(true_param, n_obs))
    }

    times_bsl[r] <- system.time(
      synthetic_gaussian_bsl(sim_stats, y_obs)
    )["elapsed"]

    times_semibsl[r] <- system.time(
      synthetic_semibsl(sim_stats, y_obs)
    )["elapsed"]
  }

  results <- rbind(
    results,
    data.frame(
      n_sim = n_sim,
      method = "BSL",
      mean_time = mean(times_bsl)
    ),
    data.frame(
      n_sim = n_sim,
      method = "SemiBSL",
      mean_time = mean(times_semibsl)
    )
  )
}

print(results)
