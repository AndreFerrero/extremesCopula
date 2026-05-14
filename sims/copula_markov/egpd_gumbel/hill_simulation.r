# ============================================================
# Simulation study for Hill and Bias-Corrected Hill estimators
# under the EGPD-Gumbel copula Markov model
# ============================================================

# --- Load model components ---
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/hill_estimator.r")

# --- Build model ---
egpd_gumbel_model <- make_copula_markov_model(
  margin_egpd,
  copula_gumbel,
  stan_mod = NULL
)

# ============================================================
# Simulation settings
# ============================================================

set.seed(123)

n <- 2000
mc_it <- 300

theta_grid <- c(1, 2, 4, 6)

# True tail index xi
true_xi <- 0.1

margin_param <- c(
  mu = 0,
  kappa = 2,
  sigma = 1,
  xi = true_xi
)

# Range of k values
k_grid <- seq(5, 300, by = 5)

# ============================================================
# Containers
# ============================================================

results <- list()

# ============================================================
# Simulation loop
# ============================================================

for (theta in theta_grid) {
  cat("\n====================================\n")
  cat("Running theta =", theta, "\n")
  cat("====================================\n")

  hill_mat <- matrix(
    NA,
    nrow = mc_it,
    ncol = length(k_grid)
  )

  hill_bc_mat <- matrix(
    NA,
    nrow = mc_it,
    ncol = length(k_grid)
  )

  # ----------------------------------------------------------
  # Monte Carlo replications
  # ----------------------------------------------------------

  for (mc in 1:mc_it) {
    cat("Replication:", mc, "\n")

    # Simulate sample
    sim_data <- egpd_gumbel_model$simulate(
      n = n,
      margin_param = margin_param,
      copula_param = theta
    )

    x <- sim_data$x

    # Ensure positivity for Hill estimator
    x <- x[x > 0]

    # --------------------------------------------------------
    # Compute estimators over k
    # --------------------------------------------------------

    for (j in seq_along(k_grid)) {
      k <- k_grid[j]

      est <- tryCatch(
        hill_bc_hat(x, k),
        error = function(e) c(Hill = NA, Hill_BC = NA)
      )

      hill_mat[mc, j] <- est["Hill"]
      hill_bc_mat[mc, j] <- est["Hill_BC"]
    }
  }

  # ==========================================================
  # Compute Bias and RMSE
  # ==========================================================

  hill_bias <- colMeans(hill_mat - true_xi, na.rm = TRUE)

  hill_rmse <- sqrt(
    colMeans((hill_mat - true_xi)^2, na.rm = TRUE)
  )

  hill_rel_bias <- abs(colMeans(hill_mat / true_xi, na.rm = TRUE) - 1)

  hill_bc_bias <- colMeans(hill_bc_mat - true_xi, na.rm = TRUE)

  hill_bc_rmse <- sqrt(
    colMeans((hill_bc_mat - true_xi)^2, na.rm = TRUE)
  )

  hill_bc_rel_bias <- abs(colMeans(hill_bc_mat / true_xi, na.rm = TRUE) - 1)

  # Store results
  results[[paste0("theta_", theta)]] <- list(
    k = k_grid,
    hill_bias = hill_bias,
    hill_rmse = hill_rmse,
    hill_rel_bias = hill_rel_bias,
    hill_bc_bias = hill_bc_bias,
    hill_bc_rmse = hill_bc_rmse,
    hill_bc_rel_bias = hill_bc_rel_bias,
    hill_samples = hill_mat,
    hill_bc_samples = hill_bc_mat
  )
}

save(results, file = "sims/copula_markov/egpd_gumbel/res/hill_sim_n2000_mc300_kappa2_xi01.Rdata")
