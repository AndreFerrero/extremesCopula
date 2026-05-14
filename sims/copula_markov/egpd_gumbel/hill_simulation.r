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
mc_it <- 1

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
k_grid <- seq(20, 300, by = 10)

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

  hill_bc_bias <- colMeans(hill_bc_mat - true_xi, na.rm = TRUE)

  hill_bc_rmse <- sqrt(
    colMeans((hill_bc_mat - true_xi)^2, na.rm = TRUE)
  )

  # Store results
  results[[paste0("theta_", theta)]] <- list(
    k = k_grid,
    hill_bias = hill_bias,
    hill_rmse = hill_rmse,
    hill_bc_bias = hill_bc_bias,
    hill_bc_rmse = hill_bc_rmse,
    hill_samples = hill_mat,
    hill_bc_samples = hill_bc_mat
  )
}

# ============================================================
# Plotting
# ============================================================

par(mfrow = c(2, 2))

# ------------------------------------------------------------
# Bias plots
# ------------------------------------------------------------

for (theta in theta_grid) {

  

  plot(
    res$k,
    res$hill_bias,
    type = "l",
    lwd = 2,
    ylim = ylim_range,
    xlab = "k",
    ylab = "Bias",
    main = paste("Bias - theta =", theta)
  )

  lines(
    res$k,
    res$hill_bc_bias,
    lwd = 2,
    lty = 2
  )

  abline(h = 0, col = "gray")

  legend(
    "topright",
    legend = c("Hill", "Hill BC"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}

# ------------------------------------------------------------
# RMSE plots
# ------------------------------------------------------------

dev.new()

par(mfrow = c(2, 2))

for (theta in theta_grid) {

  res <- results[[paste0("theta_", theta)]]

  ylim_range <- range(
    c(res$hill_rmse, res$hill_bc_rmse),
    na.rm = TRUE
  )

  plot(
    res$k,
    res$hill_rmse,
    type = "l",
    lwd = 2,
    ylim = ylim_range,
    xlab = "k",
    ylab = "RMSE",
    main = paste("RMSE - theta =", theta)
  )

  lines(
    res$k,
    res$hill_bc_rmse,
    lwd = 2,
    lty = 2
  )

  legend(
    "topright",
    legend = c("Hill", "Hill BC"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}

# ============================================================
# Optional: combine all theta values on one plot
# ============================================================

dev.new()

plot(
  results$theta_1$k,
  results$theta_1$hill_rmse,
  type = "l",
  lwd = 2,
  ylim = range(
    sapply(results, function(x) x$hill_rmse),
    na.rm = TRUE
  ),
  xlab = "k",
  ylab = "RMSE",
  main = "Hill RMSE across dependence levels"
)

for (theta in theta_grid[-1]) {

  lines(
    results[[paste0("theta_", theta)]]$k,
    results[[paste0("theta_", theta)]]$hill_rmse,
    lwd = 2
  )
}

legend(
  "topright",
  legend = paste("theta =", theta_grid),
  lwd = 2,
  bty = "n"
)

# ============================================================
# Example summary table
# ============================================================

summary_table <- do.call(
  rbind,
  lapply(theta_grid, function(theta) {

    res <- results[[paste0("theta_", theta)]]

    best_k_hill <- res$k[which.min(res$hill_rmse)]
    best_k_bc <- res$k[which.min(res$hill_bc_rmse)]

    data.frame(
      theta = theta,
      best_k_hill = best_k_hill,
      min_rmse_hill = min(res$hill_rmse, na.rm = TRUE),
      best_k_hill_bc = best_k_bc,
      min_rmse_hill_bc = min(res$hill_bc_rmse, na.rm = TRUE)
    )
  })
)

print(summary_table)