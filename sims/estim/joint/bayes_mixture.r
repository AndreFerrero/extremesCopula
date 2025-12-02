library(copula)
library(coda)
library(mvtnorm)
library(future) # For parallel processing
library(future.apply)

set.seed(123)

res_dir <- here("sims", "estim", "joint", "res")
# -----------------------------------------------------------
# 1. Simulate Data
# -----------------------------------------------------------
n <- 150
theta_true <- 3
mu_true <- 0
sigma_true <- 1

cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))
X <- qlnorm(U, meanlog = mu_true, sdlog = sigma_true)
Y <- log(X)

# -----------------------------------------------------------
# 2. Model Definitions
# -----------------------------------------------------------
d_mix <- function(y, w, m1, s1, m2, s2) {
  w * dnorm(y, m1, s1) + (1 - w) * dnorm(y, m2, s2)
}

p_mix <- function(y, w, m1, s1, m2, s2) {
  w * pnorm(y, m1, s1) + (1 - w) * pnorm(y, m2, s2)
}

log_posterior <- function(param_vec, data_y) {
  # param_vec: theta, w, mu1, s1, mu2, s2
  theta <- param_vec[1]
  w <- param_vec[2]
  mu1 <- param_vec[3]
  s1 <- param_vec[4]
  mu2 <- param_vec[5]
  s2 <- param_vec[6]

  # --- PRIORS ---
  if (theta <= 1.01 || theta > 20) {
    return(-Inf)
  }
  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

  if (w <= 0.01 || w >= 0.99) {
    return(-Inf)
  }
  if (mu1 > mu2) {
    return(-Inf)
  } # Ordering constraint

  if (s1 <= 0.05 || s2 <= 0.05) {
    return(-Inf)
  }
  lp_s <- dgamma(s1, 3, 2, log = TRUE) + dgamma(s2, 3, 2, log = TRUE)

  lp_mu <- dnorm(mu1, 0, 5, log = TRUE) + dnorm(mu2, 0, 5, log = TRUE)

  # --- LIKELIHOOD ---
  dens_vals <- d_mix(data_y, w, mu1, s1, mu2, s2)
  if (any(dens_vals <= 1e-300)) {
    return(-Inf)
  }
  ll_margin <- sum(log(dens_vals))

  u_hat <- p_mix(data_y, w, mu1, s1, mu2, s2)
  u_hat <- pmin(pmax(u_hat, 1e-6), 1 - 1e-6)

  d_cop <- dCopula(u_hat,
    copula = gumbelCopula(theta, dim = length(data_y)),
    log = TRUE
  )
  if (is.na(d_cop) || d_cop <= 0) {
    return(-Inf)
  }

  return(lp_theta + lp_s + lp_mu + ll_margin + log(d_cop))
}

# -----------------------------------------------------------
# 3. Parallel Worker Function
# -----------------------------------------------------------
run_adaptive_chain_parallel <- function(data, n_iter, init_vals, chain_id, progress_dir) {
  p_names <- c("theta", "w", "mu1", "s1", "mu2", "s2")
  n_par <- length(p_names)
  chain <- matrix(NA, nrow = n_iter, ncol = n_par)
  colnames(chain) <- p_names

  curr_par <- init_vals
  curr_lp <- log_posterior(curr_par, data)

  # Tuning vars
  prop_sd <- c(0.2, 0.05, 0.1, 0.1, 0.1, 0.1)
  batch_size <- 50
  accept_batch <- numeric(n_par)
  total_accepts <- 0

  # Define status file path for this chain
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    for (j in 1:n_par) {
      prop_val <- curr_par
      prop_val[j] <- rnorm(1, curr_par[j], prop_sd[j])
      prop_lp <- log_posterior(prop_val, data)
      ratio <- prop_lp - curr_lp

      if (is.finite(ratio) && log(runif(1)) < ratio) {
        curr_par <- prop_val
        curr_lp <- prop_lp
        accept_batch[j] <- accept_batch[j] + 1
        total_accepts <- total_accepts + 1
      }
    }
    chain[i, ] <- curr_par

    # Adaptation
    if (i <= n_iter / 2 && i %% batch_size == 0) {
      acc_rates <- accept_batch / batch_size
      for (k in 1:n_par) {
        if (acc_rates[k] < 0.2) {
          prop_sd[k] <- prop_sd[k] * 0.8
        } else if (acc_rates[k] > 0.3) prop_sd[k] <- prop_sd[k] * 1.2
      }
      accept_batch <- numeric(n_par)
    }

    # --- WRITE STATUS TO FILE (For the Main Process to read) ---
    # We update file every 50 iterations to avoid I/O bottlenecks
    if (i %% 50 == 0 || i == n_iter) {
      rate <- (total_accepts / (i * n_par)) * 100
      pct <- as.integer((i / n_iter) * 100)
      # Atomic write is hard, but simple writeLines is usually fine here
      try(
        {
          writeLines(sprintf("%d|%.1f", pct, rate), status_file)
        },
        silent = TRUE
      )
    }
  }
  return(mcmc(chain))
}

# -----------------------------------------------------------
# 4. Setup Parallel Execution & Monitoring
# -----------------------------------------------------------

# A. Configure Parallel Backend
# 'multisession' works on Windows, Mac, and Linux
plan(multisession, workers = 3)

cat("Setting up parallel workers...\n")

# B. Setup IPC (Inter-Process Communication) via Temp Files
progress_dir <- tempdir()
# Clean up old status files if any
file.remove(list.files(progress_dir, pattern = "status_", full.names = TRUE))

n_chains <- 3
n_iter <- 5000
burn_in <- 1000

# Prepare Initial Values
inits_list <- list()
for (c in 1:n_chains) {
  inits_list[[c]] <- c(
    theta = runif(1, 1.5, 4), w = 0.5,
    mu1 = mean(Y) - 0.2, s1 = sd(Y),
    mu2 = mean(Y) + 0.2, s2 = sd(Y)
  )
}

# C. Launch Workers Asynchronously
# We use future() explicitly so we can return control to main loop immediately
futures_list <- list()
for (c in 1:n_chains) {
  futures_list[[c]] <- future(
    {
      run_adaptive_chain_parallel(Y, n_iter, inits_list[[c]], c, progress_dir)
    },
    seed = TRUE
  ) # seed=TRUE ensures robust parallel RNG
}

cat("MCMC Started. Monitoring progress...\n")

# D. The "Dashboard" Loop
# This runs in the main process while workers churn in background
while (!all(resolved(futures_list))) {
  status_texts <- character(n_chains)

  for (c in 1:n_chains) {
    fpath <- file.path(progress_dir, paste0("status_", c, ".txt"))
    if (file.exists(fpath)) {
      # Read the file (Percent|Rate)
      # suppressWarnings in case file is being written to at that exact moment
      info <- suppressWarnings(try(readLines(fpath, n = 1), silent = TRUE))
      if (inherits(info, "try-error") || length(info) == 0) {
        status_texts[c] <- "Init..."
      } else {
        parts <- strsplit(info, "\\|")[[1]]
        status_texts[c] <- sprintf("C%d: %3s%% (Acc: %4s%%)", c, parts[1], parts[2])
      }
    } else {
      status_texts[c] <- sprintf("C%d: Pending...", c)
    }
  }

  # Print consolidated line with \r to overwrite
  # Example: "C1: 10% (Acc: 23.1%) | C2: 10% (Acc: 24.0%) | C3: 9% (Acc: 21.5%)"
  cat("\r", paste(status_texts, collapse = " | "))
  flush.console()

  Sys.sleep(0.5) # Refresh every half second
}

# Final aesthetic newline
cat("\nAll chains finished.\n")

# E. Collect Results
chains_list <- value(futures_list) # Retrieve data from futures
save(chains_list, file = here(res_dir, "bayes_chains_list.Rdata"))
load(here(res_dir, "bayes_chains_list.Rdata"))

mcmc_obj <- mcmc.list(chains_list)
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 10)


# -----------------------------------------------------------
# 5. Diagnostics
# -----------------------------------------------------------
cat("\n--- DIAGNOSTICS ---\n")

gelman_res <- gelman.diag(mcmc_clean)
print(gelman_res)

eff_size <- effectiveSize(mcmc_clean)
cat("\nEffective Sample Sizes (Theta):", round(eff_size["theta"]), "\n")

# -----------------------------------------------------------
# 6. Plotting
# -----------------------------------------------------------
par(mfrow = c(2, 1))
plot(mcmc_clean[, "theta"], main = "Traceplot: Theta", auto.layout = FALSE)
par(mfrow = c(1, 1))

dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(dens_est,
  main = "Posterior Density of Theta", lwd = 2, col = "blue",
  xlab = expression(theta)
)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

mean(as.matrix(mcmc_clean[, "theta"]))
median(as.matrix(mcmc_clean[, "theta"]))

##########################################################################################
# -----------------------------------------------------------
# 1. Define the JOINT Negative Log-Likelihood
# -----------------------------------------------------------
# params: c(w, mu1, s1, mu2, s2)  <-- Theta is fixed externally for profiling
negLL_joint_profile <- function(params, x_data, fixed_theta) {
  w <- params[1]
  mu1 <- params[2]
  s1 <- params[3]
  mu2 <- params[4]
  s2 <- params[5]

  # 1. Marginal Likelihood
  d_vals <- w * dnorm(x_data, mu1, s1) + (1 - w) * dnorm(x_data, mu2, s2)
  if (any(d_vals <= 0)) {
    return(1e10)
  }
  ll_margin <- sum(log(d_vals))

  # 2. Transform to U
  u_hat <- w * pnorm(x_data, mu1, s1) + (1 - w) * pnorm(x_data, mu2, s2)
  u_hat <- pmin(pmax(u_hat, 1e-6), 1 - 1e-6)

  # 3. Copula Likelihood
  # Note: Theta is fixed here!
  d_cop <- dCopula(u_hat, copula = gumbelCopula(fixed_theta, dim = length(x_data)))

  if (any(is.na(d_cop)) || any(d_cop <= 0)) {
    return(1e10)
  }
  ll_cop <- sum(log(d_cop))

  return(-(ll_margin + ll_cop))
}

# -----------------------------------------------------------
# 2. Compute Frequentist Profile Likelihood Curve
# -----------------------------------------------------------
cat("Computing Frequentist Joint Profile Likelihood...\n")

# Grid of Thetas to profile
# (Matches the range of your Bayesian posterior)
theta_seq <- seq(1.2, 5, length.out = 30)
prof_lik <- numeric(length(theta_seq))

# Initial Guesses for Nuisance Params (w, mu1, s1, mu2, s2)
# We use sensible stats from the data to give MLE a fighting chance
start_nuisance <- c(0.5, mean(Y) - 0.5, sd(Y), mean(Y) + 0.5, sd(Y))

# HARD BOUNDS (The Frequentist "Prior")
# We forbid sigma < 0.2 to prevent singularity crashes.
lower_b <- c(0.01, -10, 0.2, -10, 0.2)
upper_b <- c(0.99, 10, 5.0, 10, 5.0)

cat("Starting Profile Likelihood Calculation (Total steps: ", length(theta_seq), ")\n")

for (i in 1:length(theta_seq)) {
  t_val <- theta_seq[i]

  # Optimize nuisance params for THIS theta
  fit <- try(optim(
    par = start_nuisance,
    fn = negLL_joint_profile,
    x_data = Y,
    fixed_theta = t_val,
    method = "L-BFGS-B",
    lower = lower_b,
    upper = upper_b,
    control = list(maxit = 1000) # Safety limit on iterations
  ), silent = TRUE)

  if (inherits(fit, "try-error")) {
    prof_lik[i] <- -1e10
    # Don't update start_nuisance if it failed
  } else {
    prof_lik[i] <- -fit$value

    # Warm start: use this solution as start for the next theta
    # This speeds up convergence significantly
    start_nuisance <- fit$par
  }

  # --- PROGRESS INDICATOR ---
  # \r overwrites the line, keeping the console clean
  cat(sprintf(
    "\rStep %d/%d | Current Theta: %.2f | LogLik: %.2f",
    i, length(theta_seq), t_val, prof_lik[i]
  ))
  flush.console()
}
cat("\nDone.\n")

save(prof_lik, file = here(res_dir, "prof_lik_optim_mixture.Rdata"))
# -----------------------------------------------------------
# 3. Normalize for Comparison
# -----------------------------------------------------------
# Convert Log-Likelihood to "Likelihood"
lik_raw <- exp(prof_lik - max(prof_lik))

# Normalize area to 1 so it looks like a PDF
area <- sum(lik_raw) * (theta_seq[2] - theta_seq[1])
freq_profile_dens <- lik_raw / area

# -----------------------------------------------------------
# 4. THE VISUAL COMPARISON
# -----------------------------------------------------------
par(mfrow = c(1, 1))

# Bayesian Posterior
dens_bayes <- density(as.matrix(mcmc_clean)[, "theta"])

# Plot Setup
plot(dens_bayes,
  lwd = 3, col = "blue",
  main = "Joint Estimation: Bayesian vs. Frequentist",
  xlab = expression(theta),
  ylim = c(0, max(c(dens_bayes$y, freq_profile_dens)) * 1.2),
  xlim = c(1.2, 5)
)

# Frequentist Profile Likelihood
lines(theta_seq, freq_profile_dens, col = "red", lwd = 3, lty = 2)

# True Value
abline(v = theta_true, col = "darkgreen", lwd = 2)

legend("topright",
  legend = c("Bayesian Joint Posterior", "Freq. Joint Profile Likelihood", "True Theta"),
  col = c("blue", "red", "darkgreen"),
  lwd = c(3, 3, 2),
  lty = c(1, 2, 1)
)

# -----------------------------------------------------------
# 5. Peak Comparison
# -----------------------------------------------------------
bayes_peak <- dens_bayes$x[which.max(dens_bayes$y)]
freq_peak <- theta_seq[which.max(freq_profile_dens)]

cat("\n--- COMPARISON ---\n")
cat("True Theta:       ", theta_true, "\n")
cat("Bayesian Mode:    ", round(bayes_peak, 3), "\n")
cat("Frequentist MLE:  ", round(freq_peak, 3), "\n")
