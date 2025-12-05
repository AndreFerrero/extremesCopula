library(copula)
library(evd)
library(coda)
library(mvtnorm)
library(future)
library(future.apply)

set.seed(123)

# -----------------------------------------------------------
# 1. Simulate Data (Fréchet Margins)
# -----------------------------------------------------------
n <- 150
theta_true <- 3

# Fréchet Parameters
# shape (alpha), location (mu), scale (sigma)
alpha_true <- 3 # Shape (Tail index)
mu_true <- 10 # Location (Lower bound)
sigma_true <- 5 # Scale

cat("1. Simulating Data (Fréchet Margins, n =", n, ")...\n")

# A. Dependence
cop <- gumbelCopula(param = theta_true, dim = n)
U <- as.numeric(rCopula(1, cop))

# B. Margins (Fréchet Quantile Function)
# X = mu + sigma * (-log(U))^(-1/alpha)
X <- mu_true + sigma_true * (-log(U))^(-1 / alpha_true)

# -----------------------------------------------------------
# 2. Model Definitions (Fréchet + Gumbel)
# -----------------------------------------------------------

log_posterior <- function(param_vec, data) {
  # Unpack parameters
  theta <- param_vec[1]
  mu <- param_vec[2]
  sigma <- param_vec[3]
  alpha <- param_vec[4]

  # --- A. CONSTRAINTS (CRITICAL) ---

  # 1. Copula constraints
  if (theta <= 1.01 || theta > 20) {
    return(-Inf)
  }

  # 2. Fréchet Constraints
  if (sigma <= 0.01) {
    return(-Inf)
  }
  if (alpha <= 0.1 || alpha > 20) {
    return(-Inf)
  }

  # 3. SUPPORT CONSTRAINT:
  # Fréchet is only defined for x > mu.
  # If any data point is <= mu, likelihood is 0 (log is -Inf).
  if (mu >= min(data)) {
    return(-Inf)
  }

  # --- B. PRIORS ---

  lp_theta <- dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)

  # Location (mu): Normal, but effectively truncated by the data min
  lp_mu <- dnorm(mu, mean = mean(data) - sd(data), sd = 10, log = TRUE)

  # Scale (sigma): Boundary Avoiding LogNormal
  lp_sigma <- dlnorm(sigma, meanlog = log(sd(data) / 2), sdlog = 0.5, log = TRUE)

  # Shape (alpha): Gamma (favoring values around 2-5 typical for extremes)
  lp_alpha <- dgamma(alpha, shape = 3, rate = 1, log = TRUE)

  # --- C. MARGINAL LIKELIHOOD ---
  dens_vals <- dfrechet(data, mu, sigma, alpha)

  # Safety for log(0)
  if (any(dens_vals <= 0)) {
    return(-Inf)
  }
  ll_margin <- sum(log(dens_vals))

  # --- D. COPULA LIKELIHOOD ---
  u_hat <- pfrechet(data, mu, sigma, alpha)
  u_hat <- pmin(pmax(u_hat, 1e-8), 1 - 1e-8)

  ld_cop <- dCopula(u_hat,
    copula = gumbelCopula(theta, dim = length(data)),
    log = TRUE
  )

  if (is.na(ld_cop) || !is.finite(ld_cop)) {
    return(-Inf)
  }

  return(lp_theta + lp_mu + lp_sigma + lp_alpha + ll_margin + ld_cop)
}

# -----------------------------------------------------------
# 3. Worker: Block Adaptive Metropolis
# -----------------------------------------------------------
run_worker <- function(data, n_iter, init_vals, chain_id, progress_dir) {
  library(copula)
  library(mvtnorm)
  library(coda)
  library(evd)

  p_names <- c("theta", "mu", "sigma", "alpha")
  n_par <- length(p_names)
  chain <- matrix(NA, nrow = n_iter, ncol = n_par)
  colnames(chain) <- p_names

  curr_par <- init_vals

  # Robust Start
  curr_lp <- log_posterior(curr_par, data)
  if (!is.finite(curr_lp)) {
    stop("Worker failed: Initial values yield -Inf posterior (Likely mu >= min(x)).")
  }

  # Adaptive Settings
  cov_mat <- diag(rep(0.01, n_par))
  # Optimal Scaling for d=4
  sd_scale <- (2.38^2) / n_par
  eps <- 1e-6 * diag(n_par)

  adapt_start <- 500
  total_accepts <- 0
  status_file <- file.path(progress_dir, paste0("status_", chain_id, ".txt"))

  for (i in 1:n_iter) {
    # Block Proposal
    prop_val <- rmvnorm(1, mean = curr_par, sigma = cov_mat)[1, ]

    prop_lp <- log_posterior(prop_val, data)
    ratio <- prop_lp - curr_lp

    if (is.finite(ratio) && log(runif(1)) < ratio) {
      curr_par <- prop_val
      curr_lp <- prop_lp
      total_accepts <- total_accepts + 1
    }
    chain[i, ] <- curr_par

    # Adaptation
    if (i > adapt_start && i %% 50 == 0) {
      emp_cov <- cov(chain[1:i, ])
      cov_mat <- (sd_scale * emp_cov) + eps
    }

    # Dashboard
    if (i %% 50 == 0 || i == n_iter) {
      rate <- (total_accepts / i) * 100
      pct <- as.integer((i / n_iter) * 100)
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
# 4. Execution
# -----------------------------------------------------------
n_chains <- 3
plan(multisession, workers = n_chains)
progress_dir <- tempdir()
file.remove(list.files(progress_dir, pattern = "status_", full.names = TRUE))

n_iter <- 15000

inits_list <- list()
for (c in 1:n_chains) {
  # Careful Initialization for mu: Must be < min(X)
  data_min <- min(X)

  inits_list[[c]] <- c(
    theta = runif(1, 2, 4),
    mu    = data_min - runif(1, 0.5, 2), # Ensure below min
    sigma = runif(1, 2, 6),
    alpha = runif(1, 2, 4)
  )
}

futures_list <- list()
for (c in 1:n_chains) {
  futures_list[[c]] <- future(
    {
      run_worker(X, n_iter, inits_list[[c]], c, progress_dir)
    },
    seed = TRUE
  )
}

cat("MCMC Started (Fréchet Margins)...\n")

# Robust Monitor
while (!all(resolved(futures_list))) {
  status_texts <- character(n_chains)
  for (c in 1:n_chains) {
    fpath <- file.path(progress_dir, paste0("status_", c, ".txt"))
    if (file.exists(fpath)) {
      info <- suppressWarnings(try(readLines(fpath, n = 1), silent = TRUE))
      if (inherits(info, "try-error") || length(info) == 0) {
        status_texts[c] <- "Init..."
      } else {
        parts <- strsplit(info, "\\|")[[1]]
        status_texts[c] <- sprintf("C%d: %3s%% (Acc: %4s%%)", c, parts[1], parts[2])
      }
    } else {
      status_texts[c] <- "Init..."
    }
  }
  cat("\r", paste(status_texts, collapse = " | "))
  flush.console()
  Sys.sleep(0.5)
}
cat("\nDone.\n")

# -----------------------------------------------------------
# 5. Analysis
# -----------------------------------------------------------
res_dir <- here("sims", "estim", "joint", "res")

save(futures_list, file = here(res_dir, "frechet_bayes_chains.Rdata"))
# load(here(res_dir, "frechet_bayes_chains.Rdata"))

chains_list <- value(futures_list)
mcmc_obj <- mcmc.list(chains_list)
burn_in <- 3000
mcmc_clean <- window(mcmc_obj, start = burn_in + 1, thin = 5)

cat("\n--- DIAGNOSTICS ---\n")
print(gelman.diag(mcmc_clean))
print(effectiveSize(mcmc_clean))

# Visuals
par(mfrow = c(2, 2))
plot(mcmc_clean[, "theta"], main = "Theta", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "mu"], main = "Mu (Loc)", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "sigma"], main = "Sigma (Scale)", density = FALSE, auto.layout = FALSE)
plot(mcmc_clean[, "alpha"], main = "Alpha (Shape)", density = FALSE, auto.layout = FALSE)
par(mfrow = c(1, 1))

# DENSITIES
par(mfrow = c(2, 2))
theta_dens_est <- density(as.matrix(mcmc_clean)[, "theta"])
plot(theta_dens_est,
    main = "Posterior Theta",
    lwd = 2, col = "blue", xlab = expression(theta)
)
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

mu_dens_est <- density(as.matrix(mcmc_clean)[, "mu"])
plot(mu_dens_est,
    main = "Posterior Mu",
    lwd = 2, col = "blue", xlab = expression(mu)
)
abline(v = mu_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

sigma_dens_est <- density(as.matrix(mcmc_clean)[, "sigma"])
plot(sigma_dens_est,
    main = "Posterior Sigma",
    lwd = 2, col = "blue", xlab = expression(sigma)
)
abline(v = sigma_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)

alpha_dens_est <- density(as.matrix(mcmc_clean)[, "alpha"])
plot(alpha_dens_est,
    main = "Posterior Alpha",
    lwd = 2, col = "blue", xlab = expression(alpha)
)
abline(v = alpha_true, col = "red", lwd = 2, lty = 2)
legend("topright", c("Posterior", "True"), col = c("blue", "red"), lty = 1:2)
par(mfrow = c(1, 1))



# Blocking analysis: should I block some parameters
library(coda)
library(MASS) # kde2d
library(psych) # pairs.panels (optional, nice viz)

# convert to matrix of samples (chains combined)
mat_all <- as.matrix(mcmc_clean) # or as.matrix(mcmc_clean) if you want post-burnin
pnames <- colnames(mat_all)

# 1.A: Posterior correlation matrix
cor_mat <- cor(mat_all)
print(round(cor_mat, 3))

# 1.B: Spearman (rank) correlations (robust to nonlinearity)
spearman_mat <- cor(mat_all, method = "spearman")
print(round(spearman_mat, 3))

# 1.C: Pairs plot with densities and smoothing (nice quick view)
pairs.panels(mat_all[, c("sigma", "alpha", "mu", "theta")], ellipses = FALSE, lm = TRUE)

# 1.D: Bivariate density contour for sigma vs alpha
z <- MASS::kde2d(mat_all[, "sigma"], mat_all[, "alpha"], n = 80)
plot(mat_all[, "sigma"], mat_all[, "alpha"],
  pch = 20, cex = 0.4,
  xlab = "sigma", ylab = "alpha", main = "sigma vs alpha (posterior)"
)
contour(z, add = TRUE)

z <- MASS::kde2d(mat_all[, "theta"], mat_all[, "alpha"], n = 80)
plot(mat_all[, "theta"], mat_all[, "alpha"],
  pch = 20, cex = 0.4,
  xlab = "theta", ylab = "alpha", main = "theta vs alpha (posterior)"
)
contour(z, add = TRUE)

####################
# GETTING G
# 1. Setup Grid
y_grid <- seq(min(X), max(X) * 4, length.out = 200)
post_mat <- as.matrix(mcmc_clean)
n_sims <- nrow(post_mat)
n_data <- length(X)

# Matrices to store results
G2   <- matrix(0, n_sims, length(y_grid))
G1 <- matrix(0, n_sims, length(y_grid))

# 2. Calculation Loop
for(i in 1:n_sims) {
  # Extract params
  th <- post_mat[i, "theta"]
  mu <- post_mat[i, "mu"]
  si <- post_mat[i, "sigma"]
  al <- post_mat[i, "alpha"]
  
  # --- Calculate F(y) (Marginal) ---
  # Handle support carefully: F(y) = 0 if y <= mu
  z <- (y_grid - mu) / si
  F_y <- numeric(length(y_grid))
  valid <- z > 0
  F_y[valid] <- exp(-z[valid]^(-al))
  
  # Clamp for numerical stability
  F_y <- pmin(pmax(F_y, 1e-15), 1 - 1e-15)
  
  # --- METHOD 2: EXACT FINITE SAMPLE (RECOMMENDED) ---
  # G(y) = psi( n * psi^-1( F(y) ) )
  # psi(t) = exp(-t^(1/th)), psi^-1(u) = (-log u)^th
  
  neg_log_F <- -log(F_y)
  t_val <- neg_log_F^th
  t_sum <- n_data * t_val
  G2[i, ] <- exp(-t_sum^(1/th))
  
  # --- METHOD 1 CORRECTED: ASYMPTOTIC POWER DIAGONAL ---
  # For Gumbel, the effective sample size is n^(1/theta)
  # G_approx(y) = F(y)^(n^(1/theta))
  
  n_eff <- n_data^(1/th)
  G1[i, ] <- F_y^n_eff
}

# 3. True Distribution (for comparison)
z_true <- (y_grid - mu_true) / sigma_true
F_true <- numeric(length(y_grid))
F_true[z_true > 0] <- exp(-z_true[z_true > 0]^(-alpha_true))
G_true <- exp(-(n_data * (-log(F_true))^theta_true)^(1/theta_true))

# Summaries
G1_mean <- colMeans(G1)
G2_mean <- colMeans(G2)

G1_lower <- apply(G1, 2, quantile, 0.05)
G1_upper <- apply(G1, 2, quantile, 0.95)

G2_lower <- apply(G2, 2, quantile, 0.05)
G2_upper <- apply(G2, 2, quantile, 0.95)

# Archi diagonal
par(mfrow = c(1, 1))
plot(y_grid, G2_mean, type = "l", lwd = 1, col = "blue",
     ylim = c(0, 1), ylab = "G", xlab = "y",
     main = "Archimedean copula diagonal - Estimation of G")

polygon(c(y_grid, rev(y_grid)), c(G2_lower, rev(G2_upper)), 
        col = rgb(0, 0, 1, 0.1), border = NA)
# True G
lines(y_grid, G_true, col = "red", lwd = 2, lty = 3)


# Power diagonal

plot(y_grid, G1_mean, type = "l", lwd = 1, col = "green",
     ylim = c(0, 1), ylab = "P(Max <= y)", xlab = "y",
     main = "Comparison of Methods to Recover G")
polygon(c(y_grid, rev(y_grid)), c(G1_lower, rev(G1_upper)), 
        col = rgb(0, 0, 1, 0.1), border = NA)

lines(y_grid, G_true, col = "red", lwd = 2, lty = 3)



# ==============================================================================
# ASYMPTOTIC DIAGONAL APPROACH
# ==============================================================================

# 1. Gumbel Generator
psi_gumbel <- function(t, theta) {
  exp(-t^(1 / theta))
}

# 2. Extract Posterior Matrix
post_mat <- as.matrix(mcmc_clean)
n_sims <- nrow(post_mat)

# 3. Define Evaluation Grid (Raw Scale)
y_grid <- seq(min(X), max(X) * 10, length.out = 200)
G_curves <- matrix(0, nrow = n_sims, ncol = length(y_grid))

# 4. Compute G(y) using Asymptotic H
for (i in 1:n_sims) {
  # A. Get Parameters
  theta <- post_mat[i, "theta"]
  mu    <- post_mat[i, "mu"]
  sigma <- post_mat[i, "sigma"]
  alpha <- post_mat[i, "alpha"]

  # B. Calculate Normalizing Constants (Fréchet Domain)
  # Based on EVT theory for Fréchet(alpha, mu, sigma)
  bn <- mu
  an <- sigma * n^(1 / alpha)

  # C. Normalize the Grid
  # Transform Raw y -> Normalized z
  # y = an * z + bn  ==>  z = (y - bn) / an
  z_norm <- (y_grid - bn) / an
  
  # Handle support: Fréchet limit is 0 for z <= 0
  z_norm[z_norm <= 0] <- NA 

  # D. Calculate Theoretical Limit H(z)
  # The limit of normalized Fréchet maxima is Standard Fréchet Phi_alpha
  # H(z) = exp( -z^(-alpha) )
  # -log(H) = z^(-alpha)
  
  neg_log_H <- - log(pfrechet(z_norm, shape = alpha))
  
  # Handle NAs (out of support)
  neg_log_H[is.na(neg_log_H)] <- Inf # Probability 0 -> -log(0) = Inf

  # E. Apply Distortion: G = psi( -log H )
  G_curves[i, ] <- psi_gumbel(neg_log_H, theta)
}

# -----------------------------------------------------------
# Visualization
# -----------------------------------------------------------

G_mean <- colMeans(G_curves)
G_lower <- apply(G_curves, 2, quantile, 0.025)
G_upper <- apply(G_curves, 2, quantile, 0.975)

# Calculate TRUE Theoretical G
# True constants
bn_true <- mu_true
an_true <- sigma_true * n^(1 / alpha_true)

# True Normalized Grid
z_true <- (y_grid - bn_true) / an_true
z_true[z_true <= 0] <- NA
neg_log_H_true <- z_true^(-alpha_true)
neg_log_H_true[is.na(neg_log_H_true)] <- Inf

# True G
G_true <- psi_gumbel(neg_log_H_true, theta_true)

# PLOT
par(mfrow = c(1, 1))
plot(y_grid, G_mean,
  type = "l", lwd = 2, col = "blue",
  main = "Posterior Maxima G (Asymptotic EVT Method)",
  xlab = "Maximum Value (y)", ylab = "P(Mn <= y)",
  ylim = c(0, 1)
)

polygon(c(y_grid, rev(y_grid)), c(G_lower, rev(G_upper)), 
        col = rgb(0, 0, 1, 0.1), border = NA)

lines(y_grid, G_true, col = "red", lwd = 2, lty = 2)

legend("bottomright",
  legend = c("Estimated G (EVT)", "True G"),
  col = c("blue", "red"), lty = c(1, 2), lwd = 2
)


# -----------------------------------------------------------
# FINITE SAMPLE DIAGONAL APPROACH
# -----------------------------------------------------------
# Assume 'mcmc_clean' exists from your previous run
post_mat <- as.matrix(mcmc_clean)
n_sims <- nrow(post_mat)

# Define the Grid (Raw Scale)
y_grid <- seq(min(X), max(X) * 6, length.out = 200)
G_curves <- matrix(0, nrow = n_sims, ncol = length(y_grid))

# -----------------------------------------------------------
# 2. DEFINE THE UNIVERSAL GENERATORS
# -----------------------------------------------------------
# You can switch "FAMILY" to "Gumbel", "Clayton", or "Frank"

FAMILY <- "Gumbel" 

get_generator_functions <- function(fam) {
  if(fam == "Gumbel") {
    list(
      # psi(t) = exp(-t^(1/th))
      psi = function(t, th) exp(-t^(1/th)),
      # psi^-1(u) = (-log u)^th
      psi_inv = function(u, th) (-log(u))^th
    )
  } else if(fam == "Clayton") {
    list(
      # psi(t) = (1+t)^(-1/th)
      psi = function(t, th) (1+t)^(-1/th),
      # psi^-1(u) = u^(-th) - 1
      psi_inv = function(u, th) u^(-th) - 1
    )
  } else if(fam == "Frank") {
    list(
      # Frank is numerically tricky, use copula package or robust form
      psi = function(t, th) -log(1 - (1-exp(-th))*exp(-t))/th,
      psi_inv = function(u, th) -log((exp(-th*u)-1)/(exp(-th)-1))
    )
  }
}

GEN <- get_generator_functions(FAMILY)

# -----------------------------------------------------------
# 3. EXACT COMPUTATION LOOP (No approximations)
# -----------------------------------------------------------

for(i in 1:n_sims) {
  # A. Get Parameters
  theta <- post_mat[i, "theta"]
  mu    <- post_mat[i, "mu"]
  sigma <- post_mat[i, "sigma"]
  alpha <- post_mat[i, "alpha"] # Or 'nu' if Student-t
  
  # B. Calculate Exact Marginal CDF F(y)
  # Change this line depending on your chosen margin (Frechet or Student-t)
  
  # -- Option 1: Frechet --
  z <- (y_grid - mu) / sigma
  F_y <- numeric(length(y_grid))
  valid <- z > 0
  F_y[valid] <- exp(-(z[valid]^(-alpha)))
  
  # -- Option 2: Student-t (Uncomment if using t) --
  # F_y <- pt((y_grid - mu)/sigma, df=alpha) 
  
  # Safety Clamp (prevent log(0) or Inf)
  F_y <- pmin(pmax(F_y, 1e-10), 1 - 1e-10)
  
  # C. Apply Universal Diagonal Formula
  # G(y) = psi( n * psi_inv( F(y) ) )
  
  t_val <- GEN$psi_inv(F_y, theta) # Map to generator space
  t_sum <- n * t_val               # Sum dependence (Exchangeable)
  G_y   <- GEN$psi(t_sum, theta)   # Map back to probability
  
  G_curves[i, ] <- G_y
}

# ==============================================================================
# 7. VALIDATION: CALCULATE TRUE G (Using True Parameters)
# ==============================================================================

# 1. Calculate True Marginal CDF: F_true(y)
z_true <- (y_grid - mu_true) / sigma_true

F_true <- numeric(length(y_grid))

# Fréchet CDF Logic
valid_true <- z_true > 0
F_true[valid_true] <- pfrechet(z_true, shape = alpha_true)

# Safety clamp for numerical stability in the generator
F_true <- pmin(pmax(F_true, 1e-15), 1 - 1e-15)

# 2. Calculate True G(y) using the Universal Diagonal Formula
# G_true = psi( n * psi_inv( F_true ) ) based on TRUE theta

# Map F_true to generator space using theta_true
t_val_true <- GEN$psi_inv(F_true, theta_true)

# Apply Exchangeability (Sum n times)
t_sum_true <- n * t_val_true

# Map back to Probability space
G_true <- GEN$psi(t_sum_true, theta_true)

# ==============================================================================
# 8. VISUALIZATION
# ==============================================================================

par(mfrow = c(1, 1))

# 1. Plot the Posterior Mean (Blue)
plot(y_grid, G_mean, 
     type = "l", lwd = 2, col = "blue", 
     ylim = c(0, 1), 
     main = paste("Posterior vs True Maxima Distribution (", FAMILY, ")", sep=""),
     xlab = "Maximum Value (y)", 
     ylab = "P(Mn <= y)")

# 2. Add Credible Intervals (Shaded Blue)
polygon(c(y_grid, rev(y_grid)), c(G_lower, rev(G_upper)), 
        col = rgb(0, 0, 1, 0.1), border = NA)

# 3. Add the True Theoretical G (Red Dashed)
lines(y_grid, G_true, col = "red", lwd = 2, lty = 2)

legend("bottomright", 
       legend = c("Posterior Mean", "95% Credible Interval", "True G (Theoretical)"),
       col = c("blue", "blue", "red"), 
       lty = c(1, 1, 2), 
       lwd = c(2, 10, 2),
       fill = c(NA, rgb(0,0,1,0.1), NA), 
       border = NA)
