# Libraries for copula and tail dependence analysis
# install.packages(c("copula", "mev"))
library(copula)
library(mev)

# --- 1. BASE COPULA DEFINITIONS ---

# Tail Copula: Gaussian (Asymptotically Independent)
d_tail_cop <- function(u_star, v_star, rho) {
  cop <- normalCopula(param = rho, dim = 2)
  return(dCopula(cbind(u_star, v_star), cop))
}

# Bulk Copula: Gumbel (Asymptotically Dependent)
d_bulk_cop <- function(u_star, v_star, alpha) {
  cop <- gumbelCopula(param = alpha, dim = 2)
  return(dCopula(cbind(u_star, v_star), cop))
}

# Dynamic Weighting Function: pi(u*, v*) = (u* * v*)^theta
pi_weight <- function(u_star, v_star, theta) {
  return((u_star * v_star)^theta)
}

# --- 2. THE LATENT MIXTURE DENSITY (c*) ---
# This defines the joint density in the non-uniform latent space u*.
# K is the normalization constant required to make this a valid density.
c_star_unnorm <- function(u_star, v_star, theta, rho, alpha) {
  w <- pi_weight(u_star, v_star, theta)
  # Blending: Tail copula dominates at high values of u* and v*
  dens <- w * d_tail_cop(u_star, v_star, rho) + (1 - w) * d_bulk_cop(u_star, v_star, alpha)
  return(dens)
}

# --- 3. MODEL PARAMETERS ---
theta_param <- 1.5   # Transition weight parameter
rho_param   <- 0.5   # Gaussian parameter (Expected Limit: eta = 0.75, chi = 0)
alpha_param <- 1.2   # Gumbel parameter (Bulk behavior: chi = 0.58)

# --- 4. GLOBAL NORMALIZATION & MAPPING (F_u*) ---
cat("Calculating Global Normalization Constant K and Marginal Splines...\n")

# Compute K: The volume of the unnormalized mixture over [0,1]^2
K_val <- integrate(function(v) {
  sapply(v, function(vi) {
    integrate(function(ui) c_star_unnorm(ui, vi, theta_param, rho_param, alpha_param), 0, 1)$value
  })
}, 0, 1)$value

# Marginal Latent Density f_u*(s): normalized by K
f_u_star <- function(s) {
  val <- integrate(function(v) c_star_unnorm(s, v, theta_param, rho_param, alpha_param), 0, 1)$value
  return(val / K_val)
}

# Pre-compute Marginal CDF F_u*(s) for the mapping: u* -> u
grid_size <- 50
u_star_grid <- seq(0.001, 0.999, length.out = grid_size)
F_grid <- sapply(u_star_grid, function(limit) {
  integrate(function(s) sapply(s, f_u_star), 0, limit)$value
})

# Create Splines for fast space-switching
# F_u_star_spline: Latent u* -> Uniform u
F_u_star_spline <- splinefun(u_star_grid, F_grid, method = "monoH.FC")
# F_u_star_inv_spline: Uniform u -> Latent u*
F_u_star_inv_spline <- splinefun(F_grid, u_star_grid, method = "monoH.FC")

# --- 5. THE TRANSITION KERNEL ---

# Conditional CDF in the latent space: P(v* <= target | u* = current)
# Note: K cancels out here because it's in both numerator and denominator.
cond_cdf_latent <- function(target_v_star, current_u_star) {
  num <- integrate(function(v) c_star_unnorm(current_u_star, v, theta_param, rho_param, alpha_param), 0, target_v_star)$value
  den <- f_u_star(current_u_star) * K_val
  return(num / den)
}

# Function to propagate the chain one step: u_t -> u_{t+1}
sample_next_step <- function(u_t) {
  # 1. Map current uniform state to latent space: u_t -> u*_t
  u_star_t <- F_u_star_inv_spline(u_t)
  
  # 2. Sample next state in latent space via inverse transform sampling
  w <- runif(1)
  u_star_next <- uniroot(function(v) cond_cdf_latent(v, u_star_t) - w, 
                         lower = 0, upper = 1, tol = 1e-6)$root
  
  # 3. Map back to uniform space to ensure stationarity: u*_next -> u_{t+1}
  u_next <- F_u_star_spline(u_star_next)
  return(u_next)
}

# --- 6. SIMULATION ---
n_steps <- 200
u_chain <- numeric(n_steps)
u_chain[1] <- runif(1) # Start with a random uniform value

cat("Simulating Univariate Markov Chain...\n")
for(t in 1:(n_steps-1)) {
  u_chain[t+1] <- sample_next_step(u_chain[t])
  if(t %% 100 == 0) cat("Step", t, "completed...\n")
}

x_chain <- egpd::qegpd(u_chain, kappa = 2, sigma = 1, xi = 0.3)

# --- 7. VALIDATION ---
# Create consecutive pairs to check dependence
u_pairs <- cbind(u_chain[1:(n_steps-1)], u_chain[2:n_steps])

par(mfrow=c(1,2))
plot(u_chain, type="l", col="steelblue", main="Markov Chain Trace (u)", ylab="u_t")
plot(u_pairs, pch=20, col=rgb(0,0,0,0.2), main="Transition Copula (u_t, u_{t+1})",
     xlab="u_t", ylab="u_{t+1}")

# Check Tail Dependence (using mev package)
cat("\n--- Tail Dependence Analysis ---\n")
td <- mev::taildep(u_pairs, u = seq(0.7, 0.98, by = 0.01))

# Theoretical Values based on Gaussian Tail (rho=0.5) and Gumbel Bulk (alpha=2)
cat("Theoretical Tail Limit (r -> 1):\n")
cat("Expected chi: 0.00\n")
cat("Expected eta: ", (1 + rho_param)/2, "\n")