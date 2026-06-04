# Load necessary libraries
library(copula) # For standard copula densities
library(stats)  # For spline interpolation

# EGPD CDF
pEGPD <- function(x, kappa, sigma, xi) {
  # Standard GPD part
  gpd <- 1 - (1 + xi * x / sigma)^(-1/xi)
  gpd[x <= 0] <- 0
  return(pmax(0, pmin(1, gpd^kappa)))
}

# EGPD Density (for the marginal component of the LL)
dEGPD <- function(x, kappa, sigma, xi) {
  g_x <- 1 - (1 + xi * x / sigma)^(-1/xi)
  # Basic GPD density
  f_gpd <- (1/sigma) * (1 + xi * x / sigma)^(-1/xi - 1)
  # Chain rule for EGPD
  dens <- kappa * (g_x^(kappa - 1)) * f_gpd
  return(pmax(1e-10, dens))
}

#' 1. The Unnormalized Density (c*)
#' This is the core mixture: pi*ct + (1-pi)*cb
unnormalized_mixture <- function(u, v, theta, alpha, rho) {

  u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
  v <- pmin(pmax(v, 1e-8), 1 - 1e-8)

  out <- tryCatch({

    weight <- (u * v)^theta

    ct <- dCopula(cbind(u, v), gumbelCopula(alpha))
    cb <- dCopula(cbind(u, v), normalCopula(rho))

    val <- weight * ct + (1 - weight) * cb

    val[!is.finite(val)] <- 1e10

    pmax(val, 1e-12)

  }, error = function(e) {
    rep(1e10, length(u))
  })

  return(out)
}

#' 2. The Correction Subroutine
#' This function performs the integration and builds the spline inverse.
get_copula_correction <- function(theta, alpha, rho, grid_size = 100) {
  
  # Define a grid for numerical integration
  grid <- seq(0.001, 0.999, length.out = grid_size)
  
  # Step A: Compute the Marginal Density fU*(u) 
  # We integrate out 'v' for every 'u' on our grid.
  # Since it's symmetric, fU* = fV*
  f_star_grid <- sapply(grid, function(u_val) {
    # Numerical integration over v [0,1]
    integrate(function(v) unnormalized_mixture(u_val, v, theta, alpha, rho),
              lower = 0, upper = 1)$value
  })

  # Step B: Compute the Normalizing Constant K
  # Total volume under the surface
  K <- integrate(function(u) {
    sapply(u, function(u_i) {
      integrate(function(v) unnormalized_mixture(u_i, v, theta, alpha, rho),
                lower=0, upper=1)$value
    })
  }, lower=0, upper=1)$value
  
  # Step C: Normalize the marginal density
  f_star_norm <- f_star_grid / K
  f_star_norm <- pmax(f_star_grid / K, 1e-12)

  # Step D: Compute the warped CDF FU*(u) via cumulative integration
  # Using simple trapezoidal rule for the grid
  F_star_grid <- cumsum(f_star_norm) * (grid[2] - grid[1])
  eps_seq <- seq_along(F_star_grid) * 1e-12
  F_star_grid <- F_star_grid + eps_seq
  F_star_grid <- F_star_grid / max(F_star_grid) # Force to end at 1.0
  
  # Step E: MONOTONIC SPLINE INVERSION
  # We want a function where input is Probability (p) and output is Warped Value (w)
  # So we swap axes: x = F_star_grid, y = grid
  inv_F_spline <- splinefun(x = F_star_grid, y = grid, method = "monoH.FC")
  
  # Also need a density interpolator for the denominator of Eq 6
  dens_f_spline <- splinefun(x = grid, y = f_star_norm, method = "monoH.FC")
  
  return(list(
    inv_F = inv_F_spline,
    dens_f = dens_f_spline,
    K = K
  ))
}

# Full Joint Likelihood Function
markov_egpd_joint_ll <- function(params, x_data) {
  
  # --- Step A: Extract Parameters ---
  # We often use log-transforms to ensure positivity during optimization
  kappa <- exp(params[1])
  sigma <- exp(params[2])
  xi    <- params[3]       # xi can be negative
  theta <- exp(params[4])
  alpha <- 1 + exp(params[5]) # Gumbel alpha must be > 1
  rho   <- tanh(params[6])    # Gaussian rho must be in (-1, 1)
  
  n <- length(x_data)
  
  # --- Step B: Marginal Component ---
  # 1. Transform raw data X to Uniform U using EGPD CDF
  u_all <- pEGPD(x_data, kappa, sigma, xi)
  
  # 2. Calculate log-marginal density for all points
  log_dens_margin <- sum(log(dEGPD(x_data, kappa, sigma, xi)))
  
  # --- Step C: Copula Component (The Correction Routine) ---
  # Transition pairs: u_t and u_{t+1}
  u_t    <- u_all[1:(n-1)]
  u_next <- u_all[2:n]
  
  # Use the subroutine you provided
  # Note: get_copula_correction is called ONCE per likelihood evaluation
  # but it is expensive.
  corr <- get_copula_correction(theta, alpha, rho)
  
  # Map to warped space with your clamping
  w_u <- pmin(pmax(corr$inv_F(u_t), 1e-6), 1 - 1e-6)
  w_next <- pmin(pmax(corr$inv_F(u_next), 1e-6), 1 - 1e-6)
  
  # Equation 6 calculation
  c_star_val <- unnormalized_mixture(w_u, w_next, theta, alpha, rho) / corr$K
  denom <- corr$dens_f(w_u) * corr$dens_f(w_next)
  
  # Avoid log(0)
  final_c <- pmax(1e-10, c_star_val / denom)
  log_dens_copula <- sum(log(final_c))
  
  # --- Step D: Total Log-Likelihood ---
  total_ll <- log_dens_margin + log_dens_copula
  
  # Return negative LL for optimization
  return(-total_ll)
}

# 1. Your time series data

# 2. Initial Guesses
# It is best to start with a 'warm start' from a simple GPD fit
init_params <- c(
  log_kappa = 0,     # kappa = 1
  log_sigma = log(mean(x_obs)), 
  xi = 0.1, 
  log_theta = log(0.8), 
  log_alpha = log(1.0), # (Gives alpha = 2.0)
  atanh_rho = 0.5       # rho = 0.46
)

# 3. Optimize
fit <- optim(
  par = init_params, 
  fn = markov_egpd_joint_ll, 
  x_data = x_obs,
  method = "L-BFGS-B",
  control = list(trace = 1) # Shows progress
)

# 4. Extract results
final_kappa <- exp(fit$par[1])
final_xi    <- fit$par[3]
final_theta <- exp(fit$par[4])
final_alpha <- 1 + exp(fit$par[5])
final_rho   <- tanh(fit$par[6])