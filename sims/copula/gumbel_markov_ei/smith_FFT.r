rm(list=ls())

# -----------------------------------------------------------------------------
# 1. Model Definitions (Logistic/Gumbel)
# -----------------------------------------------------------------------------
# Density of the increments h(z)
h_logistic <- function(z, r) {
  # h(z) = (r - 1) * exp(-r * z) * (1 + exp(-r * z))^(1/r - 2)
  val <- (r - 1) * exp(-r * z) * (1 + exp(-r * z))^(1/r - 2)
  return(val)
}

# CDF of increments H(z) (for tail correction)
H_cdf_logistic <- function(z, r) {
  return((1 + exp(-r * z))^(1/r - 1))
}

# -----------------------------------------------------------------------------
# 2. FFT Solver (The "Smith Method")
# -----------------------------------------------------------------------------
solve_extremal_index_fft <- function(alpha, N_points = 2^14, x_max = 80) {
  
  # A. Setup Lattice
  # Smith used up to 2^14 points.
  dx <- x_max / (N_points - 1)
  x_grid <- seq(0, x_max, length.out = N_points)
  
  # B. Prepare Kernels for FFT
  # We need to evaluate h(x - y). 
  # For the convolution Q(y)h(x-y), if x and y are on [0, L], 
  # the argument (x-y) ranges from -L to L.
  # We construct a kernel vector that represents h(z) centered correctly for FFT.
  
  # Construct the kernel vector on the symmetric range roughly [-L, L]
  # In R's fft(), index 1 is 0, index 2 is dx, ... 
  # Negative frequencies wrap to the end of the vector.
  
  # We evaluate h at 0, dx, 2dx... and also -dx, -2dx (wrapped to end)
  h_vec <- numeric(N_points)
  
  # Positive lag part (0 to x_max)
  h_vec[1:N_points] <- h_logistic(x_grid, alpha) * dx
  
  # We need a padded vector to perform Linear Convolution (avoid circular wrap artifacts)
  # Length 2*N is standard for linear convolution via FFT.
  N_fft <- 2 * N_points
  
  # Create padded kernel 'h_pad'
  # Organization for R FFT: [0, dx, ..., xmax, 0, ..., 0, -xmax, ..., -dx]
  # But simpler approach: standard linear convolution layout.
  # Let's use the standard "convolve" logic:
  # Signal A: Q (padded with zeros)
  # Signal B: h (evaluated on full relevant range)
  
  # Let's simplify: We want to compute integral_{0}^{x_max} Q(y) h(x-y) dy
  # We can construct a toeplitz-like operation using FFT.
  
  # Construct h on the range [-x_max, x_max] to capture all shifts
  # However, for the iterative update Q(x), we only need x in [0, x_max].
  # x - y ranges from -x_max (when x=0, y=xmax) to x_max (when x=xmax, y=0).
  
  h_vals_full <- h_logistic(seq(-x_max, x_max, by=dx), alpha) * dx
  # This vector has length 2*N_points - 1
  
  # Tail correction vector (Constant for a given alpha/grid)
  # C(x) = H(x - x_max)
  tail_correction <- H_cdf_logistic(x_grid - x_max, alpha)
  
  # C. Iteration
  # Initialize Q = 1 on [0, x_max]
  Q_curr <- rep(1, N_points)
  
  iter <- 0
  diff <- 1
  tol <- 1e-8
  
  # We use the 'convolve' function which uses FFT internally, 
  # or we can do manual FFT pad. 'convolve' with type="open" is easiest 
  # to understand: it gives the full linear convolution.
  
  while(diff > tol && iter < 1000) {
    
    # 1. Convolution: Integral Q(y) h(x-y) dy
    # This effectively "smears" Q by the shape of h.
    # We want the result at indices corresponding to x in [0, x_max].
    # Using base R convolve (uses FFT). 
    # Note: convolve(x, y) in R reverses the second sequence by default for convolution definition!
    # So we pass h_vals_full as is, but we must be careful with alignment.
    
    # To be precise and fast, let's do manual FFT multiplication:
    # Pad Q to length 4N (safe bet)
    fft_len <- 4 * N_points
    Q_pad <- c(Q_curr, numeric(fft_len - N_points))
    
    # Pad h. h must be arranged so that lag 0 is at index 1.
    # We need h(z) for z in [-x_max, x_max].
    # Construct "wrapped" kernel for FFT:
    # [ h(0), h(dx), ..., h(xmax), 0...0, h(-xmax), ..., h(-dx) ]
    h_kernel <- numeric(fft_len)
    
    # 0 to x_max
    h_kernel[1:N_points] <- h_logistic(x_grid, alpha) * dx
    # -x_max to -dx (put at end of vector)
    # x_grid[2] is dx. x_grid[N] is x_max.
    neg_z <- -x_grid[2:N_points] # -dx ... -xmax
    neg_vals <- h_logistic(neg_z, alpha) * dx
    h_kernel[(fft_len - length(neg_vals) + 1):fft_len] <- rev(neg_vals) # -xmax ... -dx
    
    # Perform Convolution in Frequency Domain
    # Q_out = IFFT( FFT(Q) * FFT(h) )
    Q_conv_full <- Re(fft(fft(Q_pad) * fft(h_kernel), inverse = TRUE)) / fft_len
    
    # Extract the valid part [0, x_max]
    # Since Q starts at index 1 (y=0) and h_kernel has 0 lag at index 1,
    # the result at index k corresponds to x = (k-1)*dx.
    Q_integral <- Q_conv_full[1:N_points]
    
    # 2. Add Tail Correction
    Q_next <- Q_integral + tail_correction
    
    # Check convergence
    diff <- max(abs(Q_next - Q_curr))
    Q_curr <- Q_next
    iter <- iter + 1
  }
  
  cat(sprintf("Alpha %.1f: Converged in %d iterations (FFT).\n", alpha, iter))
  
  # D. Final Integration for Theta
  # Theta = Integral_{-inf}^{0} exp(x) * Q(x) dx
  # We need Q(x) for x < 0.
  # Q(x) = Integral_{0}^{x_max} Q(y) h(x-y) dy + H(x - x_max)
  
  # We need a grid for negative x.
  # Let's go from -30 to 0.
  x_neg_grid <- seq(-30, 0, by=dx)
  N_neg <- length(x_neg_grid)
  
  # We can compute the integral part for these negative x values using the same logic.
  # Or, since N_neg is small (~600 points) compared to N_points (16000), 
  # we can just do a direct calculation here (vectorized) to be safe and simple.
  # It's technically O(N_neg * N_points), which is fine for one step.
  
  Q_neg <- numeric(N_neg)
  
  # Vectorized calculation for negative domain
  for(i in 1:N_neg) {
    x_val <- x_neg_grid[i]
    # h(x - y_grid)
    z_vals <- x_val - x_grid
    
    # Integral part
    int_val <- sum(Q_curr * h_logistic(z_vals, alpha)) * dx
    # Correction part
    corr_val <- H_cdf_logistic(x_val - x_max, alpha)
    
    Q_neg[i] <- int_val + corr_val
  }
  
  # Integrate exp(x) * Q(x)
  integrand <- exp(x_neg_grid) * Q_neg
  theta <- sum((integrand[-1] + integrand[-length(integrand)])/2) * dx
  
  return(theta)
}

# -----------------------------------------------------------------------------
# 3. Execution
# -----------------------------------------------------------------------------

# Parameters from Smith's Table 1
alphas <- c(2.0, 3.0, 4.0, 5.0)
results_fft <- data.frame(Alpha = alphas, Theta_FFT = NA)

for(i in 1:length(alphas)) {
  results_fft$Theta_FFT[i] <- solve_extremal_index_fft(alphas[i])
}

print(results_fft)

# Comparison with Paper Table 1
cat("\n--- Comparison (Smith 1992, Table 1) ---\n")
cat("r=2.0 | Paper: 0.328 | FFT: ", round(results_fft$Theta_FFT[1], 4), "\n")
cat("r=3.0 | Paper: 0.158 | FFT: ", round(results_fft$Theta_FFT[2], 4), "\n")
cat("r=4.0 | Paper: 0.0925| FFT: ", round(results_fft$Theta_FFT[3], 4), "\n")
cat("r=5.0 | Paper: 0.0616| FFT: ", round(results_fft$Theta_FFT[4], 4), "\n")