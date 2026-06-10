library(evd)
library(tidyverse)
library(copula)

###########################################################
# 1. ROBUST CORE FUNCTIONS
###########################################################

dmegpd_biv <- function(x1, x2, kappa, xi, delta_func) {
  # Numerical stability: prevent non-positive values
  x1 <- max(x1, 1e-10)
  x2 <- max(x2, 1e-10)
  
  r <- x1 + x2
  # Radial components (unit scale)
  H_r <- if(abs(xi) < 1e-10) 1 - exp(-r) else pmax(0, 1 - (1 + xi * r)^(-1/xi))
  h_r <- if(abs(xi) < 1e-10) exp(-r) else pmax(0, (1 + xi * r)^(-1/xi - 1))
  
  # Ensure delta is not too small to cause overflow
  d <- max(delta_func(r), 0.01)
  
  # Term 1: Radial part
  term_rad <- kappa * h_r * (H_r^(kappa - 1))
  
  # Term 2: Angular part (Gaussian on log-ratios)
  log_ratio <- log(x1 / x2)
  # Stability check for the exponential
  z_sq <- (log_ratio / d)^2
  term_ang <- (1 / (sqrt(2 * pi) * d)) * exp(-0.5 * z_sq) * (r / (x1 * x2))
  
  val <- term_rad * term_ang
  return(if(is.finite(val)) val else 0)
}


##################################################
# Build grid
##################################################

make_log_grid <- function(x_prev,
                      n_grid = 1000,
                      min_ratio = 10,
                      max_ratio = 3,
                      global_min = 1e-5,
                      global_max = 10) {

  grid_min <- min(global_min, x_prev / min_ratio)
  grid_max <- max(global_max, x_prev * max_ratio)

  exp(seq(log(grid_min),
          log(grid_max),
          length.out = n_grid))
}

##################################################
# Evaluate conditional density
##################################################

conditional_pdf <- function(x_prev,
                            x_grid,
                            kappa,
                            xi,
                            delta_func) {

  sapply(
    x_grid,
    function(x)
      dmegpd_biv(
        x_prev,
        x,
        kappa,
        xi,
        delta_func
      )
  )
}

##################################################
# Numerically integrate density -> CDF
##################################################

pdf_to_cdf <- function(x_grid,
                       pdf_vals) {

  dx <- diff(x_grid)

  # trapezoidal rule for CDF
  mid_pdf <- (
    pdf_vals[-1] +
    pdf_vals[-length(pdf_vals)]
  ) / 2

  cdf <- c(
    0,
    cumsum(mid_pdf * dx)
  )

  total_mass <- tail(cdf, 1)

  if(total_mass <= 0)
    return(NULL)

  cdf <- cdf / total_mass

  # numerical safeguard
  cdf <- cummax(cdf)

  # Numerical Safeguard: Splines require unique x-values. 
  # In the flat tails, CDF values might repeat. We add a tiny jitter.
  if(any(diff(cdf) == 0)) {
    cdf <- cdf + seq(0, 1e-12, length.out = length(cdf))
  }

  list(
    cdf = cdf,
    total_mass = total_mass
  )
}

##################################################
# Inverse-CDF sampler
##################################################

sample_from_cdf <- function(x_grid,
                            cdf_vals,
                            u = runif(1)) {

  # approx(
  #   x = cdf_vals,
  #   y = x_grid,
  #   xout = u,
  #   ties = "ordered"
  # )$y

  F_inv <- splinefun(cdf_vals, x_grid, method = "monoH.FC")
  F_inv(u)
}

# ROBUST INTEGRATION
# Instead of 0 to Inf, we integrate in log-space to handle heavy tails and singularities near 0
dmegpd_marginal_integrated <- function(x_target, kappa, xi, delta_func) {
  
  # Integrand: f(x1, x2) dx2
  # We integrate over v = log(x2)
  integrand_log <- function(log_x2) {
    x2 <- exp(log_x2)
    sapply(x2, function(v) dmegpd_biv(x_target, v, kappa, xi, delta_func) * v) # *v is Jacobian
  }
  
  # Range: from log(1e-8) to log(1000) covers the bulk and tail effectively
  res <- try(integrate(integrand_log, lower = -18, upper = 7), silent=TRUE)
  
  if(inherits(res, "try-error")) return(0)
  return(res$value)
}

###########################################################
# 2. UPDATED SAMPLER
###########################################################

##################################################
# Markov transition sampler
##################################################

markov_megpd_sampler <- function(x_prev,
                                 kappa,
                                 xi,
                                 delta_func,
                                 n_grid = 1000) {

  x_grid <- make_log_grid(
    x_prev,
    n_grid = n_grid
  )

  pdf_vals <- conditional_pdf(
    x_prev,
    x_grid,
    kappa,
    xi,
    delta_func
  )

  if(sum(pdf_vals) <= 0)
    return(x_prev)

  cdf_obj <- pdf_to_cdf(
    x_grid,
    pdf_vals
  )

  if(is.null(cdf_obj)) {
    print("Warning: CDF construction failed, returning previous value.")
    return(x_prev)
  }

  sample_from_cdf(
    x_grid,
    cdf_obj$cdf
  )
}

###########################################################
# 3. SIMULATION
###########################################################

set.seed(1)
n_steps <- 2000
kappa_val <- 2
xi_val <- 0.25

# Your custom delta function (Non-monotonic)
delta_f <- function(r) {
  base <- 0.6
  valley <- dgamma(r, shape = 6, scale = 4)
  valley <- valley / max(valley)
  pmax(0.05, base - 0.2 * valley)
}

# Dependence that gets STRONGER as values get LARGER
delta_f_strong_upper <- function(r) {
  # Starts at 0.8 (weak dependence at 0) and decays to 0.1 (strong dependence at infinity)
  0.1 + 0.7 * exp(-r / 5)
}

# Run
x <- numeric(n_steps)
x[1] <- 1.0
pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
for (t in 2:n_steps) {
  x[t] <- markov_megpd_sampler(x[t-1], kappa_val, xi_val, delta_f_strong_upper)
  setTxtProgressBar(pb, t)
}
close(pb)

final_chain <- x[101:n_steps]

# Prepare bivariate pairs from the Markov Chain
x_tm1 <- final_chain[-length(final_chain)]
x_t   <- final_chain[-1]

# Transform to Uniform marginals (Empirical PIT)
u_tm1 <- copula::pobs(x_tm1)
u_t   <- copula::pobs(x_t)
