library(evd)
library(tidyverse)
library(copula)

###########################################################
# 1. ROBUST CORE FUNCTIONS
###########################################################

dmegpd_biv <- function(x1, x2, kappa, sigma, xi, delta_func) {

  r <- x1 + x2

  # Radial components (unit scale)
  term_rad <- egpd::degpd_density(r,
    kappa = kappa, sigma = sigma,
    xi = xi
  )

  # Ensure delta is not too small to cause overflow
  d <- max(delta_func(r), 0.01)

  # Term 2: Angular part (Gaussian on log-ratios)
  log_ratio <- log(x1 / x2)
  # Stability check for the exponential
  z_sq <- (log_ratio / d)^2
  term_ang <- (1 / (sqrt(2 * pi) * d)) * exp(-0.5 * z_sq) * (r / (x1 * x2))

  term_rad * term_ang
}

##################################################
# Build grid
##################################################

make_grid <- function(x_prev, n_grid = 5000) {
  seq(1e-10, x_prev * 100, length.out = n_grid)
}

##################################################
# Evaluate conditional density
##################################################

conditional_pdf <- function(x_prev,
                            x_grid,
                            kappa,
                            sigma,
                            xi,
                            delta_func) {
  sapply(
    x_grid,
    function(x) {
      dmegpd_biv(
        x_prev,
        x,
        kappa,
        sigma,
        xi,
        delta_func
      )
    }
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

  if (total_mass <= 0) {
    return(NULL)
  }

  cdf <- cdf / total_mass

  # numerical safeguard
  cdf <- cummax(cdf)

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

  keep <- c(TRUE, diff(cdf_vals) > 1e-12)

  F_inv <- splinefun(
    cdf_vals[keep],
    x_grid[keep],
    method = "monoH.FC"
  )
  F_inv(u)
}

###########################################################
# 2. UPDATED SAMPLER
###########################################################

##################################################
# Markov transition sampler
##################################################

markov_megpd_sampler <- function(x_prev,
                                 kappa,
                                 sigma,
                                 xi,
                                 delta_func,
                                 n_grid = 5000) {
  x_grid <- make_grid(
    x_prev,
    kappa = kappa,
    sigma = sigma,
    xi = xi,
    n_grid = n_grid
  )

  pdf_vals <- conditional_pdf(
    x_prev,
    x_grid,
    kappa,
    sigma,
    xi,
    delta_func
  )

  if (sum(pdf_vals) <= 0) {
    return(x_prev)
  }

  cdf_obj <- pdf_to_cdf(
    x_grid,
    pdf_vals
  )

  cat("\n Total integrated mass:", cdf_obj$total_mass, "\n")

  if (is.null(cdf_obj)) {
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

###########################################################
# Markov chain simulation
###########################################################

simulate_megpd_chain <- function(
  n_steps = 2000,
  kappa,
  sigma,
  xi,
  delta_func,
  x0 = 1,
  burn_in_prop = 0.10,
  n_grid = 5000,
  show_progress = TRUE
) {
  if (burn_in_prop < 0 || burn_in_prop >= 1) {
    stop("burn_in_prop must be in [0, 1).")
  }

  x <- numeric(n_steps)
  x[1] <- x0

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
  }

  for (t in 2:n_steps) {
    x[t] <- markov_megpd_sampler(
      x_prev = x[t - 1],
      kappa = kappa,
      sigma = sigma,
      xi = xi,
      delta_func = delta_func,
      n_grid = n_grid
    )

    if (show_progress) {
      setTxtProgressBar(pb, t)
    }
  }

  if (show_progress) {
    close(pb)
  }

  burn_in <- floor(n_steps * burn_in_prop)

  list(
    full_chain = x,
    final_chain = x[(burn_in + 1):n_steps],
    burn_in = burn_in
  )
}

set.seed(1)
n <- 2000
kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

# Dependence that gets STRONGER as values get LARGER
delta_strong_upper <- function(r) {
  # Starts at 0.8 (weak dependence at 0) and decays to 0.1 (strong dependence at infinity)
  0.2 + 0.6 * exp(-r / 5)
}
    
# Run
set.seed(1)

sim <- simulate_megpd_chain(
  n_steps = n,
  kappa = kappa_val,
  sigma = sigma_val,
  xi = xi_val,
  delta_func = delta_strong_upper,
  x0 = 1,
  burn_in_prop = 0.10,
  n_grid = 100
)

final_chain <- sim$final_chain
