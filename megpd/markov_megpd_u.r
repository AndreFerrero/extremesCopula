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

make_u_grid <- function(n_grid) {
  seq(1e-10, 1 - 1e-10, length.out = n_grid)
}

make_adaptive_grid <- function(u_prev, n_grid, phi = 200) {
  # 1. Global background (ensure we don't miss any part of the domain)
  g_uni <- seq(1e-10, 1-1e-10, length.out = n_grid * 0.2)
  
  # 2. Local concentration (zoom in on the current state)
  # Force Beta mode to match u_prev
  u_target <- max(0.001, min(0.999, u_prev))
  alpha <- u_target * (phi - 2) + 1
  beta  <- (1 - u_target) * (phi - 2) + 1
  
  g_beta <- qbeta(seq(1e-10, 1-1e-10, length.out = n_grid * 0.8), alpha, beta)
  
  # 3. Combine and remove duplicates
  grid <- sort(unique(c(g_uni, g_beta)))
  return(grid)
}

u_to_x <- function(u) {
  u / (1 - u)
}

x_to_u <- function(x) {
  x / (x + 1)
}

conditional_pdf_u <- function(x_prev,
                              u_grid,
                              kappa,
                              sigma,
                              xi,
                              delta_func) {

  x_grid <- u_to_x(u_grid)

  jacobian <- 1 / (1 - u_grid)^2

  f_x <- sapply(
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

  f_x * jacobian
}

pdf_to_cdf <- function(u_grid, pdf_vals) {

  du <- diff(u_grid)

  mid <- (pdf_vals[-1] + pdf_vals[-length(pdf_vals)]) / 2

  cdf <- c(0, cumsum(mid * du))

  total_mass <- tail(cdf, 1)

  if (total_mass <= 0) return(NULL)

  cdf <- cdf / total_mass
  cdf <- cummax(cdf)

  list(cdf = cdf, total_mass = total_mass)
}

sample_from_cdf <- function(u_grid, cdf_vals, u = runif(1)) {

  keep <- c(TRUE, diff(cdf_vals) > 1e-12)

  F_inv <- splinefun(
    cdf_vals[keep],
    u_grid[keep],
    method = "monoH.FC"
  )

  u_sample <- F_inv(u)

  u_sample
}

markov_megpd_sampler <- function(x_prev,
                                 kappa,
                                 sigma,
                                 xi,
                                 delta_func,
                                 n_grid) {


  u_grid <- make_adaptive_grid(u_prev = x_to_u(x_prev), n_grid = n_grid)
  # u_grid <- make_u_grid(n_grid = n_grid)

  pdf_vals <- conditional_pdf_u(
    x_prev = x_prev,
    u_grid = u_grid,
    kappa = kappa,
    sigma = sigma,
    xi = xi,
    delta_func = delta_func
  )

  if (sum(pdf_vals) <= 0) {
    return(x_prev)
  }

  cdf_obj <- pdf_to_cdf(u_grid, pdf_vals)

  if (is.null(cdf_obj)) {
    warning("CDF construction failed; returning previous value.")
    return(x_prev)
  }
  
  u_sample <- sample_from_cdf(
    u_grid,
    cdf_obj$cdf
  )

  # transform back to original scale
  u_to_x(u_sample)
}


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
n <- 100000
kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

# Dependence that gets STRONGER as values get LARGER
delta_strong_upper <- function(r) {
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
  n_grid = 1000
)

final_chain <- sim$final_chain

save(final_chain, file = "megpd/chain_100k_unit_integral_betagrid1k_xi05.Rdata")