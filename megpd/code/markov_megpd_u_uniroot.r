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
  log_ratio <- log(x1) - log(x2)
  # Stability check for the exponential
  z_sq <- (log_ratio / d)^2
  term_ang <- (1 / (sqrt(2 * pi) * d)) * exp(-0.5 * z_sq) * (r / (x1 * x2))

  term_rad * term_ang
}

dmegpd_biv <- function(
    x1,
    x2,
    kappa,
    sigma,
    xi,
    delta_func,
    min_delta = 1e-6
) {

    # ---------- vector handling ----------
    n <- max(length(x1), length(x2))

    if (length(x1) == 1)
        x1 <- rep(x1, n)

    if (length(x2) == 1)
        x2 <- rep(x2, n)

    out <- numeric(n)

    # ---------- support ----------
    valid <-
        is.finite(x1) &
        is.finite(x2) &
        (x1 > 0) &
        (x2 > 0)

    if (!any(valid))
        return(out)

    x1v <- x1[valid]
    x2v <- x2[valid]

    # ---------- radial component ----------
    r <- x1v + x2v

    term_rad <- egpd::degpd_density(
        r,
        kappa = kappa,
        sigma = sigma,
        xi = xi
    )

    # if radial density itself misbehaves
    good_rad <- is.finite(term_rad) & (term_rad > 0)

    if (!any(good_rad))
        return(out)

    idx <- which(valid)[good_rad]

    x1v <- x1v[good_rad]
    x2v <- x2v[good_rad]
    r <- r[good_rad]
    term_rad <- term_rad[good_rad]

    # ---------- dependence parameter ----------
    d <- pmax(delta_func(r), min_delta)

    # ---------- log angular density ----------
    log_ratio <- log(x1v) - log(x2v)

    z <- log_ratio / d

    log_term_ang <-
        -0.5 * z^2 -
        log(sqrt(2 * pi)) -
        log(d) +
        log(r) -
        log(x1v) -
        log(x2v)

    # ---------- combine on log scale ----------
    log_density <- log(term_rad) + log_term_ang

    dens <- exp(log_density)

    dens[!is.finite(dens)] <- 0

    out[idx] <- dens

    out
}

u_to_x <- function(u) {
  u / (1 - u)
}

x_to_u <- function(x) {
  x / (x + 1)
}

# 1. Define the conditional density f(x | x_prev)
conditional_pdf_u <- function(u,
                              x_prev,
                              kappa,
                              sigma,
                              xi,
                              delta_func) {

  jacobian <- 1 / (1 - u)^2

  f_x <- dmegpd_biv(
          u_to_x(u),
          x_prev,
          kappa,
          sigma,
          xi,
          delta_func
      )

  f_x * jacobian
}

# 2. Function to compute the CDF at a point 'target_x' using adaptive integration
get_cdf_u_val <- function(target_u, x_prev, kappa, sigma, xi, delta_func, norm_const) {

  # Integrate from 0 to target_x
  # Use a log-transform internally if target_x is very large
  val <- integrate(conditional_pdf_u, lower = 0, upper = target_u, x_prev=x_prev, kappa=kappa, sigma=sigma, xi=xi, 
                   delta_func=delta_func)$value
  return(val / norm_const)
}

# 3. Sampling via Uniroot
sample_conditional <- function(x_prev, kappa, sigma, xi, delta_func) {
  
  norm_const <- integrate(conditional_pdf_u, lower = 0, upper = 1, x_prev=x_prev, kappa=kappa, sigma=sigma, xi=xi, delta_func=delta_func)$value
  
  # Step B: Draw Uniform
  p_target <- runif(1)
  
  # Step C: Find x such that CDF(x) = p_target
  # search range: [lower_bound, upper_bound]
  res <- uniroot(function(u) {
    get_cdf_u_val(u, x_prev, kappa, sigma, xi, delta_func, norm_const) - p_target
  }, interval = c(1e-10, 1-1e-10), extendInt = "upX")
  
  u_to_x(res$root)
}

simulate_megpd_chain <- function(
  n_steps,
  kappa,
  sigma,
  xi,
  delta_func,
  x0 = 1,
  burn_in_prop = 0.10,
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
    x[t] <- sample_conditional(
      x_prev = x[t - 1],
      kappa = kappa,
      sigma = sigma,
      xi = xi,
      delta_func = delta_func
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