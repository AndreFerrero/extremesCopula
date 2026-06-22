library(evd)
library(tidyverse)
library(copula)

###########################################################
# 1. ROBUST CORE FUNCTIONS
###########################################################

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

# 1. Define the conditional density f(x | x_prev)
conditional_pdf_x <- function(x, x_prev, kappa, sigma, xi, delta_func) {
  # We return the joint density; normalization happens later
  dmegpd_biv(x, x_prev, kappa, sigma, xi, delta_func)
}

# 2. Function to compute the CDF at a point 'target_x' using adaptive integration
get_cdf_val <- function(target_x, x_prev, kappa, sigma, xi, delta_func, norm_const) {
  # Integrate from 0 to target_x
  # Use a log-transform internally if target_x is very large
  val <- integrate(conditional_pdf_x, lower = 0, upper = target_x, 
                   x_prev=x_prev, kappa=kappa, sigma=sigma, xi=xi, 
                   delta_func=delta_func, subdivisions = 200)$value
  return(val / norm_const)
}

# 3. Sampling via Uniroot
sample_conditional <- function(x_prev, kappa, sigma, xi, delta_func) {
  # Step A: Find the normalization constant (Total Area)
  # Integrate over a massive range or use log-range based on x_prev
  d <- delta_func(x_prev)
  lower_bound <- exp(log(x_prev) - 8*d)
  upper_bound <- exp(log(x_prev) + 8*d)
  
  norm_const <- integrate(conditional_pdf_x, lower = 1e-10, upper = Inf, 
                          x_prev=x_prev, kappa=kappa, sigma=sigma, xi=xi, 
                          delta_func=delta_func, subdivisions = 200)$value
  
  # Step B: Draw Uniform
  p_target <- runif(1)
  
  # Step C: Find x such that CDF(x) = p_target
  # search range: [lower_bound, upper_bound]
  res <- uniroot(function(x) {
    get_cdf_val(x, x_prev, kappa, sigma, xi, delta_func, norm_const) - p_target
  }, interval = c(1e-8, upper_bound * 2), extendInt = "yes")
  
  return(res$root)
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

simulate_R_quantiles <- function(
    M,
    n_steps,
    probs = seq(0.001, 0.999, length.out = 400),
    kappa,
    sigma,
    xi,
    delta_func
) {

    Q <- matrix(
        NA_real_,
        nrow = M,
        ncol = length(probs)
    )

    success <- logical(M)
    error_message <- rep(NA_character_, M)

    for (m in seq_len(M)) {

        cat(
            "Replication", m, "/", M, "\n"
        )

        res <- tryCatch(

            {

                sim <- simulate_megpd_chain(
                    n_steps = n_steps,
                    kappa = kappa,
                    sigma = sigma,
                    xi = xi,
                    delta_func = delta_func,
                    show_progress = FALSE
                )

                x <- sim$final_chain

                Rt <- x[-1] + x[-length(x)]

                Q[m, ] <- quantile(
                    Rt,
                    probs = probs,
                    names = FALSE,
                    type = 8
                )

                success[m] <- TRUE

                NULL

            },

            error = function(e) {

                error_message[m] <<- conditionMessage(e)

                cat(
                    "FAILED:",
                    conditionMessage(e),
                    "\n"
                )

                NULL

            }

        )
    }

    list(
        probs = probs,
        Q = Q,
        success = success,
        error_message = error_message
    )
}

kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

delta_strong_upper <- function(r) {
  # Starts at 0.8 (weak dependence at 0) and decays to 0.1 (strong dependence at infinity)
  0.2 + 0.6 * exp(-r / 5)
}

sim_qq <- simulate_R_quantiles(
    M = 200,
    n_steps = 10000,
    probs = seq(0.001, 0.999, length.out = 400),
    kappa = kappa_val,
    sigma = sigma_val,
    xi = xi_val,
    delta_func = delta_strong_upper
)

save(sim_qq, file = "megpd/radius_simqq_M200_n10000_xi05.Rdata")

Q <- sim_qq$Q
p <- sim_qq$probs

theoretical <- egpd::qegpd(
    p,
    kappa = kappa_val,
    sigma = sigma_val,
    xi = xi_val
)

df <- tibble(
    p = p,
    theo = theoretical,
    lower = apply(Q, 2, quantile, 0.025, na.rm = TRUE),
    median = apply(Q, 2, quantile, 0.50, na.rm = TRUE),
    upper = apply(Q, 2, quantile, 0.975, na.rm = TRUE)
)

p_qq <- ggplot(df) +

    geom_ribbon(
        aes(
            x = theo,
            ymin = lower,
            ymax = upper
        ),
        alpha = 0.25
    ) +

    geom_line(
        aes(theo, median),
        linewidth = 1
    ) +

    geom_abline(
        slope = 1,
        intercept = 0,
        colour = "red",
        linetype = 2
    ) +

    scale_x_log10() +
    scale_y_log10() +

    labs(
        title = "Radius Quantiles distribution",
        x = "Theoretical EGPD quantile",
        y = "Empirical quantile"
    ) +

    theme_bw()

# ggsave("megpd/radius_quantiles_n10000_200mc_xi05.png", p_qq)

# load("C:/Users/Andrea Ferrero/extremesCopula/megpd/radius_simqq_M100_n10000_xi05.Rdata")
# p_qq
