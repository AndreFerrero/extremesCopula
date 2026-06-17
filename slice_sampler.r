# =============================================================================
# Univariate Slice Sampler
# =============================================================================
# Neal (2003) "Slice Sampling", Annals of Statistics
#
# Usage:
#   slice_sample(log_target, n, x0, w, m, lower, upper, ...)
#
# Arguments:
#   log_target : function returning the log (unnormalised) density at x.
#                Any extra parameters are passed via ...
#   n          : number of samples to draw
#   x0         : initial value
#   w          : initial bracket width (step-out increment); default 1
#   m          : max number of step-out doublings; Inf = unlimited; default Inf
#   lower      : hard lower bound of the support; default -Inf
#   upper      : hard upper bound of the support; default  Inf
#   ...        : additional arguments forwarded to log_target
#
# Returns a list with:
#   samples    : numeric vector of length n
#   n_evals    : total number of log_target evaluations
# =============================================================================

slice_sample <- function(log_target, n, x0,
                         w     = 1,
                         m     = Inf,
                         lower = -Inf,
                         upper =  Inf,
                         ...) {

  stopifnot(is.function(log_target))
  stopifnot(n > 0, w > 0)

  samples  <- numeric(n)
  n_evals  <- 0L
  x_cur    <- x0

  for (i in seq_len(n)) {

    # ------------------------------------------------------------------
    # 1. Draw the slice level  y ~ Uniform(0, f(x_cur))
    #    equivalently:         log_y = log_target(x_cur) - Exp(1)
    # ------------------------------------------------------------------
    lf_cur  <- log_target(x_cur, ...)
    n_evals <- n_evals + 1L
    log_y   <- lf_cur - rexp(1)          # log of the slice threshold

    # ------------------------------------------------------------------
    # 2. Step-out procedure to bracket the slice
    #    (Neal 2003, Figure 3 — "stepping out")
    # ------------------------------------------------------------------
    u  <- runif(1, 0, w)
    L  <- x_cur - u
    R  <- x_cur + (w - u)

    # Clamp to hard support
    L <- max(L, lower)
    R <- min(R, upper)

    if (is.finite(m)) {
      # Bounded step-out: at most m steps left or right
      J <- floor(runif(1, 0, m))         # max steps allowed leftward
      K <- (m - 1L) - J                  # max steps allowed rightward

      while (J > 0 && L > lower &&
             log_target(L, ...) > log_y) {
        L <- max(L - w, lower)
        J <- J - 1L
        n_evals <- n_evals + 1L
      }
      while (K > 0 && R < upper &&
             log_target(R, ...) > log_y) {
        R <- min(R + w, upper)
        K <- K - 1L
        n_evals <- n_evals + 1L
      }
    } else {
      # Unlimited step-out
      while (L > lower && log_target(L, ...) > log_y) {
        L <- max(L - w, lower)
        n_evals <- n_evals + 1L
      }
      while (R < upper && log_target(R, ...) > log_y) {
        R <- min(R + w, upper)
        n_evals <- n_evals + 1L
      }
    }

    # ------------------------------------------------------------------
    # 3. Shrinkage procedure: sample x_new from [L, R] ∩ slice
    #    (Neal 2003, Figure 5 — "shrinkage")
    # ------------------------------------------------------------------
    repeat {
      x_new   <- runif(1, L, R)
      lf_new  <- log_target(x_new, ...)
      n_evals <- n_evals + 1L

      if (lf_new >= log_y) break         # accept

      # Shrink the bracket
      if (x_new < x_cur) L <- x_new else R <- x_new
    }

    x_cur     <- x_new
    samples[i] <- x_cur
  }

  list(samples = samples, n_evals = n_evals)
}


# =============================================================================
# Convenience wrapper: also returns a simple summary and trace plot
# =============================================================================
run_slice_sampler <- function(log_target, n, x0,
                              w = 1, m = Inf,
                              lower = -Inf, upper = Inf,
                              burnin = floor(n / 5),
                              plot   = TRUE,
                              ...) {

  cat(sprintf("Running slice sampler: %d iterations (burnin = %d)\n", n, burnin))
  t0  <- proc.time()
  out <- slice_sample(log_target, n, x0, w = w, m = m,
                      lower = lower, upper = upper, ...)
  elapsed <- (proc.time() - t0)["elapsed"]

  post  <- out$samples[seq(burnin + 1L, n)]
  cat(sprintf("  log-target evals : %d\n", out$n_evals))
  cat(sprintf("  elapsed          : %.3f s\n", elapsed))
  cat(sprintf("  post-burnin n    : %d\n", length(post)))
  cat(sprintf("  mean / sd        : %.4f / %.4f\n", mean(post), sd(post)))
  cat(sprintf("  2.5%% / 97.5%%  : %.4f / %.4f\n",
              quantile(post, 0.025), quantile(post, 0.975)))

  if (plot) {
    old_par <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    on.exit(par(old_par))

    plot(out$samples, type = "l", col = "steelblue",
         xlab = "Iteration", ylab = "x",
         main = "Trace (full chain)")
    abline(v = burnin, col = "firebrick", lty = 2)

    hist(post, breaks = 50, freq = FALSE,
         col = "steelblue", border = "white",
         xlab = "x", main = "Posterior histogram")

    acf(post, main = "ACF (post-burnin)", col = "steelblue")
  }

  invisible(list(samples = post, full_chain = out$samples,
                 n_evals = out$n_evals))
}


# =============================================================================
# Examples
# =============================================================================
if (FALSE) {

  ## --- 1. Standard Normal (log-target up to a constant) --------------------
  log_norm <- function(x) -0.5 * x^2

  res1 <- run_slice_sampler(log_norm, n = 5000, x0 = 0, w = 2)


  ## --- 2. Beta(2, 5) on (0, 1) — pass shape params via ... ----------------
  log_beta <- function(x, a, b) {
    if (x <= 0 || x >= 1) return(-Inf)
    (a - 1) * log(x) + (b - 1) * log(1 - x)
  }

  res2 <- run_slice_sampler(log_beta, n = 5000, x0 = 0.3,
                             w = 0.3, lower = 0, upper = 1,
                             a = 2, b = 5)


  ## --- 3. Generalised Pareto  GPD(sigma=1, xi=0.2) on (0, Inf) -------------
  log_gpd <- function(x, sigma, xi) {
    if (x <= 0) return(-Inf)
    if (xi != 0 && x > -sigma / xi) return(-Inf)
    -log(sigma) - (1 / xi + 1) * log(1 + xi * x / sigma)
  }

  res3 <- run_slice_sampler(log_gpd, n = 20000, x0 = 1,
                             w = 100, lower = 0,
                             sigma = 1, xi = 0.2)

  ## --- 4. Bimodal mixture of Gaussians (stress test step-out) --------------
  log_mix <- function(x) {
    log(0.4 * exp(-0.5 * (x + 3)^2) +
        0.6 * exp(-0.5 * (x - 3)^2))
  }

  res4 <- run_slice_sampler(log_mix, n = 10000, x0 = 0, w = 2)


  log_megpd_biv <- function(x, x_prev, kappa, sigma, xi, delta) {
    log(dmegpd_biv(x, x_prev, kappa, sigma, xi, delta))
  }

  res5 <- run_slice_sampler(log_megpd_biv,
    n = 100000,
    x0 = 1,
    w = 0.5,
    lower = 0,
    kappa = 2,
    sigma = 1,
    xi = 0.5,
    x_prev = 100,
    delta = delta_f_strong_upper
  )
}
