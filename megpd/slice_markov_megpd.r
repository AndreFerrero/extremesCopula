###########################################################
# MEGPD Markov Chain Sampler — single slice step per transition
#
# State:       X_t  (positive real)
# Transition:  X_{t+1} via one slice sampler step targeting
#              f(x2 | X_t) on the log scale y = log(x2)
#
# Theoretical justification:
#   A single slice step is a valid transition kernel leaving
#   f(x2 | x1) invariant (Neal 2003). Embedding it inside a
#   Gibbs-style sweep preserves the joint stationary distribution
#   f(x1, x2) on the extended state space (x, u) where u is the
#   slice level — Tierney (1994). No inner burn-in needed.
#
# Log scale:
#   y = log(x2), so x2 = exp(y), dy/dx2 = 1/x2
#   log_target_log_scale(y | x1) =
#       log f(exp(y), x1) + y        [+y is the log-Jacobian]
#
# Adaptive w:
#   w is updated as a running MAD of log(X_t) over a window,
#   clipped to a reasonable range. This tracks the local spread
#   of the chain on the log scale, which is the right scale for
#   the step-out procedure.
###########################################################

library(egpd)
library(ggplot2)
library(patchwork)

###########################################################
# 1.  JOINT DENSITY
###########################################################

dmegpd_biv <- function(x1, x2, kappa, sigma, xi, delta_func) {
  x1 <- max(x1, 1e-10)
  x2 <- max(x2, 1e-10)
  r  <- x1 + x2

  f_r <- egpd::degpd_density(r, kappa = kappa, sigma = sigma, xi = xi)
  if (!is.finite(f_r) || f_r <= 0) return(0)

  d      <- max(delta_func(r), 0.01)
  lrat   <- log(x1 / x2)
  f_ang  <- dnorm(lrat, mean = 0, sd = d) * (r / (x1 * x2))

  val <- f_r * f_ang
  if (is.finite(val) && val > 0) val else 0
}

# Log target on the ORIGINAL scale (used only for reference / testing)
log_megpd_biv <- function(x2, x_prev, kappa, sigma, xi, delta_func) {
  v <- dmegpd_biv(x_prev, x2, kappa, sigma, xi, delta_func)
  if (v > 0) log(v) else -Inf
}

# Log target on the LOG scale  y = log(x2)
# = log f(x_prev, exp(y)) + y    [+y is the log-Jacobian for dy = dx2/x2]
log_target_log_scale <- function(y, x_prev, kappa, sigma, xi, delta_func) {
  x2 <- exp(y)
  lv <- log_megpd_biv(x2, x_prev, kappa, sigma, xi, delta_func)
  lv + y   # Jacobian term; -Inf + y = -Inf so no special casing needed
}

###########################################################
# 2.  SLICE SAMPLER  (Neal 2003, stepping-out + shrinkage)
#     Generic — works on any log_target
###########################################################

slice_sample <- function(log_target, n, x0,
                         w     = 1,
                         m     = Inf,
                         lower = -Inf,
                         upper =  Inf,
                         ...) {
  stopifnot(is.function(log_target))
  stopifnot(n > 0, w > 0)

  samples <- numeric(n)
  n_evals <- 0L
  x_cur   <- x0

  for (i in seq_len(n)) {
    # 1. Slice level
    lf_cur  <- log_target(x_cur, ...)
    n_evals <- n_evals + 1L
    log_y   <- lf_cur - rexp(1)

    # 2. Step-out
    u <- runif(1, 0, w)
    L <- max(x_cur - u,       lower)
    R <- min(x_cur + (w - u), upper)

    if (is.finite(m)) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1L) - J
      while (J > 0 && L > lower && log_target(L, ...) > log_y) {
        L <- max(L - w, lower); J <- J - 1L; n_evals <- n_evals + 1L
      }
      while (K > 0 && R < upper && log_target(R, ...) > log_y) {
        R <- min(R + w, upper); K <- K - 1L; n_evals <- n_evals + 1L
      }
    } else {
      while (L > lower && log_target(L, ...) > log_y) {
        L <- max(L - w, lower); n_evals <- n_evals + 1L
      }
      while (R < upper && log_target(R, ...) > log_y) {
        R <- min(R + w, upper); n_evals <- n_evals + 1L
      }
    }

    # 3. Shrinkage
    repeat {
      x_new  <- runif(1, L, R)
      lf_new <- log_target(x_new, ...)
      n_evals <- n_evals + 1L
      if (lf_new >= log_y) break
      if (x_new < x_cur) L <- x_new else R <- x_new
    }

    x_cur      <- x_new
    samples[i] <- x_cur
  }

  list(samples = samples, n_evals = n_evals)
}

###########################################################
# 3.  ADAPTIVE STEP-WIDTH TRACKER
#
#     Maintains a running window of the last `window` values
#     of log(X_t) and sets w = scale_factor * MAD, clipped
#     to [w_min, w_max].
#
#     Why MAD of log(X)?
#       - The slice sampler operates on y = log(x2), so the
#         natural step unit is in log-space.
#       - MAD is robust to occasional large excursions.
#       - We multiply by 1.4826 to make MAD consistent with
#         the Gaussian SD, then by scale_factor (default 2)
#         so the step-out interval covers roughly ±1 SD.
###########################################################

make_w_adapter <- function(window       = 100,
                           scale_factor = 2,
                           w_init       = 1,
                           w_min        = 0.05,
                           w_max        = 10) {
  buf <- rep(NA_real_, window)
  ptr <- 0L
  n   <- 0L
  w   <- w_init

  list(
    update = function(log_x) {
      ptr <<- (ptr %% window) + 1L
      buf[ptr] <<- log_x
      n <<- min(n + 1L, window)
      if (n >= 10L) {
        mad_val <- median(abs(buf[!is.na(buf)] -
                               median(buf[!is.na(buf)])))
        w <<- max(w_min, min(w_max, scale_factor * 1.4826 * mad_val))
      }
    },
    get_w = function() w
  )
}

###########################################################
# 4.  MARKOV CHAIN TRANSITION  (single slice step)
#
#     Given X_t = x_prev:
#       - current position on log scale: y_cur = log(x_prev)
#         (valid start because the chain is symmetric in x1, x2,
#          so f(x_prev | x_prev) > 0 generically)
#       - run slice_sample(n = 1) on log_target_log_scale
#       - return exp(y_new) as X_{t+1}
#
#     Note: we initialise the slice at y_cur = log(x_prev).
#     This is valid: the extended-space argument requires only
#     that the initial point is in the support of the target,
#     not that it is a draw from it.
###########################################################

megpd_slice_step <- function(x_prev, kappa, sigma, xi, delta_func, w) {
  y_cur <- log(x_prev)

  res <- slice_sample(
    log_target = log_target_log_scale,
    n          = 1L,
    x0         = y_cur,
    w          = w,
    lower      = -Inf,   # log scale: x2 in (0,inf) => y in (-inf, inf)
    upper      =  Inf,
    # named args forwarded via ...
    x_prev     = x_prev,
    kappa      = kappa,
    sigma      = sigma,
    xi         = xi,
    delta_func = delta_func
  )

  exp(res$samples[1L])
}

###########################################################
# 5.  FULL CHAIN SIMULATOR
###########################################################

simulate_megpd_slice <- function(n_steps       = 5000,
                                 kappa, sigma, xi,
                                 delta_func,
                                 x0            = 1,
                                 burn_in_prop  = 0.10,
                                 w_init        = 1,
                                 w_window      = 100,
                                 w_scale       = 2,
                                 show_progress = TRUE) {
  stopifnot(burn_in_prop >= 0, burn_in_prop < 1)

  x       <- numeric(n_steps)
  x[1]    <- x0
  adapter <- make_w_adapter(window       = w_window,
                            scale_factor = w_scale,
                            w_init       = w_init)
  adapter$update(log(x0))

  w_trace <- numeric(n_steps)   # store w at each step for diagnostics
  w_trace[1] <- w_init

  if (show_progress) pb <- txtProgressBar(min = 0, max = n_steps, style = 3)

  for (t in 2:n_steps) {
    w_cur  <- adapter$get_w()
    x[t]   <- megpd_slice_step(x[t - 1], kappa, sigma, xi, delta_func,
                                w = w_cur)
    adapter$update(log(x[t]))
    w_trace[t] <- w_cur

    if (show_progress) setTxtProgressBar(pb, t)
  }

  if (show_progress) close(pb)

  burn_in <- floor(n_steps * burn_in_prop)
  list(
    full_chain  = x,
    final_chain = x[(burn_in + 1):n_steps],
    w_trace     = w_trace,
    burn_in     = burn_in
  )
}

###########################################################
# 6.  DELTA FUNCTIONS
###########################################################

# Decaying: strong upper-tail dependence
delta_strong_upper <- function(r) 0.4 + 0.4 * exp(-r / 5)
curve(delta_strong_upper, to = 20)

###########################################################
# 7.  RUN
###########################################################

set.seed(42)

kappa_val <- 2
sigma_val <- 1
xi_val    <- 0.1
n_steps   <- 10000

cat("Running MEGPD slice-sampler chain (n =", n_steps, ")...\n")

sim <- simulate_megpd_slice(
  n_steps      = n_steps,
  kappa        = kappa_val,
  sigma        = sigma_val,
  xi           = xi_val,
  delta_func   = delta_nonmonotone,
  x0           = 1,
  burn_in_prop = 0.30,
  w_init       = 1,    # initial step width on log scale
  w_window     = 100,  # MAD window
  w_scale      = 2     # w = 2 * 1.4826 * MAD
)

final_chain <- sim$final_chain
cat("Chain length after burn-in:", length(final_chain), "\n")
cat("Summary:\n"); print(summary(final_chain))


plot.ts(sim$w_trace)
