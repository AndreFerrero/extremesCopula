debug_megpd_step <- function(x_prev,
                             kappa,
                             xi,
                             delta_func,
                             n_grid = 1000) {

  par(mfrow = c(2,2))

  # -------------------------------------------------------
  # 1. GRID
  # -------------------------------------------------------
  x_grid <- make_log_grid(x_prev, n_grid = n_grid)

  plot(x_grid,
       rep(0, length(x_grid)),
       type = "l",
       main = "1. Log Grid (x-space)",
       xlab = "x",
       ylab = "",
       col = "grey")

  abline(v = x_prev, col = "red", lwd = 2)

  # -------------------------------------------------------
  # 2. CONDITIONAL DENSITY
  # -------------------------------------------------------
  pdf_vals <- conditional_pdf(x_prev, x_grid, kappa, xi, delta_func)

  plot(x_grid,
       pdf_vals,
       type = "l",
       main = "2. Conditional density f(x | x_prev)",
       xlab = "x",
       ylab = "density")

  abline(v = x_prev, col = "red", lwd = 2)

  # -------------------------------------------------------
  # 3. CDF via trapezoidal rule
  # -------------------------------------------------------
  cdf_obj <- pdf_to_cdf(x_grid, pdf_vals)

  plot(x_grid,
       cdf_obj$cdf,
       type = "l",
       main = "3. Numerical CDF F(x)",
       xlab = "x",
       ylab = "CDF")

  abline(h = seq(0,1,0.2), col = "grey90", lty = 3)

  # -------------------------------------------------------
  # 4. INVERSE CDF FUNCTION (visualization)
  # -------------------------------------------------------
  u_vals <- seq(0, 1, length.out = 1000)

  x_inv <- sapply(u_vals, function(u) sample_from_cdf(x_grid, cdf_obj$cdf, u))

  plot(u_vals,
       x_inv,
       type = "l",
       main = "4. Inverse CDF F^{-1}(u)",
       xlab = "u",
       ylab = "x")

  abline(v = c(0.25, 0.5, 0.75), col = "grey90", lty = 3)

  u_samp <- runif(10000)

  x_samp <- approx(cdf_obj$cdf, x_grid, u_samp)$y

  points(u_samp, x_samp, pch = 16, cex = 0.5, col = rgb(1,0,0,0.5))

  invisible(list(
    x_grid = x_grid,
    pdf = pdf_vals,
    cdf = cdf_obj$cdf
  ))
}

debug_megpd_step(
  x_prev = 1,
  kappa = 2,
  xi = 0.25,
  delta_func = function(r) {
    base <- 0.6
    valley <- dgamma(r, shape = 6, scale = 4)
    valley <- valley / max(valley)
    base - 0.2 * valley
  },
  n_grid = 1000
)
