debug_megpd_step_u <- function(x_prev,
                               kappa,
                               sigma,
                               xi,
                               delta_func,
                               n_grid = 1000) {

  par(mfrow = c(2,2))

  # -------------------------------------------------------
  # 1. U-GRID
  # -------------------------------------------------------

  make_hybrid_grid <- function(u_prev, n_grid = 1000) {
  # 500 points spread uniformly across (0,1)
  g_uni <- seq(0, 1, length.out = n_grid / 2)
  
  # 500 points concentrated around u_prev
  g_beta <- make_beta_grid(u_prev, n_grid = n_grid / 2, phi = 100)
  
  # Combine, sort, and remove duplicates
  sort(unique(c(g_uni, g_beta)))
}

  u_grid <- make_hybrid_grid(n_grid)

  plot(u_grid,
       rep(0, length(u_grid)),
       type = "l",
       main = "1. Grid in u-space (0,1)",
       xlab = "u",
       ylab = "",
       col = "grey")

  abline(v = c(0,1), col = "black", lwd = 2)
  abline(v = 0.5, col = "red", lwd = 2)

  # -------------------------------------------------------
  # 2. TRANSFORMED CONDITIONAL DENSITY
  # -------------------------------------------------------
  pdf_vals <- conditional_pdf_u(
    x_prev, u_grid, kappa, sigma, xi, delta_func
  )

  plot(u_grid,
       pdf_vals,
       type = "l",
       main = "2. Transformed density f_U(u | x_prev)",
       xlab = "u",
       ylab = "density")

  abline(v = 0.5, col = "red", lwd = 2)

  cat("\n--- PDF diagnostics ---\n")
  cat("Max pdf_U:", max(pdf_vals), "\n")
  cat("Min pdf_U:", min(pdf_vals), "\n")
  cat("Sum(pdf * du approx):",
      sum((pdf_vals[-1] + pdf_vals[-length(pdf_vals)]) / 2 *
            diff(u_grid)),
      "\n")

  # -------------------------------------------------------
  # 3. CDF via trapezoidal rule
  # -------------------------------------------------------
  cdf_obj <- pdf_to_cdf(u_grid, pdf_vals)

  if (is.null(cdf_obj)) {
    stop("CDF construction failed.")
  }

  plot(u_grid,
       cdf_obj$cdf,
       type = "l",
       main = "3. CDF in u-space",
       xlab = "u",
       ylab = "F(u)")

  abline(h = seq(0,1,0.2), col = "grey90", lty = 3)

  cat("\n--- CDF diagnostics ---\n")
  cat("Total mass (pre-normalization):", cdf_obj$total_mass, "\n")
  cat("CDF start:", cdf_obj$cdf[1], "\n")
  cat("CDF end:", tail(cdf_obj$cdf, 1), "\n")
  cat("Monotone:", all(diff(cdf_obj$cdf) >= -1e-10), "\n")
  cat("In [0,1]:", all(cdf_obj$cdf >= -1e-10 & cdf_obj$cdf <= 1+1e-10), "\n")

  # -------------------------------------------------------
  # 4. INVERSE CDF CHECK
  # -------------------------------------------------------
  u_test <- seq(0, 1, length.out = 200)

  x_inv <- sapply(u_test, function(u)
    sample_from_cdf(u_grid, cdf_obj$cdf, u)
  )

  plot(u_test,
       x_inv,
       type = "l",
       main = "4. Inverse CDF map u → u_sample",
       xlab = "u",
       ylab = "u_sample")

  abline(a = 0, b = 1, col = "grey60", lty = 2)

  # -------------------------------------------------------
  # 5. SAMPLING CHECK (distribution sanity)
  # -------------------------------------------------------
  par(mfrow = c(1, 2))

  u_samp <- runif(5000)

  u_draws <- sapply(u_samp, function(u)
    sample_from_cdf(u_grid, cdf_obj$cdf, u)
  )

  hist(u_draws,
       breaks = 50,
       main = "5. Samples in u-space",
       xlab = "u (should match density shape)",
       col = "grey")

  abline(v = 0.5, col = "red", lwd = 2)

  # -------------------------------------------------------
  # 6. BACK TRANSFORM CHECK
  # -------------------------------------------------------
  x_draws <- u_to_x(u_draws)

  hist(x_draws,
       breaks = 50,
       main = "6. Samples in x-space",
       xlab = "x",
       col = "grey")

  cat("\n--- Transform diagnostics ---\n")
  cat("Median u:", median(u_draws), "\n")
  cat("Median x:", median(x_draws), "\n")

  invisible(list(
    u_grid = u_grid,
    pdf_u = pdf_vals,
    cdf = cdf_obj$cdf,
    u_samples = u_draws,
    x_samples = x_draws
  ))
}

debug_megpd_step_u(10, 2, 1, 0.5, delta_strong_upper)
