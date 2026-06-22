par(mfrow = c(2, 2))

x_prev <- 100
kappa <- 2
sigma <- 1
xi <- 0.5
delta_func <- delta_strong_upper
n_grid <- 500

# -------------------------------------------------------
# 1. U-GRID
# -------------------------------------------------------

u_grid <- make_adaptive_grid(x_to_u(x_prev), n_grid, phi = 50)

hist(u_grid, main = "U-grid")

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
  ylab = "density"
)
rug(u_grid)

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
  ylab = "F(u)"
)

abline(h = seq(0, 1, 0.2), col = "grey90", lty = 3)

cat("\n--- CDF diagnostics ---\n")
cat("Total mass (pre-normalization):", cdf_obj$total_mass, "\n")
cat("CDF start:", cdf_obj$cdf[1], "\n")
cat("CDF end:", tail(cdf_obj$cdf, 1), "\n")
cat("Monotone:", all(diff(cdf_obj$cdf) >= -1e-10), "\n")
cat("In [0,1]:", all(cdf_obj$cdf >= -1e-10 & cdf_obj$cdf <= 1 + 1e-10), "\n")

# -------------------------------------------------------
# 4. INVERSE CDF CHECK
# -------------------------------------------------------
u_test <- seq(0, 1, length.out = 200)

x_inv <- sapply(u_test, function(u) {
  sample_from_cdf(u_grid, cdf_obj$cdf, u)
})

plot(u_test,
  x_inv,
  type = "l",
  main = "4. Inverse CDF map u → u_sample",
  xlab = "u",
  ylab = "u_sample"
)

abline(a = 0, b = 1, col = "grey60", lty = 2)

# -------------------------------------------------------
# 5. SAMPLING CHECK (distribution sanity)
# -------------------------------------------------------
par(mfrow = c(1, 2))

u_samp <- runif(5000)

u_draws <- sapply(u_samp, function(u) {
  sample_from_cdf(u_grid, cdf_obj$cdf, u)
})

hist(u_draws,
  breaks = 50,
  main = "5. Samples in u-space",
  xlab = "u (should match density shape)",
  col = "grey"
)

abline(v = 0.5, col = "red", lwd = 2)

# -------------------------------------------------------
# 6. BACK TRANSFORM CHECK
# -------------------------------------------------------
x_draws <- u_to_x(u_draws)

hist(x_draws,
  breaks = 50,
  main = "6. Samples in x-space",
  xlab = "x",
  col = "grey"
)

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
