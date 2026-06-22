# ---- Your density function ----
dmegpd_biv <- function(x1, x2, kappa, xi, delta_func) {
  x1 <- max(x1, 1e-10)
  x2 <- max(x2, 1e-10)
  
  r <- x1 + x2
  
  H_r <- if(abs(xi) < 1e-10) {
    1 - exp(-r)
  } else {
    pmax(0, 1 - (1 + xi * r)^(-1/xi))
  }
  
  h_r <- if(abs(xi) < 1e-10) {
    exp(-r)
  } else {
    pmax(0, (1 + xi * r)^(-1/xi - 1))
  }
  
  d <- max(delta_func(r), 0.01)
  
  term_rad <- kappa * h_r * (H_r^(kappa - 1))
  
  log_ratio <- log(x1 / x2)
  z_sq <- (log_ratio / d)^2
  
  term_ang <- (1 / (sqrt(2 * pi) * d)) * exp(-0.5 * z_sq) * (r / (x1 * x2))
  
  val <- term_rad * term_ang
  ifelse(is.finite(val), val, 0)
}

# ---- Example delta function (replace with yours) ----
delta_func <- function(r) {
  0.3 + 0.1 * log1p(r)
}

# ---- Parameters ----
kappa <- 2
xi <- 0.2

# ---- Grid setup ----
x_seq <- seq(0.01, 0.4, length.out = 800)
grid <- expand.grid(x1 = x_seq, x2 = x_seq)

# ---- Evaluate density on grid ----
grid$z <- mapply(
  dmegpd_biv,
  grid$x1,
  grid$x2,
  MoreArgs = list(kappa = kappa, xi = xi, delta_func = delta_func)
)

# reshape into matrix for plotting
z_matrix <- matrix(grid$z, nrow = length(x_seq), byrow = TRUE)

# ---- 1. Heatmap ----
filled.contour(
  x_seq, x_seq, z_matrix,
  color.palette = terrain.colors,
  xlab = "x1",
  ylab = "x2",
  main = "Joint Density Heatmap (dmegpd_biv)"
)

# ---- 2. Contour plot ----
contour(
  x_seq, x_seq, z_matrix,
  nlevels = 20,
  xlab = "x1",
  ylab = "x2",
  main = "Contour Plot of Joint Density"
)

# ---- 3. 3D surface (perspective plot) ----
persp(
  x_seq, x_seq, z_matrix,
  theta = 40,
  phi = 25,
  expand = 1,
  col = "lightblue",
  ticktype = "detailed",
  xlab = "x1",
  ylab = "x2",
  zlab = "density"
)

# ---- 4. Quick multimodality diagnostic ----
max_idx <- which(z_matrix == max(z_matrix), arr.ind = TRUE)
cat("Peak density at:\n")
cat("x1 =", x_seq[max_idx[1]], "\n")
cat("x2 =", x_seq[max_idx[2]], "\n")


# ---- Grid ----
trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

x_seq <- seq(0.001, 8, length.out = 200)

# ---- Marginal of X1: integrate over x2 ----
marg_x1 <- numeric(length(x_seq))

for (i in seq_along(x_seq)) {
  x1_fixed <- x_seq[i]
  
  vals <- sapply(x_seq, function(x2) {
    dmegpd_biv(x1_fixed, x2, kappa, xi, delta_func)
  })
  
  marg_x1[i] <- trapz(x_seq, vals)
}

# ---- Marginal of X2: integrate over x1 ----
marg_x2 <- numeric(length(x_seq))

for (j in seq_along(x_seq)) {
  x2_fixed <- x_seq[j]
  
  vals <- sapply(x_seq, function(x1) {
    dmegpd_biv(x1, x2_fixed, kappa, xi, delta_func)
  })
  
  marg_x2[j] <- trapz(x_seq, vals)
}

# ---- simple trapezoidal rule ----


# ---- Plot marginals ----
par(mfrow = c(1, 2))

plot(x_seq, marg_x1, type = "l", lwd = 2,
     main = "Marginal of X1",
     xlab = "x1", ylab = "density")

plot(x_seq, marg_x2, type = "l", lwd = 2,
     main = "Marginal of X2",
     xlab = "x2", ylab = "density")
