library(ggplot2)
library(patchwork) # For layout

# 1. THE BIVARIATE GUMBEL DENSITY FUNCTION (for visualization)
dgumbel <- function(u, v, theta) {
  # u, v, and theta can now be vectors of the same length
  # We ensure theta is slightly above 1 to avoid division by zero
  theta <- pmax(theta, 1.0001)
  
  lu <- -log(u)
  lv <- -log(v)
  term <- lu^theta + lv^theta
  r <- term^(1/theta)
  
  # Log-density calculation (more stable)
  log_dens <- -r + (theta - 1) * (log(lu) + log(lv)) + 
              (1/theta - 2) * log(term) + 
              log(r + theta - 1) - 
              (log(u) + log(v))
  
  return(exp(log_dens))
}

# 2. CREATE A GRID FOR VISUALIZATION
grid <- expand.grid(u = seq(0.01, 0.99, length.out = 50), 
                    v = seq(0.01, 0.99, length.out = 50))

# 3. PICK 3 THEAT VALUES FROM YOUR NEW PRIOR (Low, Med, High)
# Based on Gamma(2,2)+1
theta_low  <- 1.2   # Weak dependence (Tau ~ 0.16)
theta_med  <- 2.5   # Moderate dependence (Tau ~ 0.60)
theta_high <- 10.0  # Extreme dependence (Tau ~ 0.90)

plot_dens <- function(th, label) {
  grid$z <- dgumbel(grid$u, grid$v, th)
  
  ggplot(grid, aes(u, v, fill = z)) +
    geom_raster(interpolate = TRUE) +
    # FORCE the scale to be the same for all plots (e.g., 10^-2 to 10^2)
    scale_fill_viridis_c(option = "magma", 
                         trans = "log10", 
                         limits = c(0.01, 100), 
                         oob = scales::squish) + 
    labs(title = label, subtitle = paste0("theta = ", th)) +
    theme_minimal() + coord_fixed()
}

# 4. VIEW THE GEOMETRY
gridExtra::grid.arrange(plot_dens(theta_low, "Weak"), plot_dens(theta_med, "Moderate"), plot_dens(theta_high, "Extreme"), ncol = 3, nrow = 1)




# 5. BONUS: LIKELIHOOD SLICE
# Watch how the Likelihood "Sharpens" for a single observation (0.9, 0.9)
th_seq <- seq(1.001, 15, length.out = 200)
lik_slice <- dgumbel(rep(0.9, 200), rep(0.9, 200), th_seq)

plot(th_seq, lik_slice, type = "l", log="y", col="red", lwd=2,
     main="Likelihood Profile for Pair (0.9, 0.9)",
     xlab="Theta", ylab="Density (Log Scale)")
abline(v = 10, lty=2) # Showing where it becomes 'Extreme'