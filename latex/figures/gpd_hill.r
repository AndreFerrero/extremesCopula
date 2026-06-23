# =========================================================
# EXTREME VALUE DIAGNOSTICS
# Produces:
#   gpd_threshold_stability.png
#   hill_plot.png
# =========================================================

library(ismev)
library(evir)
library(ReIns)
library(evd)

# =========================================================
# GLOBAL PARAMETERS
# =========================================================

seed_value      <- 123

# Data generation
n_sample        <- 3000
frechet_alpha   <- 4
true_xi         <- 1 / frechet_alpha

# Threshold stability settings
threshold_probs <- seq(0.80, 0.97, by = 0.01)

# Hill plot settings
hill_k_min      <- 1
hill_k_max      <- 2000

# Output files
gpd_file        <- "slides/figures/gpd_threshold_stability.png"
hill_file       <- "slides/figures/hill_plot.png"

# ---------------------------------------------------------
# Plot-specific parameters (MANUAL ADJUSTMENTS HERE)
# ---------------------------------------------------------

# GPD plot
gpd_width       <- 1200
gpd_height      <- 650
gpd_res         <- 150
gpd_mar         <- c(5,5,3,1)

# Updated x-position to reflect the quantile scale
gpd_label_x     <- 0.95 
gpd_label_y_off <- 0.04

# Hill plot
hill_width      <- 1400
hill_height     <- 700
hill_res        <- 180
hill_mar        <- c(4,4,2,1)

hill_label_x     <- 300
hill_label_y_off <- 0.03

# =========================================================
# COMMON DATA (USED BY BOTH PLOTS)
# =========================================================

set.seed(seed_value)

x <- evd::rfrechet(
  n_sample,
  shape = frechet_alpha
)

# =========================================================
# GPD THRESHOLD STABILITY
# =========================================================

u_grid <- quantile(
  x,
  threshold_probs
)

xi_hat <- numeric(length(u_grid))
se_hat <- numeric(length(u_grid))

for(i in seq_along(u_grid)) {

  fit <- ismev::gpd.fit(
    x,
    threshold = u_grid[i],
    show = FALSE
  )

  xi_hat[i] <- fit$mle[2]
  se_hat[i] <- sqrt(diag(fit$cov))[2]
}

png(
  gpd_file,
  width = gpd_width,
  height = gpd_height,
  res = gpd_res
)

par(mar = gpd_mar)

plot(
  threshold_probs,
  xi_hat,
  type = "b",
  pch = 16,
  lwd = 2,
  xlab = "Threshold Quantile Level (p)",
  ylim = c(
    min(xi_hat - 1.96 * se_hat),
    max(xi_hat + 1.96 * se_hat)
  ),
  ylab = expression(hat(xi)(p)),
  main = expression(X %~% Fréchet(alpha == 4))
)

lines(
  threshold_probs,
  xi_hat + 1.96 * se_hat,
  lty = 2
)

lines(
  threshold_probs,
  xi_hat - 1.96 * se_hat,
  lty = 2
)

abline(
  h = true_xi,
  lty = 2,
  lwd = 2,
  col = "red"
)

text(
  gpd_label_x,
  true_xi + gpd_label_y_off,
  expression(xi == 0.25),
  col = "red"
)

dev.off()

# =========================================================
# HILL PLOT (Kept same as original)
# =========================================================

h <- ReIns::Hill(x)
k_full  <- h$k
xi_full <- h$gamma
idx <- k_full <= hill_k_max & k_full > hill_k_min
k  <- k_full[idx]
xi <- xi_full[idx]

png(
  hill_file,
  width = hill_width,
  height = hill_height,
  res = hill_res
)

par(mar = hill_mar, bty = "l")

plot(
  k,
  xi,
  type = "l",
  lwd = 2,
  col = "black",
  xlim = c(hill_k_min, hill_k_max),
  ylim = c(min(c(xi, true_xi - 0.05)), max(xi)),
  xlab = "k",
  ylab = expression(hat(xi)[H](k)),
  main = expression(X %~% Fréchet(alpha == 4)),
  las = 1
)

abline(h = true_xi, lty = 2, lwd = 2, col = "red")

text(
  hill_label_x,
  true_xi + hill_label_y_off,
  expression(xi == 0.25),
  col = "red"
)

dev.off()