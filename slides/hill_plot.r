# =========================================================
# Minimal Hill plot
# Focus only on first 300 k values
# =========================================================

library(evir)

set.seed(123)

# ---------------------------------------------------------
# Simulate heavy-tailed sample
# ---------------------------------------------------------
x <- abs(rt(4000, df = 2))

# ---------------------------------------------------------
# Hill estimates
# ---------------------------------------------------------
h <- ReIns::Hill(x)

k_full  <- h$k
xi_full <- h$gamma

# ---------------------------------------------------------
# Keep only first 300 k values
# ---------------------------------------------------------
idx <- k_full <= 1000 & k_full > 15

k  <- k_full[idx]
xi <- xi_full[idx]

# ---------------------------------------------------------
# Stable region
# ---------------------------------------------------------
k1 <- 100
k2 <- 180

stable_level <- mean(
  xi[k >= k1 & k <= k2],
  na.rm = TRUE
)

# =========================================================
# Plot
# =========================================================

png("slides/figures/hill_plot.png",
    width = 1200,
    height = 700,
    res = 180)

par(mar = c(4,4,2,1),
    bty = "l")

plot(k, xi,
     type = "l",
     lwd = 2,
     col = "black",
     xlim = c(15, 1000),
     ylim = range(xi),
     xlab = "",
     ylab = expression(hat(xi)[H](k)),
     main = "",
     las = 1,
     xaxt = "n")
axis(1, at = seq(15, 1000, by = 50))

# ---------------------------------------------------------
# Stable region shading
# ---------------------------------------------------------

# Replot curve
lines(k, xi,
      lwd = 2)

true_xi <- 1/2

abline(h = true_xi,
       lty = 2,
       lwd = 2,
       col = "red")

text(250,
     true_xi + 0.01,
     expression(xi == 0.5),
     col = "red")

dev.off()