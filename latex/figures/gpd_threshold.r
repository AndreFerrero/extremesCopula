# =========================================================
# GPD THRESHOLD STABILITY PLOT
#
# Produces:
#   gpd_threshold_stability.png
# =========================================================

# install.packages("ismev")

library(ismev)

set.seed(123)

# ---------------------------------------------------------
# Simulate heavy-tailed data
# ---------------------------------------------------------
x <- evd::rfrechet(3000, shape = 4)

# ---------------------------------------------------------
# Threshold grid
# ---------------------------------------------------------
u_grid <- quantile(x, seq(0.70, 0.97, by = 0.01))

xi_hat <- numeric(length(u_grid))
se_hat <- numeric(length(u_grid))

# ---------------------------------------------------------
# Fit GPD for each threshold
# ---------------------------------------------------------
for(i in seq_along(u_grid)){

  fit <- gpd.fit(x,
                 threshold = u_grid[i],
                 show = FALSE)

  xi_hat[i] <- fit$mle[2]
  se_hat[i] <- sqrt(diag(fit$cov))[2]
}

# =========================================================
# Plot
# =========================================================

png("slides/figures/gpd_threshold_stability.png",
    width = 1200,
    height = 650,
    res = 150)

par(mar = c(5,5,3,1))

plot(u_grid,
     xi_hat,
     type = "b",
     pch = 16,
     lwd = 2,
     xlab = "Threshold  u",
     ylim = c(min(xi_hat - 1.96*se_hat), max(xi_hat + 1.96*se_hat)),
     ylab = expression(hat(xi)(u)),
     main = expression(X %~% Fréchet(alpha == 4)))

# Confidence bands
lines(u_grid,
      xi_hat + 1.96*se_hat,
      lty = 2)

lines(u_grid,
      xi_hat - 1.96*se_hat,
      lty = 2)

true_xi <- 1/4

abline(h = true_xi,
       lty = 2,
       lwd = 2,
       col = "red")

text(2.2,
     true_xi + 0.04,
     expression(xi == 0.25),
     col = "red")

dev.off()