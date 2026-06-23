# =========================================================
# Gumbel copula pair plots in BASE R
# theta = 1, 2, 4, 6
# =========================================================

library(copula)

set.seed(123)

theta_vals <- c(1, 2, 4, 6)

# ---------------------------------------------------------
# Simulation function
# ---------------------------------------------------------
simulate_gumbel <- function(theta, n = 5000){

  cop <- gumbelCopula(param = theta, dim = 2)
  u <- rCopula(n, cop)

  list(U1 = u[,1], U2 = u[,2])
}

# =========================================================
# 1. FULL SAMPLE PLOTS
# =========================================================

png("slides/figures/gumbel_copula.png",
    width = 1800,
    height = 900,
    res = 150)

par(mfrow = c(2, 2),
    mar = c(2.5, 2.5, 2, 1))

for(theta in theta_vals){

  sim <- simulate_gumbel(theta)

  plot(sim$U1, sim$U2,
       pch = 16,
       cex = 0.4,
       col = rgb(0, 0, 0, 0.35),
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = expression(U[1]),
       ylab = expression(U[2]),
       main = bquote(theta == .(theta)))
}

par(mfrow = c(1, 1))
dev.off()



# =========================================================
# 2. UPPER TAIL DEPENDENCE (TOP 10%)
# =========================================================

png("slides/figures/gumbel_upper_tail.png",
    width = 1800,
    height = 900,
    res = 150)

par(mfrow = c(2, 2),
    mar = c(2.5, 2.5, 2, 1))

for(theta in theta_vals){

  sim <- simulate_gumbel(theta)

  plot(sim$U1, sim$U2,
       pch = 16,
       cex = 0.5,
       col = rgb(0, 0, 0, 0.45),

       # zoom into upper tail
       xlim = c(0.9, 1),
       ylim = c(0.9, 1),

       xlab = expression(U[1]),
       ylab = expression(U[2]),
       main = bquote(theta == .(theta)))

  # optional guide lines
  abline(v = 0.9, lty = 2, col = "red")
  abline(h = 0.9, lty = 2, col = "red")
}

par(mfrow = c(1, 1))
dev.off()



# =========================================================
# 3. LOWER TAIL DEPENDENCE (BOTTOM 10%)
# =========================================================

png("slides/figures/gumbel_lower_tail.png",
    width = 1800,
    height = 900,
    res = 150)

par(mfrow = c(2, 2),
    mar = c(2.5, 2.5, 2, 1))

for(theta in theta_vals){

  sim <- simulate_gumbel(theta)

  plot(sim$U1, sim$U2,
       pch = 16,
       cex = 0.5,
       col = rgb(0, 0, 0, 0.45),

       # zoom into lower tail
       xlim = c(0, 0.1),
       ylim = c(0, 0.1),

       xlab = expression(U[1]),
       ylab = expression(U[2]),
       main = bquote(theta == .(theta)))

  # optional guide lines
  abline(v = 0.1, lty = 2, col = "blue")
  abline(h = 0.1, lty = 2, col = "blue")
}

par(mfrow = c(1, 1))
dev.off()