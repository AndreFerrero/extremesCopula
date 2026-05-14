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
simulate_gumbel <- function(theta, n = 1000){

  cop <- gumbelCopula(param = theta, dim = 2)
  u <- rCopula(n, cop)

  list(U1 = u[,1], U2 = u[,2])
}



# ---------------------------------------------------------
# Loop over theta values
# ---------------------------------------------------------
png("slides/figures/gumbel_copula.png",
    width = 1800,
    height = 900,
    res = 150)
# ---------------------------------------------------------
# Plot layout (2x2 panels)
# ---------------------------------------------------------
par(mfrow = c(2, 2),
    mar = c(2.5, 2.5, 2, 1))  # bottom, left, top, right

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

# reset layout
par(mfrow = c(1, 1))

dev.off()
