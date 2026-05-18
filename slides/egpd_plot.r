# =========================================================
# Extended GPD illustration via power transformation
# =========================================================

library(evir)

set.seed(123)

# ---------------------------------------------------------
# 1. Simulate baseline GPD sample
# ---------------------------------------------------------
n <- 100000
xi <- 0.5
beta <- 1

u <- evir::rgpd(n = n, xi = xi)

# empirical CDF function
Fn <- ecdf(u)

# ---------------------------------------------------------
# 2. Define EGPD transformation
#    F_kappa(x) = [F(x)]^kappa
#    density changes accordingly (illustrative)
# ---------------------------------------------------------
kappa_vals <- c(2, 4)

# grid for density illustration
x_grid <- seq(0, quantile(u, 0.99), length.out = 400)

# baseline GPD density
gpd_dens <- evir::dgpd(x_grid, xi = xi)

# ---------------------------------------------------------
# approximate EGPD density via transformation rule:
# f_kappa(x) ≈ kappa * F(x)^(kappa-1) * f(x)
# ---------------------------------------------------------
egpd_density <- function(kappa) {
       Fx <- Fn(x_grid)
       kappa * (Fx^(kappa - 1)) * gpd_dens
}

dens_list <- lapply(kappa_vals, egpd_density)

# =========================================================
# PLOT 1: Density comparison
# =========================================================

png("slides/figures/egpd_plot.png",
    width = 1200,
    height = 600,
    res = 180)

par(
  mar = c(3, 4, 2, 1),
  mgp = c(2, 0.6, 0),
  las = 1
)

plot(x_grid, gpd_dens,
     type = "l",
     lwd = 2.5,
     col = "black",
     xlab = expression(x),
     ylab = expression(f(x)),
     bty = "l",
     xaxt = "n")

axis(1, tck = -0.02)
grid(col = "grey90")

lines(x_grid, dens_list[[1]],
      col = adjustcolor("steelblue4", 0.9),
      lwd = 2)

lines(x_grid, dens_list[[2]],
      col = adjustcolor("firebrick3", 0.9),
      lwd = 2)

legend("topright",
       legend = c("GPD (κ = 1)",
                  "EGPD (κ = 2)",
                  "EGPD (κ = 4)"),
       col = c("black", "steelblue4", "firebrick3"),
       lwd = 2,
       bty = "n")

dev.off()

# =========================================================
# TAIL EQUIVALENCE VISUALISATION
# log-survival comparison
# =========================================================

# grid
x_grid <- seq(quantile(u, 0.5), quantile(u, 0.999), length.out = 10000)

# empirical CDF
Fn <- ecdf(u)

# GPD survival
S_gpd <- 1 - Fn(x_grid)

# EGPD survival (power transform)
egpd_surv <- function(kappa) {
  1 - Fn(x_grid)^kappa
}

S_k2 <- egpd_surv(2)
S_k4 <- egpd_surv(4)

# avoid log(0)
eps <- 1e-10

png("slides/figures/egpd_tail_equivalence.png",
    width = 1200, height = 600, res = 180)

par(
  mar = c(2.8,4,1,1),
  mgp = c(2,0.7,0)
)

plot(log(x_grid), log(pmax(S_gpd, eps)),
     type = "l",
     lwd = 2,
     col = "black",
     xlab = expression(log(x)),
     ylab = expression(log(bar(F)(x))),
     main = "",
     bty = "l",
     xaxt = "n")
axis(1, tck = -0.02)
grid(col = "grey90")

lines(log(x_grid), log(pmax(S_k2, eps)),
      col = "steelblue4",
      lwd = 2)

lines(log(x_grid), log(pmax(S_k4, eps)),
      col = "firebrick3",
      lwd = 2)

abline(a = 0.5,                 # intercept (adjust visually)
       b = -1/xi,             # slope
       col = "gray50",
       lty = 2,
       lwd = 2)

legend("topright",
       legend = c(expression(GPD~~kappa==1),
                  expression(EGPD~~kappa==2),
                  expression(EGPD~~kappa==4),
                  expression("slope ="~~-1/xi)),
       col = c("black", "steelblue4", "firebrick3", "gray50"),
       lwd = 2,
       lty = c(1,1,1,2),
       bty = "n")

dev.off()