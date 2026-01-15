# ------------------------------
# 1. Load required packages
# ------------------------------
source("libs/packages.R")

# ------------------------------
# 2. Source your libs
# ------------------------------
source("libs/models/copulas/gumbel.R")     # defines copula_gumbel
source("libs/models/copulas/clayton.r")    # defines copula_clayton
source("libs/models/margins/egp.r")  # defines margin_egp
source("libs/models/builders/simulator.R") # defines build_simulator())


param <- c(sigma = 1, xi = 0, kappa = 2)

xg <- seq(0, 5, length.out = 1000)
plot(xg, exp(margin_egp$lpdf(xg, param_test)), col = "red", type = "l", lwd = 2)
# lines(xg, evd::dgpd(xg, scale = param_test["sigma"], shape = param_test["xi"]))
abline(0, 0)

X <- margin_egp$sample(100000, param_test)

hist(X, freq = FALSE)
xg <- seq(0, quantile(X, 0.999), length.out = 1000)
lines(xg, exp(margin_egp$lpdf(xg, param_test)), col = "red", lwd = 2)

hist(margin_egp$cdf(X, param_test))

grid <- seq(0, 100, 0.1)

cdf <- margin_egp$cdf(grid, param_test)

plot(grid, cdf, type = "l")

quant <- margin_egp$quantile(seq(0.001, 0.999, 0.01), param_test)

# ------------------------------
# 3. Define parameter mapping
# ------------------------------
# param vector will be named
param_map <- list(
  margin = c("sigma", "xi", "kappa"),
  copula = "theta"
)

# ------------------------------
# 4. Build the simulator
# ------------------------------
simulator <- build_simulator(
  copula = copula_gumbel,
  margin = margin_egp,
  param_map = param_map
)

# ------------------------------
# 5. Test simulation
# ------------------------------
param_test <- c(sigma = 1, xi = 0.5, kappa = 2, theta = 2)
n_test <- 1000

X_sim <- simulator(param_test, n_test)

# ------------------------------
# Histogram of simulated data
# ------------------------------
hist(X_sim, breaks = 30,
     main = "Gumbel + EGPD",
     xlab = "X",
     col = "skyblue",
     probability = TRUE)
