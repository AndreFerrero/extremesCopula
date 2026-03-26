library(mev)
library(evd)
library(extRemes) # For standard GPD fitting
?fit.extgp

source("code/models/margins/egp.r")

x <- margin_egp$sample(100, c(mu = 0, kappa = 1.5, sigma = 2, xi = 0.2))
y <- rfrechet(100, scale = 1, shape = 4)
x_norm <- rnorm(1e6, mean = 4, sd = 1)

p <- 0.85
plot(x_norm)
abline(h = quantile(x_norm, p), col = "red", lty = 2) # Threshold for GPD fit
sum(x_norm > quantile(x_norm, p)) # Number of exceedances
hist(x_norm[x_norm > quantile(x_norm, p)], breaks = 20, main = "Exceedances over 85th Percentile", xlab = "Exceedance Value")

plot(y)
abline(h = quantile(y, p), col = "red", lty = 2) # Threshold for GPD fit
sum(y > quantile(y, p)) # Number of exceedances
hist(y[y > quantile(y, p)], breaks = 20, main = "Exceedances over 95th Percentile", xlab = "Exceedance Value")

# 1. Fit the EGPD to ALL data (Global Fit)
# 'model = 1' in mev corresponds to your [GPD]^kappa model
fit_egpd <- fit.extgp(x_norm,
    model = 1, init = list(kappa = 1, sigma = sd(x_norm), xi = 0.1),
    method = "mle"
)

params_egpd <- fit_egpd$fit$mle

# params_egpd contains: [kappa, sigma, xi]

# 2. Fit a standard GPD to the TAIL only (Benchmark)
# Pick a high threshold (e.g., 95th percentile)
u <- quantile(x_norm, p)
fit_gpd <- fevd(x_norm, threshold = u, type = "GP", method = "MLE")
params_gpd <- list(sigma = fit_gpd$results$par[1], xi = fit_gpd$results$par[2])

plot(fit_gpd)
summary(fit_gpd)

# 3. Compare the Shape Parameter (xi)
# This is the most critical parameter for return levels!
print(paste("EGPD Global Xi:", round(params_egpd["xi"], 4)))
print(paste("GPD Tail-only Xi:", round(params_gpd$xi, 4)))

# 4. Compare Return Levels (e.g., 100-year)
# Using the formula derived earlier for a specific r and ny
r <- 100
ny <- 365 # or 8760 for hourly
p_target <- 1 - 1 / (r * ny)

# EGPD Return Level
z_egpd <-(params_egpd["sigma"] / params_egpd["xi"]) *
    ((1 - p_target^(1 / params_egpd["kappa"]))^(-params_egpd["xi"]) - 1)

# GPD Return Level (standard POT formula)
# Note: zeta_u is the probability of being above the threshold
zeta_u <- mean(x_norm > u)
z_gpd <- u + (params_gpd$sigma / params_gpd$xi) *
    ((((r * ny) * zeta_u)^params_gpd$xi) - 1)

print(paste("EGPD 100-yr Level:", round(z_egpd, 2)))
print(paste("GPD 100-yr Level:", round(z_gpd, 2)))
