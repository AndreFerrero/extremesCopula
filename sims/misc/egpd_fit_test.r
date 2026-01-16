# ------------------------------
# 1. Load required packages
# ------------------------------
source("libs/packages.R")
library(mev)

# ------------------------------
# 2. Source your libs
# ------------------------------
source("libs/models/margins/frechet.r")
source("libs/models/margins/lognormal.r")
source("libs/models/margins/egp.r")

source("libs/models/copulas/gumbel.R")
source("libs/models/builders/simulator.R")

set.seed(123)

n <- 1000
param_frechet <- c(scale = 1, shape = 4)  # Fréchet(alpha = 2)
param_lognorm <- c(mu = 0, sigma = 1)

U <- runif(n)
X_frechet <- margin_frechet$quantile(U, param_frechet)
X_lognorm <- margin_lognormal$quantile(U, param_lognorm)
X_gp <- evd::rgpd(n, shape = 0.5)

# mod_pwm <- fit.extgp(X, init = c(1, 1, 0.1), method = "pwm")

# mod_pwm$fit

mod_mle_frechet <- fit.extgp(X_frechet, init = c(1, 1, 0.1), method = "mle", confint = TRUE, R = 20)
round(mod_mle_frechet$fit$mle, 2)

mod_mle_lognorm <- fit.extgp(X_lognorm, init = c(1, 1, 0.1), method = "mle", confint = TRUE, R = 20)
round(mod_mle_lognorm$fit$mle, 2)

mod_mle_gp <- fit.extgp(X_gp, init = c(2, 1, 0.1), method = "mle", confint = TRUE, R = 20)
round(mod_mle_gp$fit$mle, 2)


# Copula part
param_map <- list(
  margin = c("mu", "sigma"),
  copula = "theta"
)

lognorm_gumbel_sim <- build_simulator(
    copula = copula_gumbel,
    margin = margin_lognormal,
    param_map = param_map
)

param_lognorm_gumb <- c(param_lognorm, theta = 2)

X_lognorm_gumb <- lognorm_gumbel_sim(param_lognorm_gumb, n)

mod_mle_lognorm_gumb <- fit.extgp(X_lognorm_gumb, init = c(1, 1, 0.1), method = "mle", confint = TRUE, R = 20)
round(mod_mle_lognorm_gumb$fit$mle, 2)
