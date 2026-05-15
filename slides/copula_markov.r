# =========================================================
# Copula Markov Chain Diagnostics
# ---------------------------------------------------------
# Top row    : latent uniforms U_t
# Bottom row : transformed series X_t
#
# Left  = histogram
# Right = lag-1 plot
# =========================================================

source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")

set.seed(123)

# ---------------------------------------------------------
# Model
# ---------------------------------------------------------

egpd_gumbel_model <- make_copula_markov_model(
  margin_egpd,
  copula_gumbel,
  stan_mod = NULL
)

margin_param <- c(
  mu    = 0,
  kappa = 2,
  sigma = 1,
  xi    = 0.1
)

copula_param <- 3
n <- 10000

# ---------------------------------------------------------
# Simulate
# ---------------------------------------------------------

egpd_gumbel_data <- egpd_gumbel_model$simulate(
  n = n,
  margin_param = margin_param,
  copula_param = copula_param
)

u <- egpd_gumbel_data$u
x <- egpd_gumbel_data$x

# lagged pairs
u_lag  <- u[-length(u)]
u_curr <- u[-1]

x_lag  <- x[-length(x)]
x_curr <- x[-1]

# =========================================================
# Plot
# =========================================================

png(
  "slides/figures/copula_markov.png",
  width = 1800,
  height = 1000,
  res = 150
)

# ---------------------------------------------------------
# Layout:
#
#  1 2
#  3 4
#
# 1 = histogram(U)
# 2 = lag plot(U_t-1, U_t)
# 3 = histogram(X)
# 4 = lag plot(X_t-1, X_t)
# ---------------------------------------------------------

par(
  mfrow = c(2, 2),
  mar = c(4, 4, 3, 1),
  mgp = c(2.2, 0.8, 0)
)

# ---------------------------------------------------------
# Histogram of U
# ---------------------------------------------------------

hist(
  u,
  breaks = 40,
  probability = TRUE,
  col = "grey80",
  border = "white",
  main = expression("Latent uniforms " ~ U[t]),
  xlab = expression(U[t])
)

abline(h = 1, lwd = 2, lty = 2, col = "red")

# ---------------------------------------------------------
# Lag plot of U
# ---------------------------------------------------------

plot(
  u_lag,
  u_curr,
  pch = 16,
  cex = 0.4,
  col = rgb(0, 0, 0, 0.35),

  xlab = expression(U[t-1]),
  ylab = expression(U[t]),
  main = expression("Copula dependence of " ~ U[t]),

  xlim = c(0, 1),
  ylim = c(0, 1)
)

# ---------------------------------------------------------
# Histogram of X
# ---------------------------------------------------------

hist(
  x,
  breaks = 40,
  probability = TRUE,
  col = "grey80",
  border = "white",
  main = expression("Observed process " ~ X[t]),
  xlab = expression(X[t])
)

# ---------------------------------------------------------
# Lag plot of X
# ---------------------------------------------------------

plot(
  x_lag,
  x_curr,
  pch = 16,
  cex = 0.4,
  col = rgb(0, 0, 0, 0.35),

  xlab = expression(X[t-1]),
  ylab = expression(X[t]),
  main = expression("Dependence of " ~ X[t])
)

par(mfrow = c(1, 1))

dev.off()