# ==============================================================
# Gaussian Copula Object
# ==============================================================
copula_gaussian <- list(
  name = "gaussian",
  
  # The h-function: P(U <= u | V = v)
  h_dist = function(u, v, rho) {
    # Standardize inputs to prevent Inf in qnorm
    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    
    # Gaussian conditional formula
    # x = qnorm(u), y = qnorm(v)
    # P(X <= x | Y = y) = Phi( (x - rho*y) / sqrt(1 - rho^2) )
    res <- pnorm( (qnorm(u) - rho * qnorm(v)) / sqrt(1 - rho^2) )
    return(res)
  },
  
  # The actual inverse h-function (for speed, though bisection works too)
  h_inv = function(w, v, rho) {
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    # Solving h(u|v) = w for u:
    # u = Phi( rho * Phi^-1(v) + sqrt(1 - rho^2) * Phi^-1(w) )
    u <- pnorm( rho * qnorm(v) + sqrt(1 - rho^2) * qnorm(w) )
    return(u)
  }
)