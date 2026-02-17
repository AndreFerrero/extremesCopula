source("libs/packages.R")

# ------------------------------
# Clayton copula functions
# ------------------------------

copula_clayton <- list(

  name = "clayton",

  # P(U <= u | V = v) = v^(-theta-1) * (u^-theta + v^-theta - 1)^(-1/theta - 1)
  h_dist = function(u, v, theta) {
    if (abs(theta) < 1e-8) return(u) # Independence case
    
    # Numerical clamping
    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    
    term <- u^(-theta) + v^(-theta) - 1
    # Check if the term inside is valid
    if (any(term <= 0)) return(rep(0, length(u))) 
    
    res <- (v^(-theta - 1)) * (term^(-1/theta - 1))
    return(pmin(pmax(res, 0), 1))
  },

  # --- 2. Inverse h-function (for conditional sampling) ---
  # Solving h(u|v) = w for u leads to:
  # u = [ (w^(-theta/(theta+1)) - 1) * v^(-theta) + 1 ]^(-1/theta)
  h_inv = function(w, v, theta) {
    if (abs(theta) < 1e-8) return(w) # Independence case
    
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    w <- pmax(pmin(w, 1 - 1e-12), 1e-12)
    
    # Analytical solution
    exponent <- -theta / (theta + 1)
    u <- ( (w^exponent - 1) * v^(-theta) + 1 )^(-1 / theta)
    
    return(pmin(pmax(u, 0), 1))
  }
) 
