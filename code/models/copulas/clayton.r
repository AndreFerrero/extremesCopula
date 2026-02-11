source("libs/packages.R")

# ------------------------------
# Clayton copula functions
# ------------------------------

copula_clayton <- list(

  name = "clayton",

  # --------------------------
  # 1. Simulate uniforms using latent gamma variable
  # --------------------------
  # Clayton copula CDF: C(u1,...,un) = (sum(u_i^-theta - 1) + 1)^(-1/theta)
  # Latent variable V ~ Gamma(1/theta, 1)
  # Then U_i = (1 + E_i / V)^(-1/theta), E_i ~ Exp(1)
  simulate = function(theta, n) {

  # --- Independence case ---
  if (abs(theta) < 1e-8) {
    return(runif(n))
  }

  # --- Invalid parameter ---
  if (theta < 0) return(NULL)

  # --- Clayton Archimedean simulation ---
  V <- rgamma(n = 1, shape = 1 / theta, rate = 1)
  if (!is.finite(V) || V <= 0) return(NULL)

  E <- rexp(n)
  (1 + E / V)^(-1 / theta)
},

  # --------------------------
  # 2. Log-density of the Clayton copula
  # --------------------------
  lpdf = function(u, theta) {
    copula::dCopula(u, copula::claytonCopula(theta, dim = length(u)), log = TRUE)
  },

  # --------------------------
  # 3. Log-prior for the copula parameter
  # --------------------------
  log_prior = function(theta, a = 2, b = 1) {
    if (theta <= 0) return(-Inf)
    dgamma(theta, shape = a, rate = b, log = TRUE)
  },

  diag = function(u, theta) {
    n <- length(u)
    
    if (theta < 0) return(NA)             # invalid parameter
    if (theta == 0) return(u^n)           # independence case
    
    (n * (u^(-theta) - 1) + 1)^(-1/theta)
  },

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
