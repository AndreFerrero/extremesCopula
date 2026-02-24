# ------------------------------
# Joe Copula Object (Stable Version)
# ------------------------------
copula_joe <- list(
  name = "joe",
  name_param = "theta",
  # -------------------------------------------------
  # Joint CDF
  # -------------------------------------------------
  cdf = function(u, v, theta){

    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)

    a <- (1 - u)^theta
    b <- (1 - v)^theta

    K <- a + b - a*b
    K <- pmax(K, 1e-15)

    C <- 1 - K^(1/theta)

    return(pmin(pmax(C,0),1))
  },

  # -------------------------------------------------
  # Diagonal CDF C(u,u)
  # -------------------------------------------------
  cdf_diag = function(u, theta){

    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)

    t <- (1-u)^theta

    K <- 2*t - t^2
    K <- pmax(K, 1e-15)

    C <- 1 - K^(1/theta)

    return(pmin(pmax(C,0),1))
  },
  # --------------------------
  # h-function: P(U <= u | V = v)
  # Derived as: K^(1/theta - 1) * (1-v)^(theta-1) * (1 - (1-u)^theta)
  # --------------------------
  h_dist = function(u, v, theta) {
    # 1. Numerical Clamping to avoid log(0) or 0^negative power
    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    
    # 2. Pre-calculate the components of the Kernel K
    # We use log1p for (1-u) if u is very small, but here standard subtraction is fine
    u_base <- 1 - u
    v_base <- 1 - v
    
    # term_u = (1-u)^theta
    # term_v = (1-v)^theta
    log_term_u <- theta * log(u_base)
    log_term_v <- theta * log(v_base)
    
    term_u <- exp(log_term_u)
    term_v <- exp(log_term_v)
    
    # 3. Calculate Kernel K = (1-u)^theta + (1-v)^theta - (1-u)^theta * (1-v)^theta
    # K can be rewritten as: 1 - (1 - (1-u)^theta) * (1 - (1-v)^theta)
    # This version is often more stable near the boundaries
    K <- term_u + term_v - (term_u * term_v)
    K <- pmax(K, 1e-15) # Safety floor for K
    
    # 4. Assemble the parts in log-space for stability
    # log(h) = (1/theta - 1)*log(K) + (theta - 1)*log(1-v) + log(1 - (1-u)^theta)
    
    log_part1 <- (1/theta - 1) * log(K)
    log_part2 <- (theta - 1) * log(v_base)
    # log1m(x) is log(1-x). More precise when term_u is very small.
    log_part3 <- log1p(-term_u) 
    
    log_h <- log_part1 + log_part2 + log_part3
    
    # 5. Return to linear space, ensuring result is in [0, 1]
    res <- exp(log_h)
    return(pmin(pmax(res, 0), 1))
  }
)