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
  # 1. Clamping for stability
  u <- pmax(pmin(u, 1 - 1e-15), 1e-15)
  v <- pmax(pmin(v, 1 - 1e-15), 1e-15)
  
  u_base <- 1 - u
  v_base <- 1 - v
  
  # 2. Use the ratio-based stable formula
  # Let ratio = ((1-u)^theta) / ((1-v)^theta)
  # log_ratio = theta * (log(1-u) - log(1-v))
  log_ratio <- theta * (log(u_base) - log(v_base))
  term_ratio <- exp(log_ratio)
  
  # Inner term: [ (1-u)^theta / (1-v)^theta + 1 - (1-u)^theta ]
  # Which is: term_ratio + 1 - (1-u)^theta
  term_u <- exp(theta * log(u_base))
  inner <- term_ratio + 1 - term_u
  
  # 3. Handle potential precision issues with inner
  inner <- pmax(inner, 1e-15)
  
  # 4. Final assembly: inner^(1/theta - 1) * (1 - (1-u)^theta)
  log_h <- (1/theta - 1) * log(inner) + log1p(-term_u)
  
  res <- exp(log_h)
  return(pmin(pmax(res, 0), 1))
}
)