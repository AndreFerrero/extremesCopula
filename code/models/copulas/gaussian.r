copula_gaussian <- list(
  name = "gaussian",
  name_param = "rho",
  
  # -------------------------------------------------
  # Joint CDF: C(u,v) = Phi_2(qnorm(u), qnorm(v); rho)
  # -------------------------------------------------
  cdf = function(u, v, rho) {
    
    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    
    x <- qnorm(u)
    y <- qnorm(v)
    
    # Vectorised evaluation
    res <- mapply(function(xi, yi) {
      mvtnorm::pmvnorm(
        lower = c(-Inf, -Inf),
        upper = c(xi, yi),
        sigma = matrix(c(1, rho, rho, 1), 2, 2)
      )[1]
    }, x, y)
    
    return(pmin(pmax(res, 0), 1))
  },
  
  # -------------------------------------------------
  # Diagonal CDF: C(u,u)
  # -------------------------------------------------
  cdf_diag = function(u, rho) {
    
    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    x <- qnorm(u)
    
    res <- sapply(x, function(xi) {
      mvtnorm::pmvnorm(
        lower = c(-Inf, -Inf),
        upper = c(xi, xi),
        sigma = matrix(c(1, rho, rho, 1), 2, 2)
      )[1]
    })
    
    return(pmin(pmax(res, 0), 1))
  },
  
  # -------------------------------------------------
  # Conditional CDF h-function
  # -------------------------------------------------
  h_dist = function(u, v, rho) {
    u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    
    res <- pnorm(
      (qnorm(u) - rho * qnorm(v)) /
      sqrt(1 - rho^2)
    )
    
    return(res)
  },
  
  # -------------------------------------------------
  # Inverse h-function
  # -------------------------------------------------
  h_inv = function(w, v, rho) {
    v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
    
    u <- pnorm(
      rho * qnorm(v) +
      sqrt(1 - rho^2) * qnorm(w)
    )
    
    return(u)
  }
)