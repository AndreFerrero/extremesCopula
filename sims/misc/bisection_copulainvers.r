library(copula)
library(microbenchmark)

# 1. THE GUMBEL H-FUNCTION (Conditional CDF: P(U <= u | V = v))
gumbel_hfunc <- function(u, v, theta) {
  # Guard against boundary values that cause log(0) or log(1)
  u <- pmax(pmin(u, 1 - 1e-12), 1e-12)
  v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
  
  ln_u <- -log(u)
  ln_v <- -log(v)
  term <- ln_u^theta + ln_v^theta
  C <- exp(-term^(1/theta))
  
  # The partial derivative dC/dv
  h <- C * (ln_v^(theta - 1)) * (term^(1/theta - 1)) / v
  return(h)
}

# 2. BISECTION SAMPLER (The Stable Alternative)
sample_bisection <- function(v_prev, theta, iterations = 20) {
  w <- runif(1)
  low <- 1e-10
  high <- 1 - 1e-10
  
  for (i in 1:iterations) {
    mid <- (low + high) / 2
    if (gumbel_hfunc(mid, v_prev, theta) < w) {
      low <- mid
    } else {
      high <- mid
    }
  }
  return((low + high) / 2)
}

# 3. COPULA PACKAGE SAMPLER (The "Gold Standard")
sample_copula_pkg <- function(u_prev, theta) {
  cop <- gumbelCopula(theta)
  # R's standard conditional sampling via Inverse Rosenblatt
  tryCatch({
    res <- cCopula(cbind(u_prev, runif(1)), copula = cop, inverse = TRUE)[2]
    return(list(val = res, error = FALSE))
  }, error = function(e) {
    return(list(val = NA, error = TRUE))
  })
}

# --- 4. THE SIMULATION STUDY ---
test_params <- c(1.1, 2, 5) # Theta range from low to very high
n_samples <- 10000
results <- data.frame()

set.seed(42)

for (theta in test_params) {
  errors_pkg <- 0
  errors_bis <- 0
  
  u_val <- runif(1) # Start value
  
  for (i in 1:n_samples) {
    # Test Package
    pkg_res <- sample_copula_pkg(u_val, theta)
    if (pkg_res$error) errors_pkg <- errors_pkg + 1
    
    # Test Bisection
    bis_res <- sample_bisection(u_val, theta)
    if (is.na(bis_res)) errors_bis <- errors_bis + 1
    
    # Update u_val for next step (Markov property)
    if (!pkg_res$error) u_val <- pkg_res$val else u_val <- runif(1)
  }
  
  results <- rbind(results, data.frame(
    Theta = theta,
    Package_Errors = errors_pkg,
    Bisection_Errors = errors_bis
  ))
}

print(results)
