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

mb <- microbenchmark(
  Bisection = sample_bisection(0.5, 2, iterations = 20),
  R_Package = cCopula(cbind(0.5, 0.5), gumbelCopula(2), inverse = TRUE),
  times = 1000
)
print(mb)

# 1. Generate samples
v_fixed <- 0.5
theta_fixed <- 2

# We need to capture the 'w' used for bisection to compare
# Let's just generate the u's and then 'un-transform' them
u_bis <- replicate(10000, sample_bisection(v_fixed, theta_fixed, iterations = 20))
u_pkg <- replicate(10000, sample_copula_pkg(v_fixed, theta_fixed)$val)

# 2. Apply the h-function to 'un-transform' them back to the probability space
w_bis <- gumbel_hfunc(u_bis, v_fixed, theta_fixed)
w_pkg <- gumbel_hfunc(u_pkg, v_fixed, theta_fixed)

# 3. NOW run the KS test on the W values
ks_bis <- ks.test(w_bis, "punif")
ks_pkg <- ks.test(w_pkg, "punif")

print(ks_bis)
print(ks_pkg)

# Generate a long Markov chain using bisection
n_long <- 5000
theta_test <- 3
u_series <- numeric(n_long)
u_series[1] <- runif(1)

for(t in 2:n_long) {
  u_series[t] <- sample_bisection(u_series[t-1], theta_test, iterations = 25)
}

# Plot ACF of the ranks
acf(u_series, main = paste("Rank ACF (Theta =", theta_test, ")"))

# Estimate empirical tau from the simulated sequence
emp_tau <- cor(u_series[-n_long], u_series[-1], method = "kendall")
theo_tau <- 1 - 1/theta_test

cat("Theoretical Tau:", theo_tau, "\n")
cat("Empirical Tau:  ", emp_tau, "\n")