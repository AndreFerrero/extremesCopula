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


# --- 1. SETUP ---
theta_values <- c(1.1, 1.5, 2.0, 3.0, 5.0, 8.0)
T_length     <- 1000 
n_mc         <- 50    # Number of Monte Carlo iterations
q_thresh     <- 0.98   
set.seed(42)

# Theoretical Lambda_U formula
theo_lambda_u <- function(theta) 2 - 2^(1/theta)

# Storage for MC results
mc_results <- list()

# --- 2. MONTE CARLO LOOP ---
for (th in theta_values) {
  
  # Temporary storage for this specific theta's 50 runs
  tau_sims    <- numeric(n_mc)
  lambda_sims <- numeric(n_mc)
  
  for (m in 1:n_mc) {
    # Generate Markov sequence
    u <- numeric(T_length)
    u[1] <- runif(1)
    for(t in 2:T_length) {
      u[t] <- sample_bisection(u[t-1], th, iterations = 25)
    }
    
    # Calculate empirical metrics for this specific run
    tau_sims[m]    <- cor(u[-T_length], u[-1], method = "kendall")
    
    excess         <- u > q_thresh
    lambda_sims[m] <- mean(excess[-1] & excess[-T_length]) / mean(excess)
  }
  
  # Calculate Theoreticals
  tau_theo    <- 1 - 1/th
  lambda_theo <- theo_lambda_u(th)
  
  # Store Aggregated Statistics
  mc_results[[as.character(th)]] <- data.frame(
    Theta       = th,
    Theo_Tau    = round(tau_theo, 4),
    Mean_Tau    = round(mean(tau_sims), 4),
    SD_Tau      = round(sd(tau_sims), 4),
    Tau_Bias    = round(mean(tau_sims) - tau_theo, 4),
    Theo_Lambda = round(lambda_theo, 4),
    Mean_Lambda = round(mean(lambda_sims, na.rm = TRUE), 4),
    SD_Lambda   = round(sd(lambda_sims, na.rm = TRUE), 4),
    Lambda_Bias = round(mean(lambda_sims, na.rm = TRUE) - lambda_theo, 4)
  )
}

# --- 3. FINAL OUTPUT TABLE ---
results_table <- do.call(rbind, mc_results)
print(results_table, row.names = FALSE)

T_plot <- 5000
u2 <- numeric(T_plot)
theta_plot <- 3

u2[1] <- runif(1)
for(t in 2:T_plot) {
u2[t] <- sample_bisection(u2[t-1], theta_plot, iterations = 25)
}

plot(u2[-T_plot], u2[-1], 
     pch = 16, col = rgb(0, 0, 1, 0.2),
     xlab = expression(U[t-1]), ylab = expression(U[t]),
     main = paste0("Pairs plot, theta = ", theta_plot))
abline(0, 1, col = "red", lty = 2)

w_t <- gumbel_hfunc(u2[-T_plot], u2[-1], theta_plot)

hist(w_t)