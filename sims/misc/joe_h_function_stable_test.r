# ---------------------------------------------------------
# COMPARISON SCRIPT: Flatlining vs. Sampling
# ---------------------------------------------------------

# Your original logic
h_dist_original <- function(u, v, theta) {
  u <- pmax(pmin(u, 1 - 1e-12), 1e-12); v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
  u_base <- 1 - u; v_base <- 1 - v
  term_u <- exp(theta * log(u_base)); term_v <- exp(theta * log(v_base))
  K <- pmax(term_u + term_v - (term_u * term_v), 1e-15)
  return(exp((1/theta - 1) * log(K) + (theta - 1) * log(v_base) + log1p(-term_u)))
}

# The stable logic
h_dist_stable <- function(u, v, theta) {
  u <- pmax(pmin(u, 1 - 1e-15), 1e-15); v <- pmax(pmin(v, 1 - 1e-15), 1e-15)
  log_ratio <- theta * (log(1 - u) - log(1 - v))
  term_u <- exp(theta * log(1 - u))
  inner <- pmax(exp(log_ratio) + 1 - term_u, 1e-15)
  return(exp((1/theta - 1) * log(inner) + log1p(-term_u)))
}

sample_bisection <- function(w, v_prev, fun, theta) {
  low <- 1e-10; high <- 1 - 1e-10
  for (i in 1:20) {
    mid <- (low + high) / 2
    if (fun(mid, v_prev, theta) < w) { low <- mid } else { high <- mid }
  }
  return(mid)
}

# --- THE TEST ---
n_sim <- 500
theta <- 4
start_v <- 0.999 # Start near the top

results_orig <- numeric(n_sim)
results_stab <- numeric(n_sim)

# Set starting values
results_orig[1] <- start_v
results_stab[1] <- start_v

set.seed(42)
ws <- runif(n_sim)

for(t in 2:n_sim) {
  results_orig[t] <- sample_bisection(ws[t], results_orig[t-1], h_dist_original, theta)
  results_stab[t] <- sample_bisection(ws[t], results_stab[t-1], h_dist_stable, theta)
}

# --- ANALYZE RESULTS ---

cat("--- Original Sampler Behavior ---\n")
print(summary(results_orig))
cat("Unique values in last 400 steps:", length(unique(results_orig[100:500])), "\n")
cat("Last 5 values:", tail(results_orig, 5), "\n\n")

cat("--- Stable Sampler Behavior ---\n")
print(summary(results_stab))
cat("Unique values in last 400 steps:", length(unique(results_stab[100:500])), "\n")
cat("Last 5 values:", tail(results_stab, 5), "\n")

# Check if the Original sampler got "Stuck" at the magic value
stuck_val <- 0.9999995231
cat("\nDid the original sampler hit the boundary limit?", 
    any(abs(results_orig - stuck_val) < 1e-8), "\n")