source("code/packages.R")
source("code/models/margins/egp.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/copula_markov_model.R")

results_path <- here("sims/copula_markov/res/bisection_copula_invers.rds")

# --- 1. SETUP ---
theta_values <- c(1.1, 1.5, 2.0, 3.0, 5.0, 8.0)
n <- 2000
n_mc <- 25 # Number of Monte Carlo iterations
q_thresh <- 0.95
set.seed(42)

# Theoretical Lambda_U formula
theo_lambda_u <- function(theta) 2 - 2^(1 / theta)

# Storage for MC results
mc_results <- list()

# --- 2. MONTE CARLO LOOP ---
for (th in theta_values) {
  # Temporary storage for this specific theta's 50 runs
  tau_sims <- numeric(n_mc)
  lambda_sims <- numeric(n_mc)

  for (m in 1:n_mc) {
    # Generate Markov sequence
    u <- numeric(n)
    u[1] <- runif(1)
    for (t in 2:n) {
      u[t] <- sample_bisection(runif(1), u[t - 1],
        copula_gumbel$h_dist,
        th,
        iterations = 25
      )
    }

    # Calculate empirical metrics for this specific run
    tau_sims[m] <- cor(u[-n], u[-1], method = "kendall")

    excess <- u > q_thresh
    lambda_sims[m] <- mean(excess[-1] & excess[-n]) / mean(excess)
  }

  # Calculate Theoreticals
  tau_theo <- 1 - 1 / th
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

set.seed(123)
mb <- microbenchmark(
  Bisection = sample_bisection(
    w = 0.5, v_prev = 0.5,
    fun = copula_gumbel$h_dist, copula_param = 2, iterations = 20
  ),
  R_Package = cCopula(cbind(0.5, 0.5), gumbelCopula(2), inverse = TRUE),
  times = 1000
)

# --- 3. SAVE ---
saveRDS(
  list(mc_results = mc_results, benchmark = mb),
  results_path
)
