## =========================================================
## 0. Setup (Unchanged)
## =========================================================
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")

egpd_gumbel_model <- make_copula_markov_model(
  margin_egpd, copula_gumbel, stan_mod = NULL
)

margin_param <- c(mu = 0, kappa = 2, sigma = 1, xi = 0.3)
copula_param <- 3

n <- 1000
mc_it <- 50
B <- 50  # Number of Bootstrap iterations per MC step

set.seed(123)

## =========================================================
## 1. Storage (Added Bootstrap SD)
## =========================================================
par_mat_gumbel  <- matrix(NA, mc_it, 4)
sd_mat_gumbel   <- matrix(NA, mc_it, 4) # Hessian-based
boot_sd_gumbel  <- matrix(NA, mc_it, 4) # Bootstrap-based

colnames(par_mat_gumbel) <- colnames(sd_mat_gumbel) <- 
  colnames(boot_sd_gumbel) <- c("sigma", "xi", "kappa", "theta")

## =========================================================
## 2. Monte Carlo loop
## =========================================================

for (m in seq_len(mc_it)) {
  
  cat("Iteration:", m, " ")
  
  # --- 2.1 Simulate "Real" Data ---
  data_gumbel <- egpd_gumbel_model$simulate(n, copula_param, margin_param)
  
  # --- 2.2 Fit Original Model ---
  fit_gumbel <- fit_egpd_gumbel(data_gumbel$x)
  
  # Store Estimates and Hessian SD
  par_mat_gumbel[m, ] <- fit_gumbel$estimate
  sd_mat_gumbel[m, ]  <- fit_gumbel$sd
  
  # --- 2.3 Parametric Bootstrap ---
  # We use the estimates from fit_gumbel to generate B new datasets
  cat("[Running Bootstrap...]")
  
  boot_estimates <- matrix(NA, B, 4)
  
  # Extract estimates for the simulator
  curr_margin <- c(
    mu    = 0, 
    fit_gumbel$estimate["sigma"],
    fit_gumbel$estimate["xi"],
    fit_gumbel$estimate["kappa"]
  )
  curr_copula <- fit_gumbel$estimate["theta"]
  
  for (b in seq_len(B)) {
    # Generate data from the ESTIMATED model
    boot_data <- egpd_gumbel_model$simulate(n, curr_copula, curr_margin)
    
    # Fit model to bootstrap data (hessian=FALSE to save time)
    # We only need the point estimates for the bootstrap SE
    fit_b <- try(fit_egpd_gumbel(boot_data$x, hessian = FALSE), silent = TRUE)
    
    if (!inherits(fit_b, "try-error") && fit_b$convergence == 0) {
      boot_estimates[b, ] <- fit_b$estimate
    }
  }
  
  # The Bootstrap SE is the SD of the bootstrap estimates
  boot_sd_gumbel[m, ] <- apply(boot_estimates, 2, sd, na.rm = TRUE)
  
  cat(" Done.\n")
}

save(
  par_mat_gumbel, sd_mat_gumbel,
  boot_sd_gumbel,
  file = "sims/copula_markov/egpd_gumbel/res/se_boot_verify.Rdata"
)

load("sims/copula_markov/egpd_gumbel/res/se_boot_verify.Rdata")

## =========================================================
## 3. Updated Diagnostic (Comparing Hessian vs Bootstrap)
## =========================================================

analyze_mc_extended <- function(par_mat, sd_hessian, sd_boot, true_par) {
  
  true_par <- true_par[colnames(par_mat)]
  mc_sd    <- apply(par_mat, 2, sd)
  
  mean_hess_se <- apply(sd_hessian, 2, mean, na.rm = TRUE)
  mean_boot_se <- apply(sd_boot, 2, mean, na.rm = TRUE)
  
  comparison <- cbind(
    True_MC_SD   = mc_sd,
    Hessian_SE   = mean_hess_se,
    Bootstrap_SE = mean_boot_se,
    Hess_Ratio   = mean_hess_se / mc_sd,
    Boot_Ratio   = mean_boot_se / mc_sd
  )
  
  # Coverage using Bootstrap SE
  # (Assuming Normality of the bootstrap distribution for simplicity)
  coverage_boot <- sapply(colnames(par_mat), function(p) {
    mean(
      (par_mat[, p] - 1.96 * sd_boot[, p] <= true_par[p]) &
      (par_mat[, p] + 1.96 * sd_boot[, p] >= true_par[p]),
      na.rm = TRUE
    )
  })
  
  list(comparison = comparison, coverage_boot = coverage_boot)
}

# Run diagnostics
true_par_vals <- c(margin_param["sigma"], margin_param["xi"], 
                   margin_param["kappa"], theta = copula_param)

res <- analyze_mc_extended(par_mat_gumbel, sd_mat_gumbel, boot_sd_gumbel, true_par_vals)

print(res$comparison)
print("Bootstrap Coverage:")
print(res$coverage_boot)
