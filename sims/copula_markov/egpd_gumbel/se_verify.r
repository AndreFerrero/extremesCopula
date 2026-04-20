## =========================================================
## 0. Setup
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

n <- 5000
mc_it <- 250

set.seed(123)


## =========================================================
## 1. Storage
## =========================================================

# Copula-Markov model
par_mat_gumbel <- matrix(NA, mc_it, 4)
sd_mat_gumbel  <- matrix(NA, mc_it, 4)

colnames(par_mat_gumbel) <- c("sigma", "xi", "kappa", "theta")
colnames(sd_mat_gumbel)  <- c("sigma", "xi", "kappa", "theta")

# IID model
par_mat_iid <- matrix(NA, mc_it, 3)
sd_mat_iid  <- matrix(NA, mc_it, 3)

colnames(par_mat_iid) <- c("sigma", "xi", "kappa")
colnames(sd_mat_iid)  <- c("sigma", "xi", "kappa")


## =========================================================
## 2. Monte Carlo loop
## =========================================================

for (m in seq_len(mc_it)) {

  cat("Iteration:", m, "\n")

  ## --- simulate data ---
  data_gumbel <- egpd_gumbel_model$simulate(
    n, copula_param, margin_param
  )

  data_iid <- egpd::regpd(
    n,
    kappa = margin_param["kappa"],
    sigma = margin_param["sigma"],
    xi    = margin_param["xi"]
  )

  ## --- fit models ---
  fit_gumbel <- fit_egpd_gumbel(data_gumbel$x)

  fit_iid <- egpd::fitegpd(
    data_iid,
    type = 1,
    family = "egpd"
  )

  ## --- store results ---
  par_mat_gumbel[m, ] <- fit_gumbel$estimate
  sd_mat_gumbel[m, ]  <- fit_gumbel$sd

  par_mat_iid[m, ] <- fit_iid$estimate
  sd_mat_iid[m, ]  <- fit_iid$sd
}


## =========================================================
## 3. Save results
## =========================================================

save(
  par_mat_gumbel, sd_mat_gumbel,
  par_mat_iid, sd_mat_iid,
  file = "sims/copula_markov/egpd_gumbel/res/se_verify.Rdata"
)

load("sims/copula_markov/egpd_gumbel/res/se_verify.Rdata")
## =========================================================
## 4. Diagnostic function
## =========================================================

analyze_mc <- function(par_mat, sd_mat, true_par) {

  ## --- sanity checks ---
  if (is.null(colnames(par_mat)) || is.null(colnames(sd_mat))) {
    stop("par_mat and sd_mat must have column names")
  }

  if (!all(colnames(par_mat) %in% names(true_par))) {
    stop("Mismatch between par_mat columns and true_par names")
  }

  ## --- reorder true parameters to match matrices ---
  true_par <- true_par[colnames(par_mat)]

  ## --- compute MC statistics ---
  mc_sd   <- apply(par_mat, 2, sd)
  mean_se <- apply(sd_mat, 2, mean)

  comparison <- cbind(
    MC_SD   = mc_sd,
    Mean_SE = mean_se,
    Ratio   = mean_se / mc_sd
  )

  ## --- MC uncertainty on SD ---
  mc_se_of_sd <- mc_sd / sqrt(2 * (nrow(par_mat) - 1))
  within_mc_error <- abs(mean_se - mc_sd) < 2 * mc_se_of_sd

  ## --- coverage (now safe) ---
  coverage <- sapply(colnames(par_mat), function(param) {
    mean(
      (par_mat[, param] - 1.96 * sd_mat[, param] <= true_par[param]) &
      (par_mat[, param] + 1.96 * sd_mat[, param] >= true_par[param])
    )
  })

  ## --- bias ---
  bias <- colMeans(par_mat) - true_par

  list(
    comparison = comparison,
    coverage   = coverage,
    bias       = bias,
    within_mc_error = within_mc_error
  )
}


## =========================================================
## 5. Run diagnostics
## =========================================================

# True parameters
true_par_gumbel <- c(
  margin_param["sigma"],
  margin_param["xi"],
  margin_param["kappa"],
  theta = copula_param
)

true_par_iid <- c(
  margin_param["sigma"],
  margin_param["xi"],
  margin_param["kappa"]
)

res_gumbel <- analyze_mc(
  par_mat_gumbel,
  sd_mat_gumbel,
  true_par_gumbel
)

res_iid <- analyze_mc(
  par_mat_iid,
  sd_mat_iid,
  true_par_iid
)


## =========================================================
## 6. Print results
## =========================================================

cat("\n=== Copula-Markov model ===\n")
print(res_gumbel$comparison)
print("Coverage:")
print(res_gumbel$coverage)

cat("\n=== IID model ===\n")
print(res_iid$comparison)
print("Coverage:")
print(res_iid$coverage)

print("SD from Hessian within Monte Carlo SE?")
print(res_gumbel$within_mc_error)
print(res_iid$within_mc_error)
