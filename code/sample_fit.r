################################################################################
## PLUG-AND-PLAY EXTREME VALUE ANALYSIS ENGINE
## Models: Hill/GenHill, GPD (Declustered), EGPD (Naveau/Bernstein), Copula-Markov
################################################################################

# --- 0. PREREQUISITES ---
# Ensure these libraries are installed
library(ReIns)     # Hill/GenHill
library(extRemes)  # GPD, Declustering
library(egpd)      # EGPD fits
library(mev)       # fit.extgp
library(dplyr)

# Load your custom Copula-Markov and local logic here
source("code/models/copula_markov/egpd_gumbel.r")
source("code/models/copula_markov/egpd_hr.r")
source("code/dep_hill.r")

# --- 1. USER INPUT ---
# Replace 'YOUR_DATA_HERE' with your actual vector (e.g., winter_hourly_gust)
vec <- wave_vector

# --- 2. CONFIGURATION ---
# Tweak these parameters based on your knowledge of the data
CONFIG <- list(
  prob_threshold = 0.90,   # Quantile for GPD (0.90 - 0.95 usually)
  run_length     = 24,     # Declustering window (e.g. 12 or 24 for hourly data)
  bernstein_m    = 5,     # Degree for Bernstein polynomials (m)
  tol_genhill    = 0.05    # Tolerance for automatic k selection
)

# --- 3. MODULE: HILL & GENHILL (Semi-parametric) ---
cat("\n[1/4] Running Hill Estimators...")
h_hat <- ReIns::Hill(vec, plot = TRUE)
gh_hat <- ReIns::genHill(vec, gamma = h_hat$gamma, plot = TRUE)
dep_hill <- plot_dep_hill(vec, k_range = 1:1000)


# --- 4. MODULE: GPD & DECLUSTERING (Tail-only) ---
cat("\n[2/4] Running GPD + Declustering...")
u <- quantile(vec, CONFIG$prob_threshold)

# Declustering to handle serial dependence
# (Assumes vec is a regular time-sequence)
dc <- extRemes::decluster(vec, threshold = u, r = CONFIG$run_length)
fit_gpd <- extRemes::fevd(dc, threshold = u, type = "GP", method = "MLE")
xi_gpd <- fit_gpd$results$par["shape"]


# --- 5. MODULE: EGPD (Whole-distribution) ---
cat("\n[3/4] Running EGPD Models (Naveau & Bernstein)...")

# 1. EGPD 
fit_egpd <- egpd::fitegpd(
     vec,
     type = 1,
     family = "egpd"
)

summary(fit_egpd)
plot(fit_egpd)

# 2. Bernstein Polynomial EGPD
fit_egpd_bern <- egpd::fitegpd(vec, type = 1, method = "bernstein", 
                               bernstein.m = CONFIG$bernstein_m)
summary(fit_egpd_bern)
plot(fit_egpd_bern)

# --- 6. MODULE: DEPENDENCE & COPULAS ---
cat("\n[4/4] Running Dependence & Copula Diagnostics...")
# Create lagged dataframe
lagged_df <- data.frame(x = vec, x_lag = dplyr::lag(vec)) %>% na.omit()

# 1. Transform to Pseudo-Observations (U, V) using ECDF
# We use (rank / (n+1)) to avoid values exactly equal to 1
u_vec <- rank(lagged_df$x) / (nrow(lagged_df) + 1)
v_vec <- rank(lagged_df$x_lag) / (nrow(lagged_df) + 1)
pseudo_obs <- cbind(u_vec, v_vec)
plot(pseudo_obs, main="Pseudo-Observations (U,V)", xlab="U", ylab="V", pch=20, col=rgb(0,0,0,0.2))

# 2. Kendall's Tau (Rank-based, invariant but calculated here on pseudo-obs)
(tau <- cor(u_vec, v_vec, method = "kendall")
)

hr_egpd_fit <- fit_egpd_hr_copula(vec)
hr_egpd_fit$estimate

t_egpd_fit <- fit_egpd_t_copula(vec, optim.method = "SANN")
t_egpd_fit$estimate
t_egpd_fit$aic
sann_t_est <- t_egpd_fit$estimate
sanninit_t_egpd <- fit_egpd_t_copula(vec, optim.method = "BFGS", init = sann_t_est)

gumbel_egpd_fit <- fit_egpd_gumbel_copula(vec, optim.method = "SANN")
gumbel_egpd_fit$estimate
round(gumbel_egpd_fit$estimate, 3)
gumbel_egpd_fit$aic

# --- 7. FINAL COMPARISON & VISUALIZATION ---
summary_table <- data.frame(
  Model = c("Hill (genHill)", "GPD (Declustered)", "EGPD (Naveau)", "EGPD (Bernstein)", "Copula-EGPD"),
  Xi_Estimate = c(xi_hill_final, xi_gpd, xi_egpd_std, xi_egpd_bern, xi_copula)
)

summary_table$Tail_Type <- ifelse(summary_table$Xi_Estimate < -0.05, "Bounded (Weibull)", 
                                  ifelse(summary_table$Xi_Estimate > 0.05, "Heavy (Frechet)", "Near-Gumbel"))

cat("\n\n--- TAIL INDEX (xi) COMPARISON ---\n")
print(summary_table)

# Multi-panel Diagnostic Plot
par(mfrow=c(2,2))
# Hill Plot
if(!is.na(k_star)) {
  plot(1:k_limit, xi_gh_vec[1:k_limit], type="l", main="GenHill Plot", xlab="k", ylab="xi")
  abline(h=xi_hill_final, col="red", lty=2)
  points(k_star, xi_hill_final, col="blue", pch=19)
}
# GPD QQ Plot
plot(fit_gpd, type="qq", main="GPD QQ-Plot")
# Threshold Stability
extRemes::threshrange.plot(vec, nint = 20, main="GPD Shape Stability")
# Dependence
plot(lagged_df$x_lag, lagged_df$x, pch=20, col=rgb(0,0,0,0.2), 
     main=paste("Lag-1 Plot (Tau =", round(tau, 3), ")"), xlab="x(t-1)", ylab="x(t)")
par(mfrow=c(1,1))