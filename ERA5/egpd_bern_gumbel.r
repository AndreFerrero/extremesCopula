# Gumbel Copula Log-Density
# u, v are the PIT values from the marginal CDF
log_dgumbel_copula <- function(u, v, theta) {
  # Ensure u, v are in (0, 1)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  
  x <- -log(u)
  y <- -log(v)
  t <- x^theta + y^theta
  
  # Archimedean formula for Gumbel
  log_c <- -t^(1/theta) + (theta-1)*log(x) + (theta-1)*log(y) + 
           (1/theta - 2)*log(t) + log(t^(1/theta) + theta - 1) - 
           (log(u) + log(v))
  return(log_c)
}

bernstein_copula_nll <- function(theta_vec, x, m) {
  # 1. Parameter Extraction (un-transforming to ensure positivity)
  kappa <- exp(theta_vec[1])
  sigma <- exp(theta_vec[2])
  xi    <- theta_vec[3]
  # Bernstein log-weights (softmax approach as used in the package)
  alpha <- theta_vec[4:(3 + m)]
  alpha_shifted <- alpha - max(alpha)
  weights <- exp(alpha_shifted) / sum(exp(alpha_shifted))
  
  # Gumbel parameter (theta > 1)
  theta_copula <- exp(theta_vec[4 + m]) + 1 
  
  # 2. Marginal Calculations (The EGPD part)
  # Using package internal functions or manual equivalents
  # .pgpd_std is effectively 1 - (1 + xi * x/sigma)^(-1/xi)
  u_std <- egpd:::.pgpd_std(x, sigma, xi)
  u_std <- pmin(pmax(u_std, 1e-10), 1 - 1e-10)
  
  # The transition variable (uk)
  uk <- u_std^kappa
  uk <- pmin(pmax(uk, 1e-10), 1 - 1e-10)
  
  # Marginal Density
  b_dens <- egpd:::.bernstein_density(uk, weights, m)
  b_dens <- pmax(b_dens, .Machine$double.xmin)
  gpd_logd <- egpd:::.dgpd_std(x, sigma, xi, log = TRUE)
  
  # Sum of marginal log-densities
  log_f <- log(b_dens) + log(kappa) + (kappa - 1) * log(u_std) + gpd_logd
  
  # 3. Copula Calculations (The Dependence part)
  # Get the marginal CDF values (U_t)
  U_pit <- egpd:::.bernstein_cdf(uk, weights, m)
  
  # Calculate log copula density for pairs (t, t-1)
  n <- length(x)
  log_c <- log_dgumbel_copula(U_pit[2:n], U_pit[1:(n-1)], theta_copula)
  
  # 4. Total NLL
  total_ll <- sum(log_f) + sum(log_c)
  
  if (!is.finite(total_ll)) return(1e+20)
  return(-total_ll)
}

source("ERA5/load_data.r")

winter_hourly_gust <- data$fg10[data$season == "Winter"]

# Suppose your data is 'winter_hourly_gust' and m 
m_degree <- 15

# Initial values from marginal-only fit
init_marg <- fitegpd(winter_hourly_gust, method = "bernstein", bernstein.m = m_degree)

# Initial theta_vec
# [log_kappa, log_sigma, xi, alpha_1...m, log(theta_copula - 1)]
start_params <- c(
  log(init_marg$estimate["kappa"]),
  log(init_marg$estimate["sigma"]),
  init_marg$estimate["xi"],
  rep(0, m_degree), # flat weights to start
  log(2) # starting theta_copula at 3 (exp(log(2))+1)
)

# Fit the model
fit_copula_egpd_bern <- optim(
  par = start_params,
  fn = bernstein_copula_nll,
  x = winter_hourly_gust,
  m = m_degree,
  method = "Nelder-Mead",
  control = list(maxit = 5000)
)

# Extract results
res_kappa <- exp(fit_copula_egpd$par[1])
res_sigma <- exp(fit_copula_egpd$par[2])
res_xi    <- fit_copula_egpd$par[3]
res_theta <- exp(fit_copula_egpd$par[4 + m_degree]) + 1

alpha_raw <- fit_copula_egpd$par[4:(3 + m_degree)]

# We subtract the max for numerical stability (prevents overflow)
alpha_shifted <- alpha_raw - max(alpha_raw)
res_weights <- exp(alpha_shifted) / sum(exp(alpha_shifted))


# Generate x sequence for the smooth red line
x_seq <- seq(min(winter_hourly_gust), max(winter_hourly_gust), length.out = 500)
y_fit <- egpd:::.bernstein_full_density(x_seq, res_sigma, res_xi, res_kappa, res_weights, m_degree)

# Plot
hist(winter_hourly_gust, breaks = 40, prob = TRUE, col = "lightblue", border = "darkgray",
     main = "Histogram and Fitted Density", xlab = "x")
lines(x_seq, y_fit, col = "red", lwd = 2)

# Theoretical quantiles
p_points <- (1:length(winter_hourly_gust) - 0.5) / length(winter_hourly_gust)
theory_q <- egpd:::.bernstein_full_quantile(p_points, res_sigma, res_xi, res_kappa, res_weights, m_degree)
sample_q <- sort(winter_hourly_gust)

# Plot
plot(theory_q, sample_q, pch = 1, cex = 0.6, col = "darkgray",
     main = "Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
abline(0, 1, col = "red", lwd = 2)
