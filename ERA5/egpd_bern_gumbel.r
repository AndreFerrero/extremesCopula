source("ERA5/load_data.r")
source("code/models/copula_markov/bernstein_egpd_gumbel.r")
library(egpd)

winter_hourly_gust <- data$fg10[data$season == "Winter"]

# Suppose your data is 'winter_hourly_gust' and m 
m_degree <- 5
fit_copula_egpd_bern <- fitegpd_berstein_gumbel(winter_hourly_gust, m = m_degree)

# Extract results
res_kappa <- exp(fit_copula_egpd_bern$par[1])
res_sigma <- exp(fit_copula_egpd_bern$par[2])
res_xi    <- fit_copula_egpd_bern$par[3]
res_theta <- exp(fit_copula_egpd_bern$par[4 + m_degree]) + 1

alpha_raw <- fit_copula_egpd_bern$par[4:(3 + m_degree)]

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
