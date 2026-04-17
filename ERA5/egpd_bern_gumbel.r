source("ERA5/load_data.r")
source("code/models/copula_markov/bernstein_egpd_gumbel.r")
library(egpd)

winter_hourly_gust <- data$fg10[data$season == "Winter"]

# Suppose your data is 'winter_hourly_gust' and m 
m_degree <- 10
fit_copula_egpd_bern <- fit_egpd_bernstein_gumbel(winter_hourly_gust, m = m_degree)

# Extract results
margin_fit <- fit_copula_egpd_bern$estimate
bern_w <- fit_copula_egpd_bern$weights

# Generate x sequence for the smooth red line
x_seq <- seq(min(winter_hourly_gust), max(winter_hourly_gust), length.out = 500)
y_fit <- egpd:::.bernstein_full_density(x_seq, margin_fit["sigma"], margin_fit["xi"], margin_fit["kappa"], bern_w, m_degree)

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
