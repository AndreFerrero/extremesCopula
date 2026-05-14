source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")

egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

margin_param <- c(mu = 0, kappa = 3, sigma = 1, xi = 0.1)
copula_param <- 3
n <- 2000

egpd_gumbel_data <- egpd_gumbel_model$simulate(
     n = n,
     margin_param = margin_param,
     copula_param = copula_param
)

res <- 150

png("slides/figures/egpd_gumbel_hist_ts.png",
     width = 800, height = 400,
     res = res
)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
hist(egpd_gumbel_data$x, breaks = 30, main = "", xlab = "", ylab = "")
plot(1:n, egpd_gumbel_data$x, type = "l", main = "", xlab = "", ylab = "", col = "black")
dev.off()

png("slides/figures/egpd_gumbel_acf_pacf.png",
     width = 800, height = 400,
     res = res
)
par(mfrow = c(1, 2), mar = c(4, 2, 3, 2))
acf(egpd_gumbel_data$x, main = "ACF", ylab = "")
pacf(egpd_gumbel_data$x, main = "Partial ACF", ylab = "")
dev.off()

png("slides/figures/egpd_gumbel_lag_plot.png",
     width = 700, height = 700,
     res = res
)
par(mar = c(6.3, 4, 4, 2) + 0.1)
plot(egpd_gumbel_data$u[-1], egpd_gumbel_data$u[-n], main = "Lag Plot", xlab = expression(U[t-1]==F(x[t-1])), ylab = expression(U[t]==F(x[t])), col = "black", pch = 16, cex = 0.5)
dev.off()

png("slides/figures/egpd_gumbel_lag_plot_upper_zoom.png",
     width = 700, height = 700,
     res = res
)
par(mar = c(6.3, 4, 4, 2) + 0.1)

plot(egpd_gumbel_data$u[-1], egpd_gumbel_data$u[-n], pch = 16, cex = 0.5, main = "Lag Plot", xlab = expression(U[t-1]==F(x[t-1])), ylab = expression(U[t]==F(x[t])))
rect(0.9, 0.9, 1, 1, border = "red", lwd = 2)
dev.off()

# Upper tail
png("slides/figures/lag_upper.png", 700, 700,
     res = res
)
par(mar = c(6.3, 4, 4, 2) + 0.1)

idx <- egpd_gumbel_data$u[-n] > 0.9 & egpd_gumbel_data$u[-1] > 0.9
plot(egpd_gumbel_data$u[-1][idx], egpd_gumbel_data$u[-n][idx],
     xlim = c(0.9, 1), ylim = c(0.9, 1), xlab = expression(U[t-1]==F(x[t-1])), ylab = expression(U[t]==F(x[t])), main = "Upper Tail Lag Plot", col = "black",
     pch = 16, cex = 0.5
)
dev.off()

png("slides/figures/egpd_gumbel_lag_plot_lower_zoom.png",
     width = 700, height = 700,
     res = res
)
par(mar = c(6.3, 4, 4, 2) + 0.1)

plot(egpd_gumbel_data$u[-1], egpd_gumbel_data$u[-n], pch = 16, cex = 0.5, main = "Lag Plot", xlab = expression(U[t-1]==F(x[t-1])), ylab = expression(U[t]==F(x[t])))
rect(0, 0, 0.1, 0.1, border = "green", lwd = 2)
dev.off()

# Lower tail
png("slides/figures/lag_lower.png", 700, 700,
     res = res
)
par(mar = c(6.3, 4, 4, 2) + 0.1)

idx <- egpd_gumbel_data$u[-n] < 0.1 & egpd_gumbel_data$u[-1] < 0.1
plot(egpd_gumbel_data$u[-1][idx], egpd_gumbel_data$u[-n][idx],
     xlim = c(0, 0.1), ylim = c(0, 0.1), xlab = expression(U[t-1]==F(x[t-1])), ylab = expression(U[t]==F(x[t])), main = "Lower Tail Lag Plot", col = "black",
     pch = 16, cex = 0.5
)
dev.off()
