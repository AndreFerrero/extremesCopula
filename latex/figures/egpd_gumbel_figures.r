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
