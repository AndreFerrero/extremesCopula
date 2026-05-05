png("slides/figures/egpd_gumbel_hist_ts.png", width = 800, height = 400)
par(mfrow = c(1, 2))
hist(egpd_gumbel_data$x, breaks = 30, main = "EGPD-Gumbel Data", xlab = "Value")
plot(1:n, egpd_gumbel_data$x, type = "l", main = "Time Series Plot", xlab = "Time", ylab = "Value", col = "black")
dev.off()

png("slides/figures/egpd_gumbel_acf_pacf.png", width = 800, height = 400)
par(mfrow = c(1, 2))
acf(egpd_gumbel_data$x, main = "ACF")
pacf(egpd_gumbel_data$x, main = "PACF")
dev.off()

png("slides/figures/egpd_gumbel_lag_plot.png", width = 350, height = 350)
plot(egpd_gumbel_data$u[-1], egpd_gumbel_data$u[-n], main = "Lag Plot", xlab = "U_(t-1)", ylab = "U_t", col = "black", pch = 16, cex = 0.5)
dev.off()

png("slides/figures/egpd_gumbel_lag_plot_upper_zoom.png", width = 350, height = 350)
plot(egpd_gumbel_data$u[-1], egpd_gumbel_data$u[-n], pch=16, cex=0.5, main = "Lag Plot", xlab="U_(t-1)", ylab="U_t")
rect(0.9, 0.9, 1, 1, border="red", lwd=2)
dev.off()

# Upper tail
png("slides/figures/lag_upper.png", 350, 350)
idx <- egpd_gumbel_data$u[-n] > 0.9 & egpd_gumbel_data$u[-1] > 0.9
plot(egpd_gumbel_data$u[-1][idx], egpd_gumbel_data$u[-n][idx], xlim=c(0.9,1), ylim=c(0.9,1), xlab="U_(t-1)", ylab="U_t", main="Upper Tail Lag Plot", col = "black",
     pch=16, cex=0.5)
dev.off()

png("slides/figures/egpd_gumbel_lag_plot_lower_zoom.png", width = 350, height = 350)
plot(egpd_gumbel_data$u[-1], egpd_gumbel_data$u[-n], pch=16, cex=0.5, main = "Lag Plot", xlab="U_(t-1)", ylab="U_t")
rect(0, 0, 0.1, 0.1, border="green", lwd=2)
dev.off()

# Lower tail
png("slides/figures/lag_lower.png", 350, 350)
idx <- egpd_gumbel_data$u[-n] < 0.1 & egpd_gumbel_data$u[-1] < 0.1
plot(egpd_gumbel_data$u[-1][idx], egpd_gumbel_data$u[-n][idx], xlim=c(0,0.1), ylim=c(0,0.1), xlab="U_(t-1)", ylab="U_t", main="Lower Tail Lag Plot", col = "black",
     pch=16, cex=0.5)
dev.off()