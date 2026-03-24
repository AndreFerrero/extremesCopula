source("winter2016/load_data.r")

# Make the plots that can be seen in figure 2 of the data section 5.1 in the paper
# The time-series plot is the first one and then the plot of consecutive pairs
plot.dat <- JJA.data$TX * 0.1
plot.dat[plot.dat < -100] <- NA
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))
plot(plot.dat,
  ylim = c(10, 40), type = "l", xaxt = "n", xlab = "Year", ylab = "Temperature",
  cex.axis = 1.2, lwd = 2, cex.lab = 1.2
)
plotseq <- seq(368, 6164, by = 920)
labs <- c(1950, 1960, 1970, 1980, 1990, 2000, 2010)
axis(1, at = plotseq, labels = labs, cex.axis = 1.2, cex.lab = 1.2)
box(lwd = 3)

#pcaf
pacf(j.data[, 1], main = "", lwd = 2, cex.axis = 1.2, cex.lab = 1.2)


#pairs plot
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))

box(lwd = 3)
plot(j.data,
  xlab = expression(X[t]), ylab = expression(X[t + 1]), cex.axis = 1.2,
  lwd = 2, cex.lab = 1.2
)
box(lwd = 3)

# Estimate F once on all observed values
all_vals <- c(j.data[, 1], j.data[, 2])
n <- length(all_vals)
F_hat <- ecdf(all_vals)

# Apply the same F_hat to both columns
u1 <- F_hat(j.data[, 1]) * n / (n + 1)  # scaled to avoid exact 1
u2 <- F_hat(j.data[, 2]) * n / (n + 1)
pobs_pairs <- cbind(u1, u2)

plot(pobs_pairs,
  xlab = expression(hat(U)[t]), ylab = expression(hat(U)[t + 1]),
  pch = 16, col = rgb(0, 0, 0, 0.15),   # semi-transparent points
  cex = 0.6, cex.axis = 1.2, cex.lab = 1.2
)
box(lwd = 3)

# ACF and PACF for the whole data
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))
acf((JJA.data$TX * 0.1)[(JJA.data$TX * 0.1) > -100], main = "", lwd = 2, cex.axis = 1.2, cex.lab = 1.2, xlim = c(0, 30))
box(lwd = 3)
pacf((JJA.data$TX * 0.1)[(JJA.data$TX * 0.1) > -100], main = "", lwd = 2, cex.axis = 1.2, cex.lab = 1.2, xlim = c(0, 30))
box(lwd = 3)

# Change the data onto Gaussian margins and replot the PACF
# to show that the margins are not having any effect on things,
# Any depednence is being captured by the copula
ecdf.fun <- ecdf((JJA.data$TX * 0.1)[(JJA.data$TX * 0.1) > -100])
dat.vals <- sapply(X = (JJA.data$TX * 0.1)[(JJA.data$TX * 0.1) > -100], FUN = ecdf.fun)
dat.vals[dat.vals == 0] <- 1e-10
dat.vals[dat.vals == 1] <- 1 - 1e-10
norm.vals <- qnorm(dat.vals)
pacf(norm.vals, main = "", lwd = 2, cex.axis = 1.2, cex.lab = 1.2, xlim = c(0, 30))
box(lwd = 3)


