library(copula)
source("megpd/megpd_fit.r")

gumbelC <- gumbelCopula(param = 5, dim = 2)
n <- 1000

set.seed(46)
u <- rCopula(n, gumbelC)
x <- egpd::qegpd(u[,1], kappa = 2, sigma = 1, xi = 0.3)
y <- egpd::qegpd(u[,2], kappa = 2, sigma = 1, xi = 0.1)

fit_results <- fit_megpd_bivariate(x, y, m_degree = 5, k = 10)

fit_radial <- fit_results$radial_fit
fit_delta <- fit_results$angular_fit

summary(fit_radial)
plot(fit_radial)

summary(fit_delta)

plot(fit_delta)

R <- x + y

V <- log(x / y)

plot(R, V,
     pch = 16, col = rgb(0,0,0,0.3),
     main = "Angular spread vs radial component",
     xlab = "R", ylab = "V")


plot_delta(fit_delta, R)

fitted_chi <- plot_fitted_chi(fit_megpd = fit_results, n_mc = 50)

abline(h = 2-2^(1/5), col = "red", lty = 2)

emp_chi <- texmex::chi(cbind(x,y))
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 2], col = "blue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 1], col = "lightblue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 3], col = "lightblue", lwd = 2)

plot(emp_chi)
