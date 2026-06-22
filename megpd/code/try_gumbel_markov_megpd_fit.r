source("megpd/megpd_fit.r")
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")

egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

set.seed(46)
margin_param <- c(mu = 0, kappa = 2, sigma = 1, xi = 0.1)
copula_param <- 5
n <- 10000

egpd_gumbel_data <- egpd_gumbel_model$simulate(
    n = n,
    margin_param = margin_param,
    copula_param = copula_param
)

x_t <- egpd_gumbel_data$x[-1]
x_tm1 <- egpd_gumbel_data$x[-n]

R <- x_t + x_tm1
V <- log(x_t / x_tm1)

fit_results <- fit_megpd_bivariate(x_t, x_tm1, m_degree = 10, k = 10)

fit_radial <- fit_results$radial_fit
fit_delta <- fit_results$angular_fit

summary(fit_radial)
plot(fit_radial)

summary(fit_delta)

plot(fit_delta)

plot_delta(fit_delta, R)

plot(R, V,
     pch = 16, col = rgb(0,0,0,0.3),
     main = "Angular spread vs radial component",
     xlab = "R", ylab = "V")

fitted_chi <- plot_fitted_chi(fit_megpd = fit_results, n_mc = 50)

abline(h = 2-2^(1/5), col = "red", lty = 2)

emp_chi <- texmex::chi(cbind(x_t,x_tm1))
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 2], col = "blue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 1], col = "lightblue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 3], col = "lightblue", lwd = 2)
