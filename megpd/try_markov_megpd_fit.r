source("megpd/megpd_fit.r")
source("megpd/markov_megpd.r")

# load("megpd/bign_100000_burnin1000_markov_megpd_upperdelta_xi01_kappa2.Rdata")

# Prepare bivariate pairs from the Markov Chain
x_tm1 <- final_chain[-length(final_chain)]
x_t   <- final_chain[-1]

# Transform to Uniform marginals (Empirical PIT)
u_tm1 <- copula::pobs(x_tm1)
u_t   <- copula::pobs(x_t)

fit_results <- fit_megpd_bivariate(x_t, x_tm1, m_degree = 10, k = 10)

fit_radial <- fit_results$radial_fit
fit_delta <- fit_results$angular_fit

summary(fit_radial)
plot(fit_radial)

summary(fit_delta)

plot(fit_delta)

R <- x_t + x_tm1

u <- ppoints(length(R))

R_theory <- qegpd(
  u,
  sigma = 1,
  xi = xi_val,
  kappa = kappa_val,
  type = 1
)

qqplot(
  sort(R_theory),
  sort(R)
)
abline(0,1,col=2)

plot_delta(fit_delta, R)
curve(delta_f_strong_upper, from = 0, to = max(R), add = TRUE, col = "blue", lwd = 2)

V <- log(x_t / x_tm1)

plot(R, V,
     pch = 16, col = rgb(0,0,0,0.3),
     main = "Angular spread vs radial component",
     xlab = "R", ylab = "V")

fitted_chi <- plot_fitted_chi(fit_radial = fit_radial, fit_delta = fit_delta, n_mc = 200)

emp_chi <- texmex::chi(cbind(x_t,x_tm1))
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 2], col = "blue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 1], col = "lightblue", lwd = 2)
lines(seq(0.8, 0.99, length.out = 20), emp_chi$chi[80:99, 3], col = "lightblue", lwd = 2)