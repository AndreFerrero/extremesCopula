library(ExtremalDep)

data(logReturns, package = "ExtremalDep")
mm_gbp_usd <- ts(logReturns$USD, start = c(1991, 3), end = c(2014, 12), frequency = 12)
mm_gbp_jpy <- ts(logReturns$JPY, start = c(1991, 3), end = c(2014, 12), frequency = 12)

seas_usd <- stl(mm_gbp_usd, s.window = "period")
seas_jpy <- stl(mm_gbp_jpy, s.window = "period")

mm_gbp_usd_filt <- mm_gbp_usd - rowSums(seas_usd$time.series[, -3])
mm_gbp_jpy_filt <- mm_gbp_jpy - rowSums(seas_jpy$time.series[, -3])

hyperparam <- list(mu.nbinom = 3.2, var.nbinom = 4.48)
mm_gbp <- cbind(as.vector(mm_gbp_usd_filt), as.vector(mm_gbp_jpy_filt))

set.seed(123)

gbp_mar <- fExtDep.np(
    x = mm_gbp, method = "Bayesian", par10 = rep(0.1, 3), par20 = rep(0.1, 3),
    sig10 = 0.0001, sig20 = 0.0001, k0 = 5, hyperparam = hyperparam, nsim = 5e+4
)

debug(fExtDep.np(
    x = mm_gbp, method = "Bayesian", par10 = rep(0.1, 3), par20 = rep(0.1, 3),
    sig10 = 0.0001, sig20 = 0.0001, k0 = 5, hyperparam = hyperparam, nsim = 5e+4
))

diagnostics(gbp_mar)

gbp_mar_sum <- summary(object = gbp_mar, burn = 30000, plot = TRUE)

mm_gbp_rg <- apply(mm_gbp, 2, quantile, c(0.9, 0.995))