source("ERA5/load_data.R")

winter_hourly_gust <- data$fg10[data$season == "Winter"]

## Hill estimator
# Increasing curve suggest weibull domain of attraction instead of frechet
h_hat <- ReIns::Hill(winter_hourly_gust)

xi_hill <- h_hat$gamma

genh_hat <- ReIns::genHill(winter_hourly_gust,
     gamma = xi_hill,
     plot = TRUE
)

xi_genhill <- genh_hat$gamma

k_max <- 900

plot(1:k_max, genh_hat$gamma[1:k_max], type = "l")

tol <- 0.05 # tolerance (e.g. 2%)
k_star <- NA

for (k in 1:k_max) {

     # relative variation in the window
     rel_range <- abs((xi_genhill[k+1] - xi_genhill[k])/xi_genhill[k])

     if (rel_range < tol) {
          k_star <- k
          break
     }
}

k_star
xi_genhill[k_star]


########
## EGPD NAVEAU
########

fit_egpd_mev <- mev::fit.extgp(winter_hourly_gust,
     model = 1, init = list(kappa = 1, sigma = sd(winter_hourly_gust), xi = 0.1),
     method = "mle"
)

fit_egpd_mev$fit$mle
# fit_egpd <- mev::fit.extgp(winter_hourly_gust,
#     model = 3, init = list(kappa = 1, delta = 1, sigma = sd(winter_hourly_gust), xi = 0.1),
#     method = "mle"
# )

########
## EGPD New package
########

fit_egpd <- egpd::fitegpd(
     winter_hourly_gust,
     type = 1,
     family = "egpd"
)

summary(fit_egpd)
plot(fit_egpd)

egpd_bic <- fit_egpd$bic
egpd_aic <- fit_egpd$aic

bern_m <- c(3, 4, 6, 8, 12)
aic_compare <- numeric(length(bern_m))
bic_compare <- numeric(length(bern_m))

for (i_m in seq_along(bern_m)) {
     fit <- egpd::fitegpd(winter_hourly_gust,
          type = 1,
          method = "bernstein", bernstein.m = bern_m[i_m]
     )
     aic_compare[i_m] <- AIC(fit)
     bic_compare[i_m] <- BIC(fit)
}

data.frame(m = bern_m, AIC = round(aic_compare, 2), BIC = round(bic_compare, 2))

fit_egpd_bern <- egpd::fitegpd(
     winter_hourly_gust,
     type = 1,
     family = "egpd",
     method = "bernstein",
     bernstein.m = 10
)

summary(fit_egpd_bern)
plot(fit_egpd_bern)
BIC(fit_egpd_bern)
BIC(fit_egpd)

########
## GPD
########
library(extRemes)

# stability plot for GPD parameters
th_plot <- threshrange.plot(winter_hourly_gust, nint = 50)

head(th_plot)


u <- quantile(winter_hourly_gust, 0.96)
fit_gpd <- extRemes::fevd(winter_hourly_gust, threshold = u, type = "GP", method = "MLE")

fit_gpd$results$par

sum(winter_hourly_gust > u)

# try declustering to see if removing dependence changes the estimates

dc_runs <- extRemes::decluster(winter_hourly_gust, threshold = u, r = 3)
dc_runs

dc_intervals <- extRemes::decluster(winter_hourly_gust,
     threshold = u,
     method = "interval"
)
dc_intervals

fit_dc <- extRemes::fevd(dc_runs, threshold = u, type = "GP")
fit_dc$results$par

fit_dc_int <- extRemes::fevd(dc_intervals, threshold = u, type = "GP")
fit_dc_int$results$par

extRemes::threshrange.plot(dc_runs)
