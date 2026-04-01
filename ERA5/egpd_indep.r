source("ERA5/load_data.R")

winter_hourly_gust <- data$fg10[data$season == "Winter"]

## Hill estimator
# Increasing curve suggest weibull domain of attraction instead of frechet
h_hat <- ReIns::Hill(winter_hourly_gust, plot = TRUE)

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


########
## GPD
########
library(extRemes)

# stability plot for GPD parameters
th_plot <- threshrange.plot(winter_hourly_gust, nint = 50)

head(th_plot)


u <- quantile(winter_hourly_gust, 0.90)
fit_gpd <- extRemes::fevd(winter_hourly_gust, threshold = u, type = "GP", method = "MLE")

fit_gpd$results$par

sum(winter_hourly_gust > u)

# try declustering to see if removing dependence changes the estimates

dc_runs <- extRemes::decluster(winter_hourly_gust, threshold = u)
dc_runs

dc_intervals <- extRemes::decluster(winter_hourly_gust,
     threshold = u,
     method = "interval"
)
dc_intervals

fit_dc <- fevd(dc_runs, threshold = u, type = "GP")
fit_dc$results$par

fit_dc_int <- fevd(dc_intervals, threshold = u, type = "GP")
fit_dc_int$results$par
