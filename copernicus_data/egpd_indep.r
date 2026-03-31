source("copernicus_data/load_data.R")

winter_hourly_gust <- data$fg10[data$season == "Winter"]

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

u <- quantile(winter_hourly_gust, 0.90)
fit_gpd <- extRemes::fevd(winter_hourly_gust, threshold = u, type = "GP", method = "MLE")

fit_gpd$results$par

sum(winter_hourly_gust > u)

library(extRemes)

th_plot <- threshrange.plot(winter_hourly_gust, nint = 50)

head(th_plot)

th_plot
