library(texmex)
library(egpd)

data(wavesurge)

hist(wavesurge$wave)

acf(wavesurge$wave)
pacf(wavesurge$wave)

h_hat <- ReIns::Hill(wavesurge$wave, plot = TRUE)


fit_egpd <- egpd::fitegpd(
     wavesurge$wave,
     type = 1,
     family = "egpd",
     method = "bernstein",
     bernstein.m = 3
)

summary(fit_egpd)
plot(fit_egpd)

AIC(fit_egpd)

aic_m <- sapply(c(3, 4, 6, 8, 12), function(m) {
  AIC(fitegpd(wavesurge$wave, type = 1, method = "bernstein", bernstein.m = m))
})
data.frame(m = c(4, 6, 8, 12), AIC = round(aic_m, 2))

library(extRemes)

# stability plot for GPD parameters
th_plot <- threshrange.plot(wavesurge$wave, nint = 50)

head(th_plot)


u <- quantile(wavesurge$wave, 0.90)
fit_gpd <- extRemes::fevd(wavesurge$wave, threshold = u, type = "GP", method = "MLE")

fit_gpd$results$par

sum(wavesurge$wave > u)



dc_runs <- extRemes::decluster(wavesurge$wave, threshold = u)
dc_runs

fit_dc <- fevd(dc_runs, threshold = u, type = "GP")
fit_dc$results$par
