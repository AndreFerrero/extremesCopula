source("winter2016/load_data.r")

# fit mle indep egpd
library(mev)
mle_indep_egpd <- fit.extgp(data-min(data), init = c(1, 1, 0.1), method = "mle", confint = TRUE, R = 20)
round(mle_indep_egpd$fit$mle, 2)

p <- 0.9
u <- quantile(data, p)
fit_gpd <- extRemes::fevd(data, threshold = u, type = "GP", method = "MLE")

plot(fit_gpd)
summary(fit_gpd)