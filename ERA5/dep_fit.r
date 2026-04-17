source("ERA5/load_data.r")

source("code/models/gpd_gumbel.r")

source("code/models/copula_markov/load_models.r")

source("code/models/copula_markov/egpd_gumbel.r")

source("code/models/copula_markov/egpd_gumbel_diag.r")

winter_hourly_gust <- data$fg10[data$season == "Winter"]

lag_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) |>
    na.omit()

  return(df)
}

winter_lag <- lag_df(winter_hourly_gust)

tau <- cor(winter_lag$x, winter_lag$x_lag, method = "kendall")

theta_tau <- 1 / (1 - tau)
theta_tau

# dependence statistics : Strong evidence for asymptotic dependence
chi <- texmex::chi(as.matrix(winter_lag))
par(mfrow = c(1, 2))
plot(chi, show = c("Chi" = TRUE, "ChiBar" = TRUE))
par(mfrow = c(1, 1))


###
# GUMBEL + GPD - Censored likelihood - Winter2016 implementation
###
u.thresh <- quantile(winter_hourly_gust, probs = 0.90) # Modelling threshold u on uniform scale

cens_gpd_gumbel <- FitGpd(dat = winter_lag, u = u.thresh)

gamma <- cens_gpd_gumbel$par[1]
sig.u <- cens_gpd_gumbel$par[2]
xi <- cens_gpd_gumbel$par[3]
lambda.u <- cens_gpd_gumbel$par[4]

gamma
1 / gamma
sig.u
xi
lambda.u

#######
# Frequentist Gumbel + EGPD Power
######

egpd_gumbel_fit <- fit_egpd_gumbel(winter_hourly_gust)
egpd_gumbel_fit$estimate
egpd_gumbel_bic <- egpd_gumbel_fit$bic
egpd_gumbel_aic <- egpd_gumbel_fit$aic

wind_pit <- get_pit_values(winter_hourly_gust, egpd_gumbel_fit$estimate, copula_gumbel$h_dist)

plot_diag_plots(wind_pit)

plot_density(winter_hourly_gust, egpd_gumbel_fit$estimate, bins = 40)
