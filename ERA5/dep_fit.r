source("ERA5/load_data.r")

source("code/models/gpd_gumbel.r")

source("code/models/copula_markov/load_models.r")

source("code/models/copula_markov/egpd_gumbel.r")
source("code/models/copula_markov/egpd_hr.r")
source("code/models/copula_markov/egpd_tev.r")

source("code/models/copula_markov/egpd_gumbel_diag.r")

winter_hourly_gust <- data$fg10[data$season == "Winter"]
spring_hourly_gust <- data$fg10[data$season == "Spring"]

lag_df <- function(x) {
  x_lag <- dplyr::lag(x)

  df <- data.frame(x = x, x_lag = x_lag) |>
    na.omit()

  return(df)
}

winter_lag <- lag_df(winter_hourly_gust)
spring_lag <- lag_df(spring_hourly_gust)

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

### =========================
### EGPD - Gumbel Copula
### =========================
egpd_gumbel_fit <- fit_egpd_gumbel_copula(winter_hourly_gust, method = "mle")

cat("\n--- EGPD Gumbel Copula ---\n")
cat("Estimates:\n"); print(egpd_gumbel_fit$estimate)
cat("AIC:", egpd_gumbel_fit$aic, "\n")

# Density plot
plot_density(
  winter_hourly_gust,
  egpd_gumbel_fit$estimate,
  title_cop = "EGPD - Gumbel",
  bins = 40
)


### =========================
### EGPD - Hüsler–Reiss Copula
### =========================
egpd_hr_fit <- fit_egpd_hr_copula(winter_hourly_gust)

cat("\n--- EGPD Hüsler–Reiss Copula ---\n")
cat("Estimates:\n"); print(egpd_hr_fit$estimate)
cat("AIC:", egpd_hr_fit$aic, "\n")

plot_density(
  winter_hourly_gust,
  egpd_hr_fit$estimate,
  title_cop = "EGPD - HR",
  bins = 40
)


### =========================
### EGPD - tEV Copula
### =========================
egpd_tev_fit <- fit_egpd_tev_copula(winter_hourly_gust)

cat("\n--- EGPD tEV Copula ---\n")
cat("Estimates:\n"); print(egpd_tev_fit$estimate)
cat("AIC:", egpd_tev_fit$aic, "\n")

plot_density(
  winter_hourly_gust,
  egpd_tev_fit$estimate,
  title_cop = "EGPD - tEV",
  bins = 40
)

m_degree <- 6
bern_egpd_tev_fit <- fit_bernstein_tev(winter_hourly_gust, m = m_degree)
cat("\n--- bernstein EGPD tEV Copula ---\n")
cat("Estimates:\n"); print(bern_egpd_tev_fit$estimate)
cat("AIC:", bern_egpd_tev_fit$aic, "\n")

# Extract results
margin_fit <- bern_egpd_tev_fit$estimate
bern_w <- bern_egpd_tev_fit$weights

# Generate x sequence for the smooth red line
x_seq <- seq(min(winter_hourly_gust), max(winter_hourly_gust), length.out = 500)
y_fit <- egpd:::.bernstein_full_density(x_seq, margin_fit["sigma"], margin_fit["xi"], margin_fit["kappa"], bern_w, m_degree)

# Plot
hist(winter_hourly_gust, breaks = 40, prob = TRUE, col = "lightblue", border = "darkgray",
     main = "Histogram and Fitted Density", xlab = "x")
lines(x_seq, y_fit, col = "red", lwd = 2)
