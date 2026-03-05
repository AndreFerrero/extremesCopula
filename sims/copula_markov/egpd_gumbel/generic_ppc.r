source("code/analysis/ppc_plot.R")
source("code/analysis/lambda_ppc.R")
source("code/analysis/extremal_index_ppc.R")

# general ppc
ppc_ts <- function(x, x_ppc) {
  N <- length(x)
  df <- data.frame(
    y = as.numeric(x),
    year = as.numeric(time(x)),
    time = 1:N
  )

  preds <- cbind(
    Estimate = colMeans(x_ppc),
    Q5 = apply(x_ppc, 2, quantile, probs = 0.05),
    Q95 = apply(x_ppc, 2, quantile, probs = 0.95)
  )

  ggplot(cbind(df, preds), aes(x = year, y = Estimate)) +
    geom_smooth(aes(ymin = Q5, ymax = Q95), stat = "identity", linewidth = 0.5) +
    geom_point(aes(y = y)) + theme_minimal()
}

ppc_ts(gaussian_sim$x, gaussian_fit$ppc)

ppc_dens_overlay(gaussian_sim$x, gaussian_fit$ppc[1:100, ]) +
  coord_cartesian(xlim = c(0, 30))
ppc_dens_overlay(gaussian_sim$x, ggumbel_fit$ppc[1:100, ]) +
  coord_cartesian(xlim = c(0, 30))
ppc_dens_overlay(gaussian_sim$x, gjoe_fit$ppc[1:100, ]) +
  coord_cartesian(xlim = c(0, 30))

ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "max")
ppc_stat(gaussian_sim$x, ggumbel_fit$ppc, stat = "max")
ppc_stat(gaussian_sim$x, gjoe_fit$ppc, stat = "max")

ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "min")
ppc_stat(gaussian_sim$x, ggumbel_fit$ppc, stat = "min")
ppc_stat(gaussian_sim$x, gjoe_fit$ppc, stat = "min")

q90 <- function(x) quantile(x, probs = 0.90)
ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "q90")
ppc_stat(gaussian_sim$x, ggumbel_fit$ppc, stat = "q90")
ppc_stat(gaussian_sim$x, gjoe_fit$ppc, stat = "q90")

q97 <- function(x) quantile(x, probs = 0.97)
ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "q97")
ppc_stat(gaussian_sim$x, ggumbel_fit$ppc, stat = "q97")
ppc_stat(gaussian_sim$x, gjoe_fit$ppc, stat = "q97")

# Extremal Index
b_vals <- seq(10, 200, 10)
b_choice <- exdex::choose_b(gaussian_sim$x, b_vals)
plot(b_choice)

(ppc_theta_gauss  <- ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "stat_theta_b40"))
(ppc_theta_gumbel <- ppc_stat(gaussian_sim$x, ggumbel_fit$ppc,  stat = "stat_theta_b40"))

# lambda evaluation
q_prob <- 0.95
z <- quantile(gaussian_sim$x, probs = q_prob)

ggumbel_lambda_posterior <- gumbel_model$lambda_z_posterior(z, ggumbel_fit)
gaussian_lambda_posterior <- gaussian_model$lambda_z_posterior(z, gaussian_fit)

par(mfrow = c(1, 2))
hist(ggumbel_lambda_posterior)
abline(v = lambda_z_empirical(gaussian_sim$x, prob = q_prob), col = "red")

hist(gaussian_lambda_posterior)
abline(v = lambda_z_empirical(gaussian_sim$x, prob = q_prob), col = "red")
par(mfrow = c(1, 1))

ggumbel_lambda_ppc <- ppc_lambda(gaussian_sim$x, ggumbel_fit$ppc, prob = 0.9)
gaussian_lambda_ppc <- ppc_lambda(gaussian_sim$x, gaussian_fit$ppc, prob = 0.9)

ppc_plot(ggumbel_lambda_ppc)
ppc_plot(gaussian_lambda_ppc)

ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "lambda_85")
ppc_stat(gaussian_sim$x, ggumbel_fit$ppc, stat = "lambda_85")

ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "lambda_90")
ppc_stat(gaussian_sim$x, ggumbel_fit$ppc, stat = "lambda_90")

ppc_stat(gaussian_sim$x, gaussian_fit$ppc, stat = "lambda_95")
ppc_stat(gaussian_sim$x, ggumbel_fit$ppc, stat = "lambda_95")

ppc_stat_2d(gaussian_sim$x, gaussian_fit$ppc, stat = c("lambda_85", "lambda_90"))
ppc_stat_2d(gaussian_sim$x, ggumbel_fit$ppc, stat = c("lambda_85", "lambda_90"))

print(ggumbel_fit$fit, pars = c("mu", "kappa", "sigma", "xi", "theta"))

print(gaussian_fit$fit, pars = c("mu", "kappa", "sigma", "xi", "rho"))
mcmc_trace(gaussian_fit$fit, pars = c("mu", "kappa", "sigma", "xi", "rho"))

jgumbel_lambda_ppc <- ppc_lambda(joe_sim$x, jgumbel_fit$ppc, prob = 0.95)
joe_lambda_ppc <- ppc_lambda(joe_sim$x, joe_fit$ppc, prob = 0.95)

ppc_plot(jgumbel_lambda_ppc)
ppc_plot(joe_lambda_ppc)

par(mfrow = c(1, 2))
gumbel_lambda_posterior <- gumbel_model$lambda_z_posterior(z, jgumbel_fit)
hist(gumbel_lambda_posterior)
abline(v = lambda_z_empirical(joe_sim$x, prob = q_prob))

joe_lambda_posterior <- joe_model$lambda_z_posterior(quantile(joe_sim$x, probs = q_prob), joe_fit)
hist(joe_lambda_posterior)
abline(v = lambda_z_empirical(joe_sim$x, prob = q_prob), col = "red")
par(mfrow = c(1, 1))