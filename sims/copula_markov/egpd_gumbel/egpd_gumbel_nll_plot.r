source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/margins/frechet.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")


egpd_sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)
frechet_sim_model <- make_copula_markov_model(margin_frechet, copula_gumbel, stan_mod = NULL)

true_egpd <- c(mu = 0, kappa = 6, sigma = 1, xi = 0.1)
true_frechet <- c(scale = 9, shape = 2.5)


# 1. Setup True Parameters
mod <- "EGPD"

if (mod == "EGPD") {
  sim_model <- egpd_sim_model
  true_margin <- true_egpd
  true_tail_index <- true_margin["xi"]
} else if (mod == "Frechet") {
  sim_model <- frechet_sim_model
  true_margin <- true_frechet
  true_tail_index <- round(1 / true_margin["shape"], 3)
}

true_theta <- 3.0
n_sim <- 1000

data <- sim_model$simulate(
  n = n_sim,
  copula_param = true_theta,
  margin_param = true_margin,
  seed = 123
)$x

x <- data

# Assume you have your data 'x' and initial estimates for kappa and sigma
kappa_fixed <- 10000       # fix kappa
sigma_fixed <- 1.0       # fix sigma
theta_vec_fixed <- c(log(kappa_fixed), log(sigma_fixed), NA, log(1))  # NA for xi, theta

# Grid over xi and theta
xi_vals <- seq(-0.15, 0.8, length.out = 100)          # tail parameter
theta_vals <- seq(1.01, 8, length.out = 100)        # copula dependence (>1 for Gumbel)

loglik_matrix <- matrix(NA, nrow = length(xi_vals), ncol = length(theta_vals))

# Loop over grid
for (i in seq_along(xi_vals)) {
  for (j in seq_along(theta_vals)) {
    theta_vec <- theta_vec_fixed
    theta_vec[3] <- xi_vals[i]              # xi
    theta_vec[4] <- log(theta_vals[j] - 1)  # theta transformation
    loglik_matrix[i, j] <- -egpd_gumbel_nll(theta_vec, x)  # remember NLL, so negate
  }
}

fit_res_egpd_gumbel <- fit_egpd_gumbel(x, init_par = c(log(kappa_fixed), log(sigma_fixed), 0.1, log(3.0 - 1)))

res_kappa <- exp(fit_res_egpd_gumbel$par[1])
res_sigma <- exp(fit_res_egpd_gumbel$par[2])
res_xi <- fit_res_egpd_gumbel$par[3]
res_theta <- exp(fit_res_egpd_gumbel$par[4]) + 1

max_loglik <- max(loglik_matrix, na.rm = TRUE)
plot_matrix <- loglik_matrix
plot_matrix[plot_matrix < (max_loglik - 500)] <- max_loglik - 500 

# 2. Create the filled contour plot
filled.contour(theta_vals, xi_vals, plot_matrix,
               color.palette = function(n) hcl.colors(n, "Viridis"),
               xlab = expression(theta), ylab = expression(xi),
               main = "Log-Likelihood Gradient Surface",
               plot.axes = {
                 axis(1); axis(2)
                 contour(theta_vals, xi_vals, plot_matrix, add = TRUE, col = "white", alpha = 0.5)
                 points(true_theta, true_tail_index, col = "red", pch = 19, cex = 1.5)
                 points(res_theta, res_xi, col = "blue", pch = 19, cex = 1.5)
               })
