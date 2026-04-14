source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")
source("code/models/copula_markov/egpd_gumbel.r")

sim_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

true_margin <- c(mu = 0, kappa = 6, sigma = 1, xi = 0.1)

true_theta <- 2.0
n_sim <- 5000

data <- sim_model$simulate(
  n = n_sim,
  copula_param = true_theta,
  margin_param = true_margin,
  seed = 123
)$x

x <- data

fit_res_egpd_gumbel <- fit_egpd_gumbel(x)

res_kappa <- exp(fit_res_egpd_gumbel$par[1])
res_sigma <- exp(fit_res_egpd_gumbel$par[2])
res_xi <- fit_res_egpd_gumbel$par[3]
res_theta <- exp(fit_res_egpd_gumbel$par[4]) + 1

# LIKELIHOOD VISUALISATION

# Assume you have your data 'x' and initial estimates for kappa and sigma
kappa_fixed <- true_margin["kappa"]       # fix kappa
sigma_fixed <- true_margin["sigma"]       # fix sigma
theta_vec_fixed <- c(log(kappa_fixed), log(sigma_fixed), NA, log(1))  # NA for xi, theta

# Grid over xi and theta
xi_vals <- seq(-0.15, 0.8, length.out = 100)          # tail parameter
theta_vals <- seq(1.01, 8, length.out = 100)        # copula dependence (>1 for Gumbel)

loglik_theta_xi <- matrix(NA, nrow = length(xi_vals), ncol = length(theta_vals))

# Loop over grid
for (i in seq_along(xi_vals)) {
  for (j in seq_along(theta_vals)) {
    theta_vec <- theta_vec_fixed
    theta_vec[3] <- xi_vals[i]              # xi
    theta_vec[4] <- log(theta_vals[j] - 1)  # theta transformation
    loglik_theta_xi[i, j] <- -egpd_gumbel_nll(theta_vec, x)  # remember NLL, so negate
  }
}

max_loglik_theta_xi <- max(loglik_theta_xi, na.rm = TRUE)
plot_matrix_theta_xi <- loglik_theta_xi
plot_matrix_theta_xi[plot_matrix_theta_xi < (max_loglik_theta_xi - 500)] <- max_loglik_theta_xi - 500 

# 2. Create the filled contour plot
filled.contour(xi_vals, theta_vals, plot_matrix_theta_xi,
               color.palette = function(n) hcl.colors(n, "Viridis"),
               xlab = expression(xi), ylab = expression(theta),
               main = "Log-Likelihood Gradient Surface",
               plot.axes = {
                 axis(1); axis(2)
                 contour(xi_vals, theta_vals, plot_matrix_theta_xi, add = TRUE, col = "white", alpha = 0.5)
                 points(true_margin["xi"], true_theta, col = "red", pch = 19, cex = 1.5)
                 points(res_xi, res_theta, col = "blue", pch = 19, cex = 1.5)
               })


# --- Parameters for Grid ---
kappa_vals <- seq(0.5, 15, length.out = 100) # Range around true_kappa = 6
xi_vals    <- seq(-0.1, 0.4, length.out = 100)

# Fix sigma and theta at True values to see the "clean" surface
sigma_fixed <- true_margin["sigma"]
theta_fixed <- true_theta

loglik_kappa_xi <- matrix(NA, nrow = length(kappa_vals), ncol = length(xi_vals))

for (i in seq_along(kappa_vals)) {
  for (j in seq_along(xi_vals)) {
    # theta_vec: [log_kappa, log_sigma, xi, log_theta_minus_1]
    t_vec <- c(log(kappa_vals[i]), log(sigma_fixed), xi_vals[j], log(theta_fixed - 1))
    loglik_kappa_xi[i, j] <- -egpd_gumbel_nll(t_vec, x)
  }
}

max_loglik_kappa_xi <- max(loglik_kappa_xi, na.rm = TRUE)
plot_matrix_kappa_xi <- loglik_kappa_xi
plot_matrix_kappa_xi[plot_matrix_kappa_xi < (max_loglik_kappa_xi - 500)] <- max_loglik_kappa_xi - 500 

# Visualization
filled.contour(kappa_vals, xi_vals, plot_matrix_kappa_xi,
               xlab = expression(kappa), ylab = expression(xi),
               main = "Log-Likelihood Surface: kappa vs xi",
               color.palette = function(n) hcl.colors(n, "Viridis"),
               plot.axes = {
                 axis(1); axis(2)
                 contour(kappa_vals, xi_vals, plot_matrix_kappa_xi, add = TRUE, col = "white")
                 points(true_margin["kappa"], true_margin["xi"], col = "red", pch = 19, cex = 1.5)
                 points(res_kappa, res_xi, col = "blue", pch = 19, cex = 1.5)
               })



profile_loglik_xi_theta <- function(x, xi_vals, theta_vals, init_par = c(log(5), log(1))) {
  
  n_xi <- length(xi_vals)
  n_theta <- length(theta_vals)
  
  prof_mat <- matrix(NA, nrow = n_xi, ncol = n_theta)
  
  for (i in seq_along(xi_vals)) {
    for (j in seq_along(theta_vals)) {
      
      xi_fixed <- xi_vals[i]
      theta_fixed <- theta_vals[j]
      
      # Objective: optimize over log(kappa), log(sigma)
      obj_fn <- function(par_sub) {
        t_vec <- c(
          par_sub[1],                      # log kappa
          par_sub[2],                      # log sigma
          xi_fixed,
          log(theta_fixed - 1)
        )
        
        egpd_gumbel_nll(t_vec, x)
      }
      
      opt <- try(
        optim(init_par, obj_fn, method = "BFGS"),
        silent = TRUE
      )
      
      if (!inherits(opt, "try-error")) {
        prof_mat[i, j] <- -opt$value
        init_par <- opt$par  # warm start (VERY important)
      }
    }
  }
  
  return(prof_mat)
}

xi_vals <- seq(-0.15, 0.6, length.out = 60)
theta_vals <- seq(1.05, 6, length.out = 60)

prof_surface <- profile_loglik_xi_theta(x, xi_vals, theta_vals)

max_val <- max(prof_surface, na.rm = TRUE)

plot_surface <- prof_surface
plot_surface[plot_surface < (max_val - 200)] <- max_val - 200

filled.contour(
  xi_vals, theta_vals, plot_surface,
  color.palette = function(n) hcl.colors(n, "Viridis"),
  xlab = expression(xi),
  ylab = expression(theta),
  main = "Profile Log-Likelihood: xi vs theta",
  plot.axes = {
    axis(1); axis(2)
    contour(xi_vals, theta_vals, plot_surface, add = TRUE, col = "white")
    
    points(true_margin["xi"], true_theta, col = "red", pch = 19, cex = 1.5)
    points(res_xi, res_theta, col = "blue", pch = 19, cex = 1.5)
  }
)


profile_loglik_kappa_xi <- function(x, kappa_vals, xi_vals, init_par = c(log(1), log(1))) {
  
  n_kappa <- length(kappa_vals)
  n_xi <- length(xi_vals)
  
  prof_mat <- matrix(NA, nrow = n_kappa, ncol = n_xi)
  
  for (i in seq_along(kappa_vals)) {
    for (j in seq_along(xi_vals)) {
      
      kappa_fixed <- kappa_vals[i]
      xi_fixed <- xi_vals[j]
      
      obj_fn <- function(par_sub) {
        # par_sub = (log_sigma, log_theta_minus_1)
        t_vec <- c(
          log(kappa_fixed),
          par_sub[1],
          xi_fixed,
          par_sub[2]
        )
        
        egpd_gumbel_nll(t_vec, x)
      }
      
      opt <- tryCatch(
        optim(init_par, obj_fn, method = "BFGS"),
        error = function(e) NULL
      )
      
      if (!is.null(opt)) {
        prof_mat[i, j] <- -opt$value
        init_par <- opt$par  # warm start
      }
    }
  }
  
  return(prof_mat)
}

kappa_vals <- seq(0.5, 15, length.out = 60)
xi_vals <- seq(-0.1, 0.5, length.out = 60)

prof_kappa_xi <- profile_loglik_kappa_xi(x, kappa_vals, xi_vals)

max_val <- max(prof_kappa_xi, na.rm = TRUE)

plot_surface <- prof_kappa_xi
plot_surface[plot_surface < (max_val - 200)] <- max_val - 200

filled.contour(
  kappa_vals, xi_vals, plot_surface,
  color.palette = function(n) hcl.colors(n, "Viridis"),
  xlab = expression(kappa),
  ylab = expression(xi),
  main = "Profile Log-Likelihood: kappa vs xi",
  plot.axes = {
    axis(1); axis(2)
    contour(kappa_vals, xi_vals, plot_surface, add = TRUE, col = "white")
    
    points(true_margin["kappa"], true_margin["xi"], col = "red", pch = 19, cex = 1.5)
    points(res_kappa, res_xi, col = "blue", pch = 19, cex = 1.5)
  }
)

profile_kappa <- sapply(kappa_vals, function(kappa_fixed) {
  
  obj_fn <- function(par) {
    # par = (log_sigma, xi, log_theta_minus_1)
    t_vec <- c(
      log(kappa_fixed),
      par[1],
      par[2],
      par[3]
    )
    
    egpd_gumbel_nll(t_vec, x)
  }
  
  opt <- optim(c(log(1), 0.1, log(1)), obj_fn)
  -opt$value
})

plot(kappa_vals, profile_kappa, type = "l",
     main = "Profile Likelihood for kappa",
     xlab = expression(kappa), ylab = "logLik")

abline(v = true_margin["kappa"], col = "red")