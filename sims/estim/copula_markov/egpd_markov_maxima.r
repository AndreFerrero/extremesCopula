############################################################
## Frequentist Copula–Markov Bootstrap for Maxima
## Gumbel copula + EGPD margins (using margin_egp object)
## WITH numerical-failure tracking
############################################################

library(copula)
library(pbapply)

margin_egp <- list(
  name = "egp",

  G_dist = function(u, param) u^param["kappa"],
  G_inv  = function(u, param) u^(1 / param["kappa"]),
  g_dist = function(u, param) {
    param["kappa"] * u^(param["kappa"] - 1)
  },

  cdf = function(x, param) {
    u <- evd::pgpd(x / param["sigma"], shape = param["xi"])
    margin_egp$G_dist(u, param)
  },

  lpdf = function(x, param) {
    u <- evd::pgpd(x / param["sigma"], shape = param["xi"])

    log(margin_egp$g_dist(u, param)) +
      log(evd::dgpd(x / param["sigma"], shape = param["xi"])) -
      log(param["sigma"])
  },

  quantile = function(p, param) {
    param["sigma"] * evd::qgpd(
      margin_egp$G_inv(p, param),
      shape = param["xi"]
    )
  },

  sample = function(n, param) {
    U <- runif(n)
    margin_egp$quantile(U, param)
  }
)
set.seed(123)

############################################################
## 1. SAFE INVERSION WITH FAILURE TRACKING
############################################################

make_safe_step <- function() {
  counter <- 0L
  
  step <- function(u_prev, cop, max_tries = 10, eps = 1e-10) {
    for (k in 1:max_tries) {
      v <- runif(1, eps, 1 - eps)
      
      u_next <- tryCatch(
        cCopula(
          cbind(u_prev, v),
          copula = cop,
          inverse = TRUE
        )[2],
        error = function(e) NA_real_
      )
      
      if (is.finite(u_next) && u_next > 0 && u_next < 1) {
        return(u_next)
      }
      
      counter <<- counter + 1L
    }
    
    stop("Copula inversion failed after repeated attempts")
  }
  
  list(
    step = step,
    get_failures = function() counter
  )
}

############################################################
## 2. SIMULATE GUMBEL COPULA MARKOV WITH EGPD MARGINS
############################################################

simulate_gumbel_markov_egpd <- function(n, alpha, param_egp, burn = 100) {
  cop <- gumbelCopula(alpha)
  tracker <- make_safe_step()
  
  # simulate uniform Markov chain
  U <- numeric(n + burn)
  U[1] <- runif(1, 1e-10, 1 - 1e-10)
  
  for (t in 1:(n + burn - 1)) {
    U[t + 1] <- tracker$step(U[t], cop)
  }
  
  list(
    X = margin_egp$quantile(U[(burn + 1):(n + burn)], param_egp),
    n_failures = tracker$get_failures()
  )
}

############################################################
## 3. FIT MARKOV COPULA MODEL (FREQUENTIST)
############################################################

fit_gumbel_markov <- function(X) {
  Uhat <- rank(X) / (length(X) + 1)
  
  fit <- fitCopula(
    gumbelCopula(dim = 2),
    cbind(Uhat[-length(Uhat)], Uhat[-1]),
    method = "ml"
  )
  
  coef(fit)
}

############################################################
## 4. SINGLE-EXPERIMENT DEMO
############################################################

n <- 500
alpha_true <- 2
param_egp <- c(sigma = 1, xi = 0.2, kappa = 0.5)

sim_obs <- simulate_gumbel_markov_egpd(n, alpha_true, param_egp)
X_obs <- sim_obs$X

par(mfrow = c(2,1))
plot(1:n, X_obs, type = "l")
hist(X_obs)
par(mfrow = c(1,1))

cat("Numerical failures during simulation:", sim_obs$n_failures, "\n")

alpha_hat <- fit_gumbel_markov(X_obs)
cat("Fitted alpha:", alpha_hat, "\n")

B_boot <- 100

maxima_boot <- pbapply::pbreplicate(
  B_boot,
  max(simulate_gumbel_markov_egpd(n, alpha_hat, param_egp)$X)
)

maxima_true <- pbapply::pbreplicate(
  B_boot,
  max(simulate_gumbel_markov_egpd(n, alpha_true, param_egp)$X)
)

plot(density(maxima_true), col = "black", lwd = 2,
     main = "Block Maxima with EGPD Margins",
     xlab = "Maximum")
lines(density(maxima_boot), col = "red", lwd = 2)
legend("topright", legend = c("True model", "Fitted model"),
       col = c("black", "red"), lwd = 2)

############################################################
## 5. OUTER MONTE CARLO STUDY
############################################################

R <- 50  # fewer reps for demo
B_boot <- 500

alpha_hat_vec <- numeric(R)
q95_vec <- numeric(R)
failures_vec <- numeric(R)

for (r in 1:R) {
  sim <- simulate_gumbel_markov_egpd(n, alpha_true, param_egp)
  X <- sim$X
  failures_vec[r] <- sim$n_failures
  
  alpha_hat <- fit_gumbel_markov(X)
  alpha_hat_vec[r] <- alpha_hat
  
  maxima_boot <- replicate(
    B_boot,
    max(simulate_gumbel_markov_egpd(n, alpha_hat, param_egp)$X)
  )
  
  q95_vec[r] <- quantile(maxima_boot, 0.95)
}

############################################################
## 6. DIAGNOSTICS
############################################################

par(mfrow = c(1, 3))

hist(alpha_hat_vec, breaks = 20, main = expression(hat(alpha)), xlab = expression(hat(alpha)), col = "grey")
hist(q95_vec, breaks = 20, main = "Bootstrap 95% quantile of M_n", xlab = "q_0.95", col = "grey")
hist(failures_vec, breaks = 20, main = "Numerical failures per path", xlab = "# failures", col = "grey")

par(mfrow = c(1, 1))

############################################################
## END
############################################################
