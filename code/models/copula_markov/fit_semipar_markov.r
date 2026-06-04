############################################################
# SEMIPARAMETRIC EGPD + KDE / GUMBEL COPULA MARKOV MODEL
# FIT VIA EM ALGORITHM (GEM FOR KDE COMPONENT)
############################################################

library(egpd)
library(copula)
library(kdecopula)

############################################################
# 1. Utilities
############################################################

clamp01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

############################################################
# 2. EGPD transformation
############################################################

egpd_to_u <- function(x, sigma, xi, kappa) {
  u <- egpd:::pegpd(
    x,
    sigma = sigma,
    xi = xi,
    kappa = kappa,
    type = 1
  )
  clamp01(u)
}

############################################################
# 3. Build KDE copula
############################################################

build_kde <- function(U) {
  kdecop(U, method = "TLL2nn")
}

############################################################
# 4. E-step (responsibilities)
############################################################

e_step <- function(U, alpha, theta, kde_fit) {

  c_kde <- dkdecop(U, kde_fit)
  c_kde <- pmax(c_kde, 1e-12)

  gum <- gumbelCopula(theta)
  c_gum <- dCopula(U, gum)
  c_gum <- pmax(c_gum, 1e-12)

  tau <- alpha * c_kde /
    (alpha * c_kde + (1 - alpha) * c_gum)

  tau <- clamp01(tau)

  list(
    tau = tau,
    c_kde = c_kde,
    c_gum = c_gum
  )
}

############################################################
# 5. M-step updates
############################################################

m_step_alpha <- function(tau) {
  mean(tau)
}

m_step_theta <- function(U, w, theta_init) {

  obj <- function(log_theta) {

    theta <- 1 + exp(log_theta)

    dens <- dCopula(U, gumbelCopula(theta))
    dens <- pmax(dens, 1e-12)

    -(sum(w * log(dens)))
  }

  opt <- optim(
    par = log(theta_init - 1),
    fn = obj,
    method = "BFGS"
  )

  1 + exp(opt$par)
}

############################################################
# 6. EM algorithm
############################################################

fit_semiparametric_markov_em <- function(
  x,
  outer_iter = 10,
  verbose = TRUE
) {

  ##########################################################
  # Initial EGPD fit
  ##########################################################

  init <- egpd::fitegpd(x, type = 1, family = "egpd")

  sigma <- init$estimate[1]
  xi <- init$estimate[2]
  kappa <- init$estimate[3]

  ##########################################################
  # Transform to uniforms
  ##########################################################

  u <- egpd_to_u(x, sigma, xi, kappa)

  U <- cbind(u[-length(u)], u[-1])

  ##########################################################
  # Initial KDE
  ##########################################################

  kde_fit <- build_kde(U)

  ##########################################################
  # Initialize parameters
  ##########################################################

  alpha <- 0.9
  theta <- 2

  ##########################################################
  # EM LOOP
  ##########################################################

  for (k in 1:outer_iter) {

    if (verbose) {
      cat("\n=====================================\n")
      cat("EM ITERATION:", k, "\n")
      cat("=====================================\n")
    }

    ######################################################
    # E-step
    ######################################################

    e <- e_step(U, alpha, theta, kde_fit)

    tau <- e$tau

    ######################################################
    # M-step: mixture weight
    ######################################################

    alpha <- m_step_alpha(tau)

    ######################################################
    # M-step: Gumbel parameter
    ######################################################

    theta <- m_step_theta(U, 1 - tau, theta)

    ######################################################
    # M-step: KDE update (GEM approximation)
    ######################################################

    idx <- sample(
      1:nrow(U),
      size = nrow(U),
      replace = TRUE,
      prob = tau + 1e-8
    )

    kde_fit <- build_kde(U[idx, , drop = FALSE])

    ######################################################
    # log-likelihood
    ######################################################

    c_kde <- dkdecop(U, kde_fit)
    c_kde <- pmax(c_kde, 1e-12)

    c_gum <- dCopula(U, gumbelCopula(theta))
    c_gum <- pmax(c_gum, 1e-12)

    ll <- sum(log(alpha * c_kde + (1 - alpha) * c_gum))

    if (verbose) {
      cat("alpha  =", alpha, "\n")
      cat("theta  =", theta, "\n")
      cat("loglik =", ll, "\n")
    }
  }

  ##########################################################
  # Final output
  ##########################################################

  list(
    alpha = alpha,
    theta = theta,
    sigma = sigma,
    xi = xi,
    kappa = kappa,
    kde_fit = kde_fit,
    uniforms = u,
    loglik = ll
  )
}


############################################################
# Fit model
############################################################

fit <- fit_semiparametric_markov_em(
  mix_data$x,
  outer_iter = 5
)

############################################################
# Results
############################################################

print(fit$estimate)

cat("\n")

cat(
  "Estimated upper tail dependence:",
  fit$upper_tail_dependence,
  "\n"
)

############################################################
# 6. Diagnostics
############################################################

u <- fit$uniforms

plot(
  u[-length(u)],
  u[-1],
  pch = 16,
  cex = 0.4,
  main = "Lagged EGPD Uniforms"
)

############################################################
# Empirical upper tail dependence
############################################################

empirical_lambda <- function(u, q = 0.95) {

  n <- length(u)

  ind <- (
    u[-length(u)] > q &
      u[-1] > q
  )

  num <- sum(ind)

  den <- sum(u[-length(u)] > q)

  return(num / den)
}

cat(
  "\nEmpirical lambda_U (0.95): ",
  empirical_lambda(u, 0.95),
  "\n"
)