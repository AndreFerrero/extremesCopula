# ==============================================================================
# DIAGNOSTIC FUNCTIONS
# ==============================================================================

source("R/asymptotic_model.R")
source("R/simulation.R")

#' Compute QQ plot data
#'
#' @param M Vector of block maxima
#' @param fit optim object
#' @return Data frame with theoretical and empirical quantiles
qq_data <- function(M, fit) {

  n_obs <- length(M)
  prob_grid <- (1:n_obs) / (n_obs + 1)

  q_theory <- model_quantile(prob_grid, fit$mu, fit$sigma, fit$xi, fit$theta, fit$family)
  q_emp <- sort(M)

  data.frame(theoretical = q_theory, empirical = q_emp)
}

#' Plot QQ diagnostic
#'
#' @param M Vector of block maxima
#' @param fit Optim object
#' @param main Plot title
plot_qq <- function(M, fit,
                    main = "QQ Plot") {

  dat <- qq_data(M, fit$mu, fit$sigma, fit$xi, fit$theta, family)

  plot(dat$theoretical, dat$empirical,
       pch = 16, cex = 0.5,
       main = main,
       xlab = "Theoretical quantiles",
       ylab = "Empirical quantiles")
  abline(0, 1, col = "red")
}

#' Bootstrap QQ confidence bands
#'
#' @param M Observed block maxima
#' @param fit Optim object
#' @param n Dimension
#' @param k Number of replications
#' @param B Number of bootstrap replications
#' @param margin_fn Marginal transformation
#' @param conf_level Confidence level (default 0.95)
plot_qq_bootstrap <- function(M, fit,
                              n, k, B = 200,
                              conf_level = 0.95) {

  n_obs <- length(M)
  prob_grid <- (1:n_obs) / (n_obs + 1)

  # theoretical quantiles from fitted asymptotic model
  q_theory <- model_quantile(prob_grid,
                             fit$mu, fit$sigma, fit$xi,
                             fit$theta, fit$family)

  q_emp <- sort(M)

  ##########################################################
  # BOOTSTRAP FROM ASYMPTOTIC MODEL
  ##########################################################

  q_boot <- matrix(NA, nrow = B, ncol = n_obs)

  for (b in 1:B) {

    # simulate directly from fitted max distribution
    u_sim <- runif(n_obs)

    m_sim <- model_quantile(
      u_sim,
      fit$mu, fit$sigma, fit$xi,
      fit$theta, fit$family
    )

    q_boot[b, ] <- sort(m_sim)
  }

  alpha <- 1 - conf_level

  q_low  <- apply(q_boot, 2, quantile, alpha / 2)
  q_high <- apply(q_boot, 2, quantile, 1 - alpha / 2)
  q_med  <- apply(q_boot, 2, quantile, 0.5)

  ##########################################################
  # PLOT
  ##########################################################

  plot(q_theory, q_emp,
       pch = 16, cex = 0.5,
       main = sprintf("QQ plot with %d%% asymptotic bands",
                      round(conf_level * 100)),
       xlab = "Theoretical quantiles",
       ylab = "Empirical maxima")

  lines(q_theory, q_low, col = "gray", lwd = 2)
  lines(q_theory, q_high, col = "gray", lwd = 2)
  lines(q_theory, q_med, col = "blue", lwd = 2)

  abline(0, 1, col = "red")
}
