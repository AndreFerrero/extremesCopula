# Load required libraries
library(mgcv)
library(egpd) # Ensure your specific EGPD package is loaded
source("code/fixed_bernstein_egpd.r")

fit_megpd_bivariate <- function(x, y, m_degree = NULL, k = 10) {
  # 0. Data Preparation
  # The model is defined for positive intensities (e.g., river discharge)
  n <- length(x)
  R <- x + y # Radial component (sum-norm)
  V <- log(x / y) # Angular component (log-ratio)
  V <- V - mean(V) # Centering V to have mean 0 (as per the Z ~ zero-mean assumption)

  cat("Step 1: Fitting Radial EGPD via Bernstein Polynomials...\n")

  # 1. First Step: Radial component EGPD
  # The paper suggests m = floor(0.5 * n / log(n)) for Bernstein degree
  if (is.null(m_degree)) {
    m_degree <- floor(0.5 * n / log(n))
  }

  fit_r <- egpd::fitegpd(
    R,
    type = 1,
    method = "bernstein",
    bernstein.m = m_degree
  )

  cat("Step 2: Fitting Heteroscedastic Angular Component via mgcv...\n")

  # 2. Second Step: Logistic-heteroscedastic modeling
  # We model V | R=r ~ N(0, delta(r)^2)
  # According to Eq (19), log(delta(r)) is a linear combination of basis functions.
  # We use mgcv's 'gauslss' family to model the scale (sigma) specifically.

  # In 'gauslss', the first formula is for the mean (mu), the second for log(sigma).
  # The paper assumes a zero-mean exchangeable vector Z (Eq 11).
  # Therefore, we fix intercept to 0 and no covariates for the mean.

  # K = 10 is the default basis dimension suggested in the paper (Section 4.1).
  fit_delta <- gam(
    list(
      V ~ 1, # Mean fixed at 0 (as per the Z ~ zero-mean assumption)
      ~ s(R, k = k, bs = "cr") # log(sigma) modeled as a cubic spline of R
    ),
    family = gaulss(),
    data = data.frame(V = V, R = R)
  )

  # 3. Return Results
  results <- list(
    radial_fit = fit_r,
    angular_fit = fit_delta
  )

  return(results)
}

plot_delta <- function(fit_delta, R) {
  R_grid <- seq(min(R), max(R), length.out = 500)

  newdat <- data.frame(R = R_grid, V = 0)

  pred <- predict(fit_delta, newdata = newdat, type = "link", se.fit = TRUE)

  eta <- pred$fit[, 2]
  se <- pred$se.fit[, 2]
  eta_lo <- eta - 1.96 * se
  eta_hi <- eta + 1.96 * se

  sigma_hat <- exp(eta) + 0.01
  sigma_lo <- exp(eta_lo) + 0.01
  sigma_hi <- exp(eta_hi) + 0.01

  plot(R_grid, sigma_hat,
    type = "l", lwd = 2,
    xlab = "R",
    ylab = expression(delta(R)),
    ylim = c(0, max(sigma_hi) * 1.2)
  )

  lines(R_grid, sigma_lo, col = "red", lty = 2)
  lines(R_grid, sigma_hi, col = "red", lty = 2)
}

fitted_chi <- function(fit_mc, p, fit_radial, fit_delta) {
  # fit_mc is the number of Monte Carlo samples to simulate from the fitted model

  kappa_fit <- fit_radial$estimate["kappa"]
  sigma_fit <- fit_radial$estimate["sigma"]
  xi_fit <- fit_radial$estimate["xi"]

  # simulate R from fitted EGPD
  Rsim <- egpd::regpd(
    fit_mc,
    kappa = kappa_fit,
    sigma = sigma_fit,
    xi = xi_fit
  )

  # fitted delta(R)
  eta <- predict(fit_delta,
    newdata = data.frame(R = Rsim),
    type = "link"
  )[, 2]

  delta <- exp(eta) + 0.01

  Vsim <- rnorm(fit_mc, sd = delta)

  Xsim <- Rsim * exp(Vsim) / (1 + exp(Vsim))
  Ysim <- Rsim / (1 + exp(Vsim))

  ux <- quantile(Xsim, p)
  uy <- quantile(Ysim, p)

  chi_hat <- mean(Xsim > ux & Ysim > uy) /
    mean(Xsim > ux)

  return(chi_hat)
}

plot_fitted_chi <- function(fit_chi = NULL, fit_megpd,
                            n_mc = 100, p = seq(0.8, 0.99, length.out = 20), fit_mc = 10000) {
  if (is.null(fit_chi)) {
    set.seed(123)

    fit_radial <- fit_megpd$radial_fit
    fit_delta <- fit_megpd$angular_fit

    fit_chi <- matrix(NA, nrow = n_mc, ncol = length(p))

    for (l in seq_len(n_mc)) {
      for (i in seq_along(p)) {
        fit_chi[l, i] <- fitted_chi(fit_mc, p[i], fit_radial, fit_delta)
      }
    }
  }

  ci <- apply(fit_chi, 2, quantile, probs = c(0.025, 0.975))
  plot(p, colMeans(fit_chi), type = "l", main = "Fitted Chi vs Threshold", xlab = "Threshold (p)", ylab = "Fitted Chi", ylim = c(0, 1))
  lines(p, ci[1, ], col = "grey20", lty = 2)
  lines(p, ci[2, ], col = "grey20", lty = 2)

  return(fit_chi)
}
