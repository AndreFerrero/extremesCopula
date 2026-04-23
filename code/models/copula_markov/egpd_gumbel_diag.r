plot_density <- function(x, margin_est, title_cop = NULL, bins = NULL) {
  sigma <- margin_est[1]
  xi <- margin_est[2]
  kappa <- margin_est[3]


  hist(x, freq = FALSE, breaks = bins, main = title_cop)

  # Create a grid for smooth curve
  x_grid <- seq(min(x), max(x), length.out = 200)

  egpd_dens <- egpd::degpd_density(
    x_grid,
    kappa = kappa,
    sigma = sigma,
    xi = xi
  )

  # Overlay density
  lines(x_grid, egpd_dens, col = "red", lwd = 2)
}

#' Calculate Marginal and Conditional PIT values
#' @param data The observed vector x
#' @param theta_vec The estimated parameter vector [log_kappa, log_sigma, xi, log_theta_minus_1]
#' @param h_dist_fn The h-function (conditional CDF) from your copula object
get_pit_values <- function(x, theta_vec, h_dist_fn) {
  # 1. Parameter Extraction
  sigma <- theta_vec[1]
  xi <- theta_vec[2]
  kappa <- theta_vec[3]
  theta_c <- theta_vec[4]

  # 2. Marginal PIT (u_t)
  u <- egpd:::pegpd(x, sigma = sigma, xi = xi, kappa = kappa, type = 1)

  # 3. Conditional PIT (w_t)
  n <- length(u)
  w <- h_dist_fn(u[2:n], u[1:(n - 1)], theta_c)

  return(list(u = u, w = w))
}

#' Generate Side-by-Side QQ-plots
#' @param pit_list A list containing 'u' and 'w' from get_pit_values
plot_diag_plots <- function(pit_list) {
  u <- pit_list$u
  w <- pit_list$w

  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par)) # Restores plot settings after function finishes

  # Plot A: Marginal Check
  plot(stats::ppoints(length(u)), sort(u),
    main = "Marginal QQ-plot (EGPD)",
    xlab = "Theoretical Uniform", ylab = "Empirical u_t",
    pch = 20, col = "grey60"
  )
  abline(0, 1, col = "firebrick", lwd = 2)

  # Plot B: Conditional Check
  plot(stats::ppoints(length(w)), sort(w),
    main = "Conditional QQ-plot (Gumbel)",
    xlab = "Theoretical Uniform", ylab = "Empirical w_t",
    pch = 20, col = "royalblue"
  )
  abline(0, 1, col = "firebrick", lwd = 2)
}