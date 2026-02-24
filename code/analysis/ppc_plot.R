ppc_plot <- function(ppc_list, ci = c(0.025, 0.975), bins = 30, title = "Posterior Predictive") {
  # ppc_list: list(ppc_vals = numeric vector, obs_val = numeric)
  # ci: credible interval quantiles
  # bins: histogram bins
  # title: plot title

  ppc_vals <- ppc_list$ppc_stat
  obs_val  <- ppc_list$obs_stat

  median_val <- median(ppc_vals, na.rm = TRUE)
  ci_vals    <- quantile(ppc_vals, ci, na.rm = TRUE)

  df <- data.frame(ppc_vals = ppc_vals)

  ggplot(df, aes(x = ppc_vals)) +
    geom_histogram(bins = bins, fill = "skyblue", color = "black", alpha = 0.6) +
    # Credible interval shading
    annotate("rect", xmin = ci_vals[1], xmax = ci_vals[2],
             ymin = 0, ymax = Inf, alpha = 0.2, fill = "blue") +
    # Observed value
    geom_vline(xintercept = obs_val, color = "red", linetype = "dashed", size = 1) +
    # Posterior median
    geom_vline(xintercept = median_val, color = "blue", linetype = "solid", size = 1) +
    labs(title = title, x = "Value", y = "Frequency") +
    theme_minimal()
}
