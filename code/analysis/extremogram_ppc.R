# Extremogram Calculation (Lag-dependent Tail dependence)
calc_extremogram <- function(x, prob, max_lag) {
    u <- quantile(x, prob)
    n <- length(x)
    ext_vec <- numeric(max_lag)
    is_ext <- x > u
    denom <- sum(is_ext)
    for (h in 1:max_lag) {
        num <- sum(is_ext[1:(n - h)] & is_ext[(h + 1):n])
        ext_vec[h] <- num / denom
    }
    return(ext_vec)
}

ppc_extremogram <- function(
  x, x_rep,
  n_rep = 100, prob = 0.95, max_lag = 10
) {
    x_ext <- calc_extremogram(x, prob = prob, max_lag = max_lag)
    x_rep_ext <- t(apply(x_rep[1:n_rep, ], 1, calc_extremogram, prob, max_lag))

    plot_df <- data.frame(
        lag = 1:10, obs = x_ext,
        low = apply(x_rep_ext, 2, quantile, 0.025),
        high = apply(x_rep_ext, 2, quantile, 0.975)
    )

    ggplot(plot_df, aes(x = lag)) +
        geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
        geom_line(aes(y = obs), color = "red", size = 1) +
        labs(title = "Extremogram PPC", subtitle = "95% Model Credible Interval") +
        geom_point(aes(y = obs), color = "red") +
        # Force integer breaks on the x-axis
        scale_x_continuous(breaks = 1:10) +
        theme_minimal()
}