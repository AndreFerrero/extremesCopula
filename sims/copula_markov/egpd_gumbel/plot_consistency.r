#' Plot Simulation Consistency and Distribution
#'
#' @param sim_results Data frame of results.
#' @param target_param The parameter to plot ("xi", "sigma", "kappa", "theta").
#' @param type Plot type: "median" for lines/points, "boxplot" for distributions.
#' @param scale Scale of the y-axis: "bias" (estimate - true) or "raw" (the estimate itself).
#' @param target_threshold Probability threshold for GPD models (e.g., 0.90).
#' @param models_to_show Vector of model names to include (NULL shows all).
#' @param free_y Logical, if TRUE facets have independent y-scales.
#' @param auto_zoom Logical, if TRUE focuses the y-axis on the 5th-95th percentiles (useful for boxplots).
#' Unified Study Plotting Function
#'
#' @param sim_results Data frame of simulation results
#' @param target_param Parameter to visualize ("xi", "sigma", "kappa", "theta")
#' @param type "median" for line plots, "boxplot" for distributions
#' @param scale "bias" for (est - true) or "raw" for absolute values
#' @param target_threshold The threshold probability to display for GPD models
#' @param models_to_show Character vector of models to include (e.g., c("IID_EGPD", "EGPD_GUMBEL"))
#' @param free_y Logical. If TRUE, facets have independent y-axis scales
#' @param auto_zoom Logical. If TRUE, focuses on the 5th-95th percentile bulk

plot_consistency <- function(sim_results,
                             target_param = "xi",
                             type = c("median", "boxplot"),
                             scale = c("bias", "raw"),
                             target_threshold = 0.90,
                             models_to_show = NULL,
                             free_y = FALSE,
                             auto_zoom = TRUE,
                             title = NULL,
                             subtitle = NULL,
                             caption = NULL,
                             xlab = NULL,
                             ylab = NULL,
                             base_size = 12,
                             legend_position = "top") {
  type <- match.arg(type)
  scale <- match.arg(scale)

  # ============================================================================
  # GLOBAL STYLE DEFINITION
  # ============================================================================
  model_colors <- c(
    "IID_EGPD"            = "#00BFC4", # Cyan
    "EGPD_GUMBEL"         = "#F8766D", # Coral
    "EGPD_GUMBEL_IFM"     = "yellow",
    "EGPD_JOE"            = "#7CAE00", # Olive Green
    "IID_GPD"             = "#C77CFF", # Purple
    "GPD_Declustering"    = "#00BA38", # Green
    "Censored_GPD_GUMBEL" = "#DE8C00", # Orange
    "Hill"                = "#619CFF", # Light Blue
    "Hill_BC"             = "#FF61CC" # Pink
  )

  model_ltypes <- c(
    "IID_EGPD"            = "dashed",
    "EGPD_GUMBEL"         = "solid",
    "EGPD_GUMBEL_IFM"     = "solid",
    "EGPD_JOE"            = "solid",
    "IID_GPD"             = "dotted",
    "GPD_Declustering"    = "dotdash",
    "Censored_GPD_GUMBEL" = "longdash",
    "Hill"                = "solid",
    "Hill_BC"             = "solid"
  )

  model_order_fixed <- c(
  "EGPD_GUMBEL",
  "EGPD_GUMBEL_IFM",
  "IID_EGPD",
  "IID_GPD",
  "GPD_Declustering",
  "Censored_GPD_GUMBEL",
  "Hill",
  "Hill_BC"
  )
  # ============================================================================
  # FILTER DATA
  # ============================================================================
  plot_data <- sim_results %>%
    filter(parameter == target_param) %>%
    filter(is.na(threshold) | threshold == target_threshold)

  if (!is.null(models_to_show)) {
    plot_data <- plot_data %>% filter(model %in% models_to_show)
  }

  if (nrow(plot_data) == 0) {
    stop("No data found for the specified filters.")
  }

  all_models <- unique(plot_data$model)

  model_levels <- c(
    model_order_fixed,
    setdiff(all_models, model_order_fixed)
  )

  plot_data$model <- factor(plot_data$model, levels = model_levels)

  # ============================================================================
  # DEFAULT LABELS
  # ============================================================================

  # Auto subtitle only if user does not provide one
  if (is.null(subtitle)) {
    is_threshold_dependent <- any(!is.na(plot_data$threshold))

    subtitle <- if (is_threshold_dependent) {
      sprintf("Threshold: %.2f", target_threshold)
    } else {
      NULL
    }
  }

  # Default title
  if (is.null(title)) {
    title <- sprintf(
      "%s Analysis: %s",
      tools::toTitleCase(type),
      target_param
    )
  }

  # Default x label
  if (is.null(xlab)) {
    xlab <- "Sample Size (n)"
  }

  # Default y label
  if (is.null(ylab)) {
    ylab <- if (scale == "bias") {
      "Median Bias"
    } else {
      "Parameter Estimate"
    }
  }

  # ============================================================================
  # Y VARIABLE
  # ============================================================================
  y_var <- if (scale == "bias") "bias" else "estimate"

  # ============================================================================
  # REFERENCE LINES
  # ============================================================================
  ref_lines <- plot_data %>%
    group_by(theta_true) %>%
    summarise(
      intercept = if (scale == "bias") 0 else unique(true_value)[1],
      .groups = "drop"
    )

  # ============================================================================
  # DATA TRANSFORMATION
  # ============================================================================
  if (type == "median") {
    plot_data <- plot_data %>%
      group_by(model, theta_true, n) %>%
      summarise(
        y_val = median(.data[[y_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(n_factor = factor(n))
  } else {
    plot_data <- plot_data %>%
      rename(y_val = !!sym(y_var)) %>%
      mutate(n_factor = factor(n))
  }

  # ============================================================================
  # BASE PLOT
  # ============================================================================
  p <- ggplot(
    plot_data,
    aes(x = n_factor, y = y_val)
  ) +
    geom_hline(
      data = ref_lines,
      aes(yintercept = intercept),
      color = "grey40",
      linetype = "dashed",
      linewidth = 0.6
    )

  # ============================================================================
  # GEOMS
  # ============================================================================
  if (type == "median") {
    p <- p +
      geom_line(
        aes(
          group = model,
          color = model,
          linetype = model
        ),
        linewidth = 1
      ) +

      geom_point(aes(color = model),
        size = 2
      ) +

      scale_color_manual(values = model_colors) +

      scale_linetype_manual(values = model_ltypes)
  } else {
    p <- p +
      geom_boxplot(aes(fill = model),
        outlier.size = 0.5,
        alpha = 0.7,
        lwd = 0.3
      ) +

      scale_fill_manual(values = model_colors)
  }

  # ============================================================================
  # FACETING + THEME
  # ============================================================================
  p <- p +
    facet_wrap(
      ~theta_true,
      labeller = as_labeller(
        function(x) {
          sapply(x, function(val) {
            bquote(theta == .(val))
          })
        },
        default = label_parsed
      ),
      scales = if (free_y) "free_y" else "fixed"
    ) +

    theme_bw(base_size = base_size) +

    theme(
      legend.position = legend_position,
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )  +
    guides(
      color = guide_legend(nrow = 1),
      fill = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
    ) +
    labs(
      title = if (nzchar(title)) title else NULL,
      subtitle = if (nzchar(subtitle)) subtitle else NULL,
      caption = caption,
      x = xlab,
      y = ylab
    )

  # ============================================================================
  # SMART ZOOMING
  # ============================================================================
  if (auto_zoom && type == "boxplot") {
    y_lims <- quantile(
      plot_data$y_val,
      probs = c(0.01, 0.99),
      na.rm = TRUE
    )

    p <- p +
      coord_cartesian(
        ylim = c(
          y_lims[1] - 0.15 * abs(y_lims[1]),
          y_lims[2] + 0.15 * abs(y_lims[2])
        )
      )
  }

  return(p)
}

# ==============================================================================
# EXECUTION
# ==============================================================================
library(dplyr)
library(ggplot2)
library(here)

load("sims/copula_markov/egpd_gumbel/res/consistency_mc300_u9095_neldermead_gumbeldata_ifm_joe_hill_kn05_kappa2_sigma1_xi01.RData")

width = 1200
height = 650

png("slides/figures/sim_consistency1.png",
  width = width, height = height, res = 100
)

plot_consistency(
  results,
  type = "boxplot",
  target_threshold = 0.95,
  scale = "raw",
  models_to_show = c("IID_EGPD", "EGPD_GUMBEL", "IID_GPD"),
  title = "",
  subtitle = "",
  xlab = NULL,
  ylab = expression(xi)
)
dev.off()

png("slides/figures/sim_consistency2.png",
  width = width, height = height, res = 100
)

plot_consistency(
  results,
  type = "boxplot",
  target_threshold = 0.95,
  scale = "raw",
  models_to_show = c("IID_EGPD", "EGPD_GUMBEL", "IID_GPD", "GPD_Declustering", "Censored_GPD_GUMBEL", "Hill", "Hill_BC"),
  title = "",
  subtitle = "",
  xlab = NULL,
  ylab = expression(xi)
)
dev.off()

png("slides/figures/sim_consistency3.png",
  width = width, height = height, res = 100
)

plot_consistency(results,
  type = "boxplot",
  scale = "raw",
  target_param = "theta", target_threshold = 0.95,
  models_to_show = c("Censored_GPD_GUMBEL", "EGPD_GUMBEL", "EGPD_GUMBEL_IFM"),
  title = "",
  subtitle = "",
  xlab = NULL,
  ylab = expression(theta),
)
dev.off()

png("slides/figures/sim_consistency4.png",
  width = width, height = height, res = 100
)

plot_consistency(
  results,
  type = "boxplot",
  target_threshold = 0.95,
  scale = "raw",
  models_to_show = c("EGPD_GUMBEL", "EGPD_JOE"),
  title = "",
  subtitle = "",
  xlab = NULL,
  ylab = expression(xi)
)
dev.off()

# plot_consistency(results,
#   type = "boxplot",
#   scale = "raw",
#   target_param = "sigma", target_threshold = 0.95,
#   models_to_show = c("IID_EGPD", "EGPD_GUMBEL")
# )

# plot_consistency(results,
#   type = "boxplot",
#   scale = "raw",
#   target_param = "kappa", target_threshold = 0.95,
#   models_to_show = c("IID_EGPD", "EGPD_GUMBEL")
# )

