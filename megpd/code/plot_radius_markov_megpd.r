# source("megpd/markov_megpd_u_uniroot.r")
library(ggplot2)
library(patchwork)

###########################################################
#  RADIAL CHECK  (X_t + X_{t+1} ~ EGPD?)
###########################################################

#----------------------------------------------------------
# Construct radial series
#----------------------------------------------------------

Rt <- final_chain[-1] + final_chain[-length(final_chain)]

cat("\nRadial diagnostics\n")
print(summary(Rt))

#----------------------------------------------------------
# 1. Trace plot
#----------------------------------------------------------

df_trace <- data.frame(
  t = seq_along(Rt),
  R = Rt
)

p_trace <- ggplot(df_trace, aes(t, R)) +
  geom_line(
    colour = "steelblue",
    linewidth = 0.3,
    alpha = 0.8
  ) +
  labs(
    title = expression(R[t] == X[t] + X[t-1]),
    x = "Iteration",
    y = expression(R[t])
  ) +
  theme_bw(base_size = 13)

#----------------------------------------------------------
# 2. Density check
#----------------------------------------------------------

p_density <- ggplot(data.frame(R = Rt), aes(R)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 80,
    fill = "steelblue",
    colour = "white"
  ) +
  stat_function(
    fun = function(x)
      egpd::degpd_density(
        x,
        kappa = kappa_val,
        sigma = sigma_val,
        xi    = xi_val
      ),
    colour = "tomato",
    linewidth = 1.2
  ) +
  labs(
    title = "Radial density",
    x = expression(R[t]),
    y = "Density"
  ) +
  theme_bw(base_size = 13)

#----------------------------------------------------------
# 3. CDF check
#----------------------------------------------------------

ecdf_R <- ecdf(Rt)

R_grid <- seq(
  min(Rt),
  quantile(Rt, 0.995),
  length.out = 500
)

true_cdf <- egpd::pegpd(
  R_grid,
  kappa = kappa_val,
  sigma = sigma_val,
  xi    = xi_val
)

df_cdf <- data.frame(
  R = R_grid,
  empirical = ecdf_R(R_grid),
  theoretical = true_cdf
)

p_cdf <- ggplot(df_cdf) +
  geom_line(
    aes(R, theoretical,
        colour = "Theoretical EGPD"),
    linewidth = 1
  ) +
  geom_line(
    aes(R, empirical,
        colour = "Empirical"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_colour_manual(
    values = c(
      "Theoretical EGPD" = "tomato",
      "Empirical" = "steelblue"
    )
  ) +
  labs(
    title = "CDF comparison",
    x = expression(R[t]),
    y = "CDF",
    colour = NULL
  ) +
  theme_bw(base_size = 13)

#----------------------------------------------------------
# 4. Tail behaviour
#----------------------------------------------------------

tail_threshold <- egpd::qegpd(
  0.90,
  kappa = kappa_val,
  sigma = sigma_val,
  xi    = xi_val
)

R_tail <- sort(Rt[Rt > tail_threshold])

n_tail <- length(R_tail)

emp_surv <- (n_tail:1) / length(Rt)

true_surv <- 1 - egpd::pegpd(
  R_tail,
  kappa = kappa_val,
  sigma = sigma_val,
  xi    = xi_val
)

df_tail <- data.frame(
  R = R_tail,
  empirical = emp_surv,
  theoretical = true_surv
)

p_tail <- ggplot(df_tail) +
  geom_line(
    aes(R, theoretical,
        colour = "Theoretical EGPD"),
    linewidth = 1
  ) +
  geom_point(
    aes(R, empirical,
        colour = "Empirical"),
    alpha = 0.6,
    size = 0.8
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(
    values = c(
      "Theoretical EGPD" = "tomato",
      "Empirical" = "steelblue"
    )
  ) +
  labs(
    title = "Tail survival (log-log)",
    x = expression(R[t]),
    y = expression(P(R[t] > r)),
    colour = NULL
  ) +
  theme_bw(base_size = 13)

#----------------------------------------------------------
# 5. EGPD QQ plot
#----------------------------------------------------------

nR <- length(Rt)

u <- ppoints(nR)

q_theoretical <- egpd::qegpd(
  u,
  kappa = kappa_val,
  sigma = sigma_val,
  xi    = xi_val
)

df_qq <- data.frame(
  theoretical = sort(q_theoretical),
  empirical   = sort(Rt)
)

p_qq <- ggplot(
  df_qq,
  aes(theoretical, empirical)
) +
  geom_point(
    colour = "steelblue",
    alpha = 0.25,
    size = 1
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "tomato",
    linetype = "dashed"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "EGPD QQ plot",
    x = "Theoretical quantiles",
    y = "Empirical quantiles"
  ) +
  theme_bw(base_size = 13)

#----------------------------------------------------------
# 6. Radial dependence
#----------------------------------------------------------

df_dep <- data.frame(
  Rt  = Rt[-length(Rt)],
  Rt1 = Rt[-1]
)

p_dep <- ggplot(
  df_dep,
  aes(Rt, Rt1)
) +
  geom_point(
    alpha = 0.15,
    size = 0.6,
    colour = "steelblue"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw(base_size = 13)

#----------------------------------------------------------
# Hill tail-index estimate
#----------------------------------------------------------

# hill <- function(x, k = 100) {
#   xs <- sort(x, decreasing = TRUE)
#   mean(log(xs[1:k])) - log(xs[k + 1])
# }

# cat(
#   "\nHill xi estimate (k=50): ",
#   round(hill(Rt, 50), 3),
#   "\n"
# )

# source("code/models/hill_estimator.r")

# # hill_bc_hat(Rt, 50)

# cat(
#   "Hill xi estimate (k=100): ",
#   round(hill(Rt, 100), 3),
#   "\n"
# )

# cat(
#   "True xi:",
#   xi_val,
#   "\n"
# )

#----------------------------------------------------------
# Final panel
#----------------------------------------------------------

radial_panel <-
  (p_trace | p_density) /
  (p_cdf   | p_qq) /
  (p_tail  | p_dep)

plot(radial_panel)
# ggsave("megpd/radial_numfun_n10k_xi05.png", radial_panel, dpi = 600, height = 10, width = 15)
