library(ggplot2)
library(dplyr)
library(tidyr)

load("sims/copula_markov/egpd_gumbel/res/hill_sim_n2000_mc300_kappa2_xi01.Rdata")

theta_grid <- c(1, 2, 4, 6)

plot_df <- bind_rows(lapply(theta_grid, function(theta) {

  res <- results[[paste0("theta_", theta)]]

  data.frame(
    theta = theta,
    k = res$k,

    Hill_bias = res$hill_bias,
    Hill_BC_bias = res$hill_bc_bias,

    Hill_rmse = res$hill_rmse,
    Hill_BC_rmse = res$hill_bc_rmse,

    Hill_rel_bias = res$hill_rel_bias,
    Hill_BC_rel_bias = res$hill_bc_rel_bias
  )
}))

bias_df <- plot_df %>%
  select(theta, k, Hill_bias, Hill_BC_bias) %>%
  pivot_longer(
    cols = c(Hill_bias, Hill_BC_bias),
    names_to = "estimator",
    values_to = "value"
  ) %>%
  mutate(estimator = sub("_bias", "", estimator))

rmse_df <- plot_df %>%
  select(theta, k, Hill_rmse, Hill_BC_rmse) %>%
  pivot_longer(
    cols = c(Hill_rmse, Hill_BC_rmse),
    names_to = "estimator",
    values_to = "value"
  ) %>%
  mutate(estimator = sub("_rmse", "", estimator))

rel_bias_df <- plot_df %>%
  select(theta, k, Hill_rel_bias, Hill_BC_rel_bias) %>%
  pivot_longer(
    cols = c(Hill_rel_bias, Hill_BC_rel_bias),
    names_to = "estimator",
    values_to = "value"
  ) %>%
  mutate(estimator = sub("_rel_bias", "", estimator))

# dev.new()
# ggplot(bias_df, aes(x = k, y = value, color = estimator)) +
#   geom_line(linewidth = 1) +
#   facet_wrap(
#       ~theta,
#       labeller = as_labeller(
#         function(x) {
#           sapply(x, function(val) {
#             bquote(theta == .(val))
#           })
#         },
#         default = label_parsed
#       ), scales = "free_y") +
#   labs(
#     x = "k",
#     y = "Bias",
#     color = "Estimator"
#   ) +
#   theme_bw() +
#   theme(legend.position = "top")

png("slides/figures/hill_rmse.png",
  width = 1200, height = 650, res = 100)

ggplot(rmse_df, aes(x = k, y = value, color = estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(
      ~theta,
      labeller = as_labeller(
        function(x) {
          sapply(x, function(val) {
            bquote(theta == .(val))
          })
        },
        default = label_parsed
      ), scales = "free_y") +
  labs(
    x = "k",
    y = "RMSE",
    color = "Estimator"
  ) +
  theme_bw() +
  theme(legend.position = "top")
dev.off()

png("slides/figures/hill_relbias.png",
  width = 1200, height = 650, res = 100)

ggplot(rel_bias_df, aes(x = k, y = value, color = estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(
      ~theta,
      labeller = as_labeller(
        function(x) {
          sapply(x, function(val) {
            bquote(theta == .(val))
          })
        },
        default = label_parsed
      ), scales = "free_y") +
  labs(
    x = "k",
    y = "Absolute Relative Bias",
    color = "Estimator"
  ) +
  theme_bw() +
  theme(legend.position = "top")

dev.off()