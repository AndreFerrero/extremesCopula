source("code/models/copula_markov/simulate.r")
source("code/models/copula_markov/egpd_gumbel_diag.r")

test_fit <- fit_egpd_gumbel(egpd_gumbel_data$x)

test_fit$estimate

test_pit <- get_pit_values(egpd_gumbel_data$x, test_fit$estimate, copula_gumbel$h_dist)

plot_diag_plots(test_pit)
