library(JuliaCall)
library(tidyverse)

julia_setup()

julia_source("megpd/code/markov_megpd.jl")

n <- 100000
kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

x <- julia_call("simulate_megpd_chain", as.integer(n), kappa_val, sigma_val, xi_val)

final_chain <- x$final_chain

quantile(final_chain)
quantile(final_chain, 0.95)

julia_source("megpd/code/markov_megpd_cluster.jl")

clust_sim <- julia_call("clust_sim", 1.9, as.integer(100000), 2.0, 1.0, 0.2)

clust_sim$pi_dist
clust_sim$piC_dist


julia_source("megpd/code/quantile_megpd_markov.jl")

m <- 20
quants <- as.vector(seq(0.001, 0.999, length.out = 400))

sim_qq <- julia_call("simulate_megpd_quantiles",
                        as.integer(m),
                        as.integer(n),
                        quants,
                        kappa_val,
                        sigma_val,
                        xi_val
)

Q <- sim_qq$Q
p <- sim_qq$probs

theoretical <- egpd::qegpd(
    p,
    kappa = kappa_val,
    sigma = sigma_val,
    xi = xi_val
)

df <- tibble(
    p = p,
    theo = theoretical,
    lower = apply(Q, 2, quantile, 0.025, na.rm = TRUE),
    median = apply(Q, 2, quantile, 0.50, na.rm = TRUE),
    upper = apply(Q, 2, quantile, 0.975, na.rm = TRUE)
)

p_qq <- ggplot(df) +

    geom_ribbon(
        aes(
            x = theo,
            ymin = lower,
            ymax = upper
        ),
        alpha = 0.25
    ) +

    geom_line(
        aes(theo, median),
        linewidth = 1
    ) +

    geom_abline(
        slope = 1,
        intercept = 0,
        colour = "red",
        linetype = 2
    ) +

    scale_x_log10() +
    scale_y_log10() +

    labs(
        title = "Radius Quantile distribution",
        x = "Theoretical EGPD quantile",
        y = "Empirical quantile"
    ) +

    theme_bw()

ggsave("megpd/figures/julia_radius_quantiles_M20_n100000_xi05.png", p_qq)
