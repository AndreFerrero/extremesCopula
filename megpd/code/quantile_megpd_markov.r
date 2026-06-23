library(evd)
library(tidyverse)
library(copula)

###########################################################
# 1. ROBUST CORE FUNCTIONS
###########################################################

# source("megpd/code/markov_megpd_uniroot.r")
source("megpd/code/markov_megpd_u_uniroot.r")

set.seed(42)
simulate_R_quantiles <- function(
    M,
    n_steps,
    probs = seq(0.001, 0.999, length.out = 400),
    kappa,
    sigma,
    xi,
    delta_func
) {

    Q <- matrix(
        NA_real_,
        nrow = M,
        ncol = length(probs)
    )

    success <- logical(M)
    error_message <- rep(NA_character_, M)

    for (m in seq_len(M)) {

        cat(
            "Replication", m, "/", M, "\n"
        )

        res <- tryCatch(

            {

                sim <- simulate_megpd_chain(
                    n_steps = n_steps,
                    kappa = kappa,
                    sigma = sigma,
                    xi = xi,
                    delta_func = delta_func,
                    show_progress = FALSE
                )

                x <- sim$final_chain

                Rt <- x[-1] + x[-length(x)]

                Q[m, ] <- quantile(
                    Rt,
                    probs = probs,
                    names = FALSE,
                    type = 8
                )

                success[m] <- TRUE

                NULL

            },

            error = function(e) {

                error_message[m] <<- conditionMessage(e)

                cat(
                    "FAILED:",
                    conditionMessage(e),
                    "\n"
                )

                NULL

            }

        )
    }

    list(
        probs = probs,
        Q = Q,
        success = success,
        error_message = error_message
    )
}

kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

delta_strong_upper <- function(r) {
  # Starts at 0.8 (weak dependence at 0) and decays to 0.1 (strong dependence at infinity)
  0.2 + 0.6 * exp(-r / 5)
}

sim_qq <- simulate_R_quantiles(
    M = 100,
    n_steps = 100000,
    probs = seq(0.001, 0.999, length.out = 400),
    kappa = kappa_val,
    sigma = sigma_val,
    xi = xi_val,
    delta_func = delta_strong_upper
)

dim(sim_qq$Q)

# save(sim_qq, file = "megpd/res/radius_simqq_uvar_M100_n100000_xi05.Rdata")

Q <- sim_qq$Q[1:100,]
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

ggsave("megpd/figures/radius_quantiles_M100_n10000_xi05.png", p_qq)

load("C:/Users/Andrea Ferrero/extremesCopula/megpd/res/radius_simqq_M200_n10000_xi05.Rdata")
# p_qq

sim_qq$success
mean(sim_qq$success)
