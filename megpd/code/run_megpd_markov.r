# source("megpd/markov_megpd_uniroot.r")
source("megpd/markov_megpd_u_uniroot.r")

set.seed(1)
n <- 10000
kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

# Dependence that gets STRONGER as values get LARGER
delta_strong_upper <- function(r) {
  0.2 + 0.6 * exp(-r / 5)
}

# Run
set.seed(4)

sim <- simulate_megpd_chain(
  n_steps = n,
  kappa = kappa_val,
  sigma = sigma_val,
  xi = xi_val,
  delta_func = delta_strong_upper,
  x0 = 1,
  burn_in_prop = 0.10,
  show_progress = TRUE
)

final_chain <- sim$final_chain
