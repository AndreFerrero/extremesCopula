library(JuliaCall)

julia_setup()

julia_source("megpd/code/markov_megpd.jl")

n <- 1000000
kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

x <- julia_call("simulate_megpd_chain", as.integer(n), kappa_val, sigma_val, xi_val)

final_chain <- x$final_chain
