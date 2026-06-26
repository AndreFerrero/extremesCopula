library(JuliaCall)
# julia_setup()
# cluster location
julia_setup("julia/1.9.4/julia-1.9.4/bin")

julia_source("dyn_mix_copula/markov_dyn_mix_copula.jl")

n_steps <- 1000

u_chain <- julia_call("simulate_chain", 
                       as.integer(n_steps), 1.5, 0.5, 2)

x_chain <- egpd::qegpd(u_chain, kappa = 2, sigma = 1, xi = 0.3)

# --- 7. VALIDATION ---
# Create consecutive pairs to check dependence
u_pairs <- cbind(u_chain[1:(n_steps-1)], u_chain[2:n_steps])

# par(mfrow=c(1,1))
# plot(u_chain, type="l", col="steelblue", main="Markov Chain Trace (u)", ylab="u_t")
# plot(u_pairs, pch=20, col=rgb(0,0,0,0.2), main="Transition Copula (u_t, u_{t+1})",
#      xlab="u_t", ylab="u_{t+1}")

# Check Tail Dependence (using mev package)
cat("\n--- Tail Dependence Analysis ---\n")
td <- mev::taildep(u_pairs, u = seq(0.7, 0.98, by = 0.01))


julia_source("dyn_mix_copula/fit_markov_dyn_mix_copula.jl")

fit <- julia_call("fit_markov_egpd", x_chain)

fit$convergence
fit$estimates
