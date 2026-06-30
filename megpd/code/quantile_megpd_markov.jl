using Statistics
using StatsBase
using Base.Threads

include("markov_megpd.jl")
const delta = Main.delta

"""
Run M replications of the MEGPD chain simulation and compute quantiles of Rt = x[t] + x[t-1].
"""
function simulate_megpd_quantiles(
    M::Int,
    n_steps::Int,
    probs::Vector{Float64},
    kappa::Float64,
    sigma::Float64,
    xi::Float64,
    x0::Float64 = 1.0,
    burn_in_prop::Float64 = 0.1,
)

    P = length(probs)

    Q = Matrix{Float64}(undef, M, P)
    success = falses(M)
    error_message = Vector{Union{String, Nothing}}(undef, M)

    Threads.@threads for m in 1:M
        try

            sim = simulate_megpd_chain(
                n_steps,
                kappa,
                sigma,
                xi,
                x0,
                burn_in_prop,
                false
            )

            x = sim.final_chain

            # Rt = x_t + x_{t-1}
            Rt = x[2:end] .+ x[1:end-1]

            # quantiles (StatsBase)
            Q[m, :] .= quantile(Rt, probs)

            success[m] = true
            error_message[m] = nothing

        catch e
            success[m] = false
            error_message[m] = sprint(showerror, e)
            @info "FAILED replication $m" exception=e
        end
    end

    return (
        probs = probs,
        Q = Q,
        success = success,
        error_message = error_message
    )
end