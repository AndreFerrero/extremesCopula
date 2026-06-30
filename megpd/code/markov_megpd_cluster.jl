include("markov_megpd.jl")

function clust_sim(
    u::Float64,
    n_sims::Int,
    kappa::Float64,
    sigma::Float64,
    xi::Float64;
    m::Int = 3,
    max_steps::Int = 100,
)
    gpd = GeneralizedPareto(sigma, xi)

    n_counts = Vector{Int}(undef, n_sims)
    max_consecutive_runs = Vector{Int}(undef, n_sims)

    for sim in 1:n_sims

        # Initial exceedance
        x_current = u + rand(gpd)

        exceedances = 1

        current_run = 1
        longest_run = 1

        consecutive_lows = 0

        for _ in 2:max_steps

            x_next = sample_conditional(x_current, kappa, sigma, xi)

            if x_next > u

                exceedances += 1
                current_run += 1
                consecutive_lows = 0

            else

                longest_run = max(longest_run, current_run)
                current_run = 0
                consecutive_lows += 1

            end

            if consecutive_lows >= m
                break
            end

            x_current = x_next
        end

        longest_run = max(longest_run, current_run)

        n_counts[sim] = exceedances
        max_consecutive_runs[sim] = longest_run
    end

    #
    # ---- Total cluster size distribution ----
    #

    max_n = maximum(n_counts)

    theta = zeros(max_n + 1)

    for i in 1:max_n
        theta[i] = count(>=(i), n_counts) / n_sims
    end

    ext_index_theta = theta[1]

    pi_dist = [
        (theta[i] - theta[i + 1]) / ext_index_theta
        for i in 1:max_n
    ]

    #
    # ---- Longest consecutive run distribution ----
    #

    max_r = maximum(max_consecutive_runs)

    thetaC = zeros(max_r + 1)

    for i in 1:max_r
        thetaC[i] = count(>=(i), max_consecutive_runs) / n_sims
    end

    piC_dist = [
        (thetaC[i] - thetaC[i + 1]) / thetaC[1]
        for i in 1:max_r
    ]

    avg_run = sum(i * piC_dist[i] for i in eachindex(piC_dist))
    thetaC_index = 1 / avg_run

    return (
        pi_dist = pi_dist,
        theta = ext_index_theta,
        piC_dist = piC_dist,
        thetaC = thetaC_index,
        raw = (n_counts, max_consecutive_runs)
    )
end

function bootstrap_clust(n_counts, max_consecutive_runs; B::Int = 500)

    n_sims = length(n_counts)

    theta_boot = Float64[]
    thetaC_boot = Float64[]

    for _ in 1:B

        idx = rand(1:n_sims, n_sims)

        nc = n_counts[idx]
        mr = max_consecutive_runs[idx]

        #
        # theta (extremal index part)
        #
        max_n = maximum(nc)
        theta = zeros(max_n + 1)

        for i in 1:max_n
            theta[i] = count(>=(i), nc) / n_sims
        end

        theta1 = theta[1]

        #
        # thetaC
        #
        max_r = maximum(mr)
        thetaC = zeros(max_r + 1)

        for i in 1:max_r
            thetaC[i] = count(>=(i), mr) / n_sims
        end

        push!(theta_boot, theta1)
        push!(thetaC_boot, thetaC[1])
    end

    return (
        theta_mean = mean(theta_boot),
        theta_sd   = std(theta_boot),
        thetaC_mean = mean(thetaC_boot),
        thetaC_sd   = std(thetaC_boot),
    )
end
