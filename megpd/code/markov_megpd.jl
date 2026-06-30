using QuadGK, Roots, ProgressMeter
using SpecialFunctions
using Distributions

# function delta(r::Float64) 0.2 + 0.6*exp(-r/5) end
function delta(r::Float64) 0.8 end

# gd = Gamma(2, 10/3)
# function delta(r::Float64) 8 * pdf(gd, r) + 0.2 end
    
function egpd_cdf(x::Float64, kappa::Float64, sigma::Float64, xi::Float64)
    sigma = max(sigma, 1e-12)
    kappa = max(kappa, 1e-12)

    if xi == 0.0
        return (1 - exp(-x/sigma))^kappa
    end

    base = 1 - (1 + xi*x/sigma)^(-1/xi)
    base = max(base, 0.0)

    return base^kappa
end


function egpd_logdens(x::Float64, kappa::Float64, sigma::Float64, xi::Float64)
    if xi == 0.0
        return log(kappa) - log(sigma) - x/sigma + (kappa-1)*log(1 - exp(-x/sigma))
    end
    t = 1 + xi*x/sigma
    t <= 0.0 && return -Inf
    log_gpd_dens = -log(sigma) + (-1/xi - 1)*log(t)
    F = max(1 - exp(-1/xi * log(t)), 1e-300)
    return log(kappa) + log_gpd_dens + (kappa-1)*log(F)
end


function dmegpd_biv(
    x1::Float64,
    x2::Float64,
    kappa::Float64,
    sigma::Float64,
    xi::Float64
)

    # Support is x1,x2 > 0
    if x1 <= 0.0 || x2 <= 0.0
        return 0.0
    end

    r = x1 + x2

    if r <= 0.0
        return 0.0
    end

    # Radial component
    log_term_rad = egpd_logdens(r, kappa, sigma, xi)

    if !isfinite(log_term_rad)
        return 0.0
    end

    # Angular scale
    delta_r = max(delta(r), 0.01)

    # log(x1/x2)
    log_ratio = log(x1) - log(x2)

    # Log Gaussian density
    log_term_ang =
        -0.5 * log(2π) -
        log(delta_r) -
        0.5 * (log_ratio / delta_r)^2 +
        log(r) -
        log(x1) -
        log(x2)

    return exp(log_term_rad + log_term_ang)
end

function get_cdf_val(
    target_x::Float64,
    x_prev::Float64,
    kappa::Float64,
    sigma::Float64,
    xi::Float64,
    norm_const::Float64,
)

    if target_x <= 0
        return 0.0
    end

    val, _ = quadgk(
        x -> dmegpd_biv(
            x,
            x_prev,
            kappa,
            sigma,
            xi
        ),
        0.0,
        target_x,
    )

    return val / norm_const
end

function sample_conditional(
    x_prev::Float64,
    kappa::Float64,
    sigma::Float64,
    xi::Float64,
)

    # Normalizing constant
    norm_const, _ = quadgk(
        x -> dmegpd_biv(
            x,
            x_prev,
            kappa,
            sigma,
            xi
        ),
        1e-10,
        Inf,
    )

    # Uniform draw
    p_target = rand()

    lower = 1e-10
    upper = max(10*x_prev, 1.0)

    f(x) = get_cdf_val(
        x,
        x_prev,
        kappa,
        sigma,
        xi,
        norm_const,
    ) - p_target

    # Expand search interval until it brackets the root
    while f(upper) < 0
        upper *= 2
    end

    return find_zero(f, (lower, upper), Bisection())
end

function simulate_megpd_chain(
    n_steps::Int,
    kappa::Float64,
    sigma::Float64,
    xi::Float64,
    x0::Float64=1.0,
    burn_in_prop::Float64=0.10,
    show_progress::Bool=true,
)

    0 <= burn_in_prop < 1 || throw(ArgumentError(
        "burn_in_prop must be in 0,1."
    ))

    x = Vector{Float64}(undef, n_steps)
    x[1] = x0

    if show_progress
        p = Progress(n_steps; desc="Simulating")
    end

    for t in 2:n_steps

        x[t] = sample_conditional(
            x[t-1],
            kappa,
            sigma,
            xi
        )

        if show_progress
            next!(p)
        end
    end

    burn_in = floor(Int, n_steps * burn_in_prop)

    return (
        full_chain=x,
        final_chain=x[(burn_in+1):end],
        burn_in=burn_in,
    )
end

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

clust_res = clust_sim(3.0, 1000, 2.0, 1.0, 0.2)
boot = bootstrap_clust(clust_res.raw[1], clust_res.raw[2])