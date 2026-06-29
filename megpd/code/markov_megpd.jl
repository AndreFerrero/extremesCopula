using QuadGK, Roots, ProgressMeter

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

using SpecialFunctions

function dmegpd_biv(
    x1::Float64,
    x2::Float64,
    kappa::Float64,
    sigma::Float64,
    xi::Float64,
    delta_func::F
) where {F<:Function}

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
    δ = max(delta_func(r), 0.01)

    # log(x1/x2)
    log_ratio = log(x1) - log(x2)

    # Log Gaussian density
    log_term_ang =
        -0.5 * log(2π) -
        log(δ) -
        0.5 * (log_ratio / δ)^2 +
        log(r) -
        log(x1) -
        log(x2)

    return exp(log_term_rad + log_term_ang)
end

function conditional_pdf_x(
    x,
    x_prev,
    kappa,
    sigma,
    xi,
    delta_func,
)
    return dmegpd_biv(
        x,
        x_prev,
        kappa,
        sigma,
        xi,
        delta_func,
    )
end

function get_cdf_val(
    target_x,
    x_prev,
    kappa,
    sigma,
    xi,
    delta_func,
    norm_const,
)

    if target_x <= 0
        return 0.0
    end

    val, _ = quadgk(
        x -> conditional_pdf_x(
            x,
            x_prev,
            kappa,
            sigma,
            xi,
            delta_func,
        ),
        0.0,
        target_x,
    )

    return val / norm_const
end

function sample_conditional(
    x_prev,
    kappa,
    sigma,
    xi,
    delta_func,
)

    # Normalizing constant
    norm_const, _ = quadgk(
        x -> conditional_pdf_x(
            x,
            x_prev,
            kappa,
            sigma,
            xi,
            delta_func,
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
        delta_func,
        norm_const,
    ) - p_target

    # Expand search interval until it brackets the root
    while f(upper) < 0
        upper *= 2
    end

    return find_zero(f, (lower, upper), Bisection())
end

using ProgressMeter

function simulate_megpd_chain(
    n_steps::Int,
    kappa,
    sigma,
    xi,
    delta_func;
    x0 = 1.0,
    burn_in_prop = 0.10,
    show_progress = true,
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
            xi,
            delta_func,
        )

        if show_progress
            next!(p)
        end
    end

    burn_in = floor(Int, n_steps * burn_in_prop)

    return (
        full_chain = x,
        final_chain = x[(burn_in+1):end],
        burn_in = burn_in,
    )
end