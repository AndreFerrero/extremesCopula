# markov_chain.jl
using QuadGK, Roots, Distributions, Interpolations
using DataFrames
using Random

# ------------------------------------------------------------------
# 1. COPULA DENSITIES
# ------------------------------------------------------------------

"""Bivariate Gaussian copula density via the conditional approach."""
function d_gaussian_copula(u, v, rho)
    # Transform to normal scores
    x = quantile(Normal(), u)
    y = quantile(Normal(), v)
    # Bivariate normal density / product of marginal normals
    z = (x^2 - 2*rho*x*y + y^2) / (1 - rho^2)
    return exp(-z/2 + (x^2 + y^2)/2) / sqrt(1 - rho^2)
end

"""Bivariate Gumbel copula density."""
function d_gumbel_copula(u, v, alpha)
    # Avoid boundary issues
    u = clamp(u, 1e-10, 1 - 1e-10)
    v = clamp(v, 1e-10, 1 - 1e-10)
    lu = -log(u); lv = -log(v)
    s  = (lu^alpha + lv^alpha)^(1/alpha)
    # Log-scale for numerical stability
    log_c = (alpha - 1)*log(lu) + (alpha - 1)*log(lv) -
            (lu^alpha + lv^alpha)^(1/alpha) +
            log(s^(1/alpha) + alpha - 1) -
            (1 + 1/alpha)*log(lu^alpha + lv^alpha) -
            log(u) - log(v)
    return exp(log_c)
end

# ------------------------------------------------------------------
# 2. MIXTURE DENSITY
# ------------------------------------------------------------------

pi_weight(u, v, theta) = (u * v)^theta

function c_star_unnorm(u, v, theta, rho, alpha)
    w   = pi_weight(u, v, theta)
    d_t = d_gaussian_copula(u, v, rho)
    d_b = d_gumbel_copula(u, v, alpha)
    return w * d_t + (1 - w) * d_b
end

# ------------------------------------------------------------------
# 3. PRECOMPUTATION (K, marginal CDF spline)
# ------------------------------------------------------------------

struct ModelPrecomp
    K          :: Float64
    u_star_grid:: Vector{Float64}
    F_grid     :: Vector{Float64}
    theta      :: Float64
    rho        :: Float64
    alpha      :: Float64
    F_spline   :: Any   # add spline object
    Finv_spline:: Any   # inverse spline
end

function build_precomp(theta, rho, alpha; grid_size=100)

    K, _ = quadgk(v -> quadgk(u -> c_star_unnorm(u, v, theta, rho, alpha),
                               0.001, 0.999)[1],
                   0.001, 0.999, rtol=1e-5)

    f_u_star(s) = quadgk(v -> c_star_unnorm(s, v, theta, rho, alpha),
                          0.001, 0.999)[1] / K

    u_grid = range(0.001, 0.999, length=grid_size) |> collect
    F_vals = similar(u_grid)

    for (i, lim) in enumerate(u_grid)
        F_vals[i], _ = quadgk(f_u_star, 0.001, lim, rtol=1e-5)
    end

    F_spline = LinearInterpolation(u_grid, F_vals, extrapolation_bc=Flat())
    Finv_spline = LinearInterpolation(F_vals, u_grid, extrapolation_bc=Flat())

    return ModelPrecomp(K, u_grid, F_vals, theta, rho, alpha, F_spline, Finv_spline)
end

# ------------------------------------------------------------------
# 4. SPLINE HELPERS  (linear interp fallback — swap for DataInterpolations)
# ------------------------------------------------------------------

# function interp1d(x, y, xi)
#     # Simple monotone linear interpolation (replace with spline if needed)
#     i = searchsortedfirst(x, xi) - 1
#     i = clamp(i, 1, length(x)-1)
#     t = (xi - x[i]) / (x[i+1] - x[i])
#     return y[i] + t * (y[i+1] - y[i])
# end

# F_u_star(p::ModelPrecomp, s)    = interp1d(p.u_star_grid, p.F_grid, s)
# F_u_star_inv(p::ModelPrecomp, u) = interp1d(p.F_grid, p.u_star_grid, u)

F_u_star(p::ModelPrecomp, s) = p.F_spline(s)
F_u_star_inv(p::ModelPrecomp, u) = p.Finv_spline(u)

# ------------------------------------------------------------------
# 5. TRANSITION KERNEL
# ------------------------------------------------------------------

function cond_cdf_latent(target_v, current_u_star, p::ModelPrecomp)
    num, _ = quadgk(v -> c_star_unnorm(current_u_star, v, p.theta, p.rho, p.alpha),
                     0.001, target_v, rtol=1e-6)
    den, _ = quadgk(v -> c_star_unnorm(current_u_star, v, p.theta, p.rho, p.alpha),
                     0.001, 0.999, rtol=1e-6)
    return num / den  # K cancels
end

function sample_next_step(u_t, p::ModelPrecomp)
    u_star_t    = F_u_star_inv(p, u_t)
    w           = rand()
    # Inverse-transform: find v* such that CDF(v*) = w
    f(v)        = cond_cdf_latent(v, u_star_t, p) - w
    u_star_next = find_zero(f, (0.001, 0.999), Bisection())
    return F_u_star(p, u_star_next)
end

# ------------------------------------------------------------------
# 6. SIMULATION ENTRY POINT
# ------------------------------------------------------------------

function simulate_chain(n_steps, theta, rho, alpha; seed=nothing)
    isnothing(seed) || Random.seed!(seed)
    println("Building precomputation...")
    p = build_precomp(theta, rho, alpha)
    println("  K = $(p.K)")

    chain = Vector{Float64}(undef, n_steps)
    chain[1] = rand()
    println("Simulating $n_steps steps...")
    for t in 1:(n_steps-1)
        chain[t+1] = sample_next_step(chain[t], p)
        t % 100 == 0 && println("  Step $t done")
    end
    return chain
end