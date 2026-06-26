# markov_fitting.jl
using QuadGK
using Optim
using ForwardDiff
using Interpolations        # replaces Interpolations.jl
using Distributions
using LinearAlgebra
using CSV, DataFrames

const EPS_DEFAULT = 1e-6

# ================================================================
# 1. COPULA COMPONENTS
# ================================================================

function d_gaussian_copula(u::Float64, v::Float64, rho::Float64)
    rho = clamp(rho, -0.999, 0.999)

    u = clamp(u, 1e-12, 1 - 1e-12)
    v = clamp(v, 1e-12, 1 - 1e-12)

    x = quantile(Normal(), u)
    y = quantile(Normal(), v)

    det = 1 - rho^2
    quad = (x^2 - 2*rho*x*y + y^2) / det

    # numerically stable exponent form
    return exp(-0.5 * quad + 0.5*(x^2 + y^2)) / sqrt(det)
end

function d_gumbel_copula(u::Float64, v::Float64, alpha::Float64)
    α = max(alpha, 1.0 + 1e-6)

    u = clamp(u, 1e-12, 1 - 1e-12)
    v = clamp(v, 1e-12, 1 - 1e-12)

    lu = max(-log(u), 1e-12)
    lv = max(-log(v), 1e-12)

    a = lu^α
    b = lv^α
    s = (a + b)^(1/α)

    # protected logs
    log_lu = log(lu)
    log_lv = log(lv)
    log_sum = log(a + b)

    log_c =
        (α - 1) * (log_lu + log_lv) +
        log(s + α - 1) -
        (1 + 1/α) * log_sum -
        log(u) - log(v) -
        s

    return exp(clamp(log_c, -700, 700))
end

pi_weight(u, v, theta) = begin
    u = clamp(u, 1e-12, 1 - 1e-12)
    v = clamp(v, 1e-12, 1 - 1e-12)
    return clamp((u * v)^theta, 0.0, 1.0)
end

function c_star_unnorm(u::Float64, v::Float64,
                       theta::Float64, rho::Float64, alpha::Float64)

    w  = pi_weight(u, v, theta)

    dt = d_gaussian_copula(u, v, rho)
    db = d_gumbel_copula(u, v, alpha)

    val = w * dt + (1 - w) * db

    return isfinite(val) && val > 0 ? val : 1e-300
end


# ================================================================
# 2. EGPD MARGINS
# ================================================================

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


# ================================================================
# 3. MARKOV COPULA BUILD  (always Float64 — no AD through here)
# ================================================================

struct MarkovCopula{I1, I2}
    theta  :: Float64
    rho    :: Float64
    alpha  :: Float64
    K      :: Float64
    grid   :: Vector{Float64}
    Fgrid  :: Vector{Float64}
    F_fwd  :: I1    # u*  → u   (FritschCarlson ≈ monoH.FC)
    F_inv  :: I2    # u   → u*
end

function build_markov_copula(theta::Float64, rho::Float64, alpha::Float64;
                              grid_size::Int = 120,
                              eps::Float64   = EPS_DEFAULT,
                              rtol::Float64  = 1e-5)

    K_val, _ = quadgk(eps, 1-eps, rtol=rtol) do v
        quadgk(eps, 1-eps, rtol=rtol) do u
            c_star_unnorm(u, v, theta, rho, alpha)
        end[1]
    end

    u_grid = collect(range(eps, 1-eps, length=grid_size))

    f_grid = [
        quadgk(eps, 1-eps, rtol=rtol) do v
            c_star_unnorm(s, v, theta, rho, alpha)
        end[1] / K_val
        for s in u_grid
    ]

    du     = diff(u_grid)
    F_vals = vcat(0.0, cumsum(0.5 .* (f_grid[1:end-1] .+ f_grid[2:end]) .* du))
    F_vals ./= F_vals[end]
    F_vals  = clamp.(F_vals, 0.0, 1.0)

    for i in 2:length(F_vals)
        F_vals[i] = max(F_vals[i], F_vals[i-1] + 1e-12)
    end

    # FritschCarlsonMonotonicInterpolation = R's monoH.FC
    # extrapolate=true handles boundary queries gracefully
    F_fwd = interpolate(u_grid, F_vals, SteffenMonotonicInterpolation())

    F_inv = interpolate(F_vals, u_grid, SteffenMonotonicInterpolation())

    return MarkovCopula(theta, rho, alpha, K_val, u_grid, F_vals, F_fwd, F_inv)
end


# ================================================================
# 4.  TRANSITION LOG-COPULA DENSITY
# ================================================================

function log_copula_density(mc::MarkovCopula, u::Float64, v::Float64)
    u = clamp(u, EPS_DEFAULT, 1-EPS_DEFAULT)
    v = clamp(v, EPS_DEFAULT, 1-EPS_DEFAULT)

    u_star = clamp(mc.F_inv(u), EPS_DEFAULT, 1-EPS_DEFAULT)
    v_star = clamp(mc.F_inv(v), EPS_DEFAULT, 1-EPS_DEFAULT)

    cstar = c_star_unnorm(u_star, v_star, mc.theta, mc.rho, mc.alpha) / mc.K

    fu = quadgk(EPS_DEFAULT, 1-EPS_DEFAULT) do t
        c_star_unnorm(u_star, t, mc.theta, mc.rho, mc.alpha)
    end[1] / mc.K

    fv = quadgk(EPS_DEFAULT, 1-EPS_DEFAULT) do t
        c_star_unnorm(v_star, t, mc.theta, mc.rho, mc.alpha)
    end[1] / mc.K

    return log(max(cstar, floatmin(Float64))) -
           log(max(fu,    floatmin(Float64))) -
           log(max(fv,    floatmin(Float64)))
end


# ================================================================
# 5.  NEGATIVE LOG-LIKELIHOOD
#     par is always Float64 here — the copula precomputation is
#     intentionally kept outside the AD graph (see Section 6).
# ================================================================

function negloglik(par::Vector{Float64}, x::Vector{Float64};
                   grid_size::Int    = 120,
                   eps::Float64      = EPS_DEFAULT,
                   penalty::Float64  = 1e12)

    theta = exp(par[1])
    rho   = tanh(par[2])
    alpha = 1.0 + exp(par[3])
    kappa, sigma, xi = par[4], par[5], par[6]

    (kappa <= 0.0 || sigma <= 0.0) && return penalty

    logf = egpd_logdens.(x, kappa, sigma, xi)
    any(!isfinite, logf) && return penalty

    u = clamp.(egpd_cdf.(x, kappa, sigma, xi), eps, 1-eps)

    mc = try
        build_markov_copula(theta, rho, alpha;
                             grid_size=grid_size, eps=eps)
    catch e
        @warn "copula failed" e
        return penalty
    end

    logc = [log_copula_density(mc, u[t], u[t+1]) for t in 1:length(u)-1]
    any(!isfinite, logc) && return penalty

    return -(sum(logf) + sum(logc))
end


# ================================================================
# 6.  FITTING
#     The Hessian is computed via finite differences on negloglik
#     (not ForwardDiff) because build_markov_copula contains
#     quadrature with Float64-only internals that block AD.
#     For a 6-parameter model this is 2*6^2 = 72 extra likelihood
#     evaluations — fast and numerically stable.
# ================================================================

function hessian_fd(f, p; h=1e-4)
    n = length(p)
    H = zeros(n, n)
    f0 = f(p)
    for i in 1:n, j in i:n
        ei = zeros(n); ei[i] = h
        ej = zeros(n); ej[j] = h
        if i == j
            H[i,i] = (f(p .+ ei) - 2*f0 + f(p .- ei)) / h^2
        else
            H[i,j] = H[j,i] =
                (f(p .+ ei .+ ej) - f(p .+ ei .- ej) -
                 f(p .- ei .+ ej) + f(p .- ei .- ej)) / (4*h^2)
        end
    end
    return H
end

function fit_markov_egpd(x::Vector{Float64};
                          start_theta::Float64  = 1.5,
                          start_rho::Float64    = 0.3,
                          start_alpha::Float64  = 1.2,
                          start_margin::Vector{Float64} = [2.0, 1.0, 0.2],
                          grid_size::Int        = 120,
                          iterations::Int       = 500,
                          show_trace::Bool      = true)

    par0 = [
        log(start_theta),
        atanh(clamp(start_rho, -0.999, 0.999)),
        log(start_alpha - 1.0),
        start_margin...
    ]

    obj(p) = negloglik(p, x; grid_size=grid_size)

    opt = optimize(obj, par0, LBFGS(),
                   Optim.Options(iterations=iterations,
                                 show_trace=show_trace,
                                 g_tol=1e-6))

    par_hat = Optim.minimizer(opt)

    # finite-difference Hessian (safe with quadrature internals)
    H = hessian_fd(obj, par_hat)

    vcov_mat = nothing
    se_vec   = nothing
    try
        invH = inv(Symmetric(H))
        d    = diag(invH)
        if all(isfinite, d) && all(d .>= 0)
            vcov_mat = invH
            se_vec   = sqrt.(d)
        end
    catch e
        @warn "Hessian inversion failed: $e"
    end

    return (
        par         = par_hat,
        estimates   = (
            theta = exp(par_hat[1]),
            rho   = tanh(par_hat[2]),
            alpha = 1.0 + exp(par_hat[3]),
            kappa = par_hat[4],
            sigma = par_hat[5],
            xi    = par_hat[6],
        ),
        convergence = Optim.converged(opt),
        loglik      = -Optim.minimum(opt),
        hessian     = H,
        vcov        = vcov_mat,
        se          = se_vec,
    )
end