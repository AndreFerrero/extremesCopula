using QuadGK
using Printf
using Statistics

f(x) = sin(x^2)

function batch_integrate(N::Int)
    results = Vector{Float64}(undef, N)

    for i in 1:N
        results[i] = quadgk(f, 0.0, 10.0)[1]
    end

    return results
end

N = 100000

t0 = time()

results = batch_integrate(N)

t1 = time()

@printf("Julia result\n")

@printf("Mean: %.10f\n", mean(results))
@printf("Elapsed time: %.6f seconds\n", t1 - t0)