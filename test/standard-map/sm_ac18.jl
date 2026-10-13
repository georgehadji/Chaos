# Criterion 18, map-level part only: pure steps are deterministic across serial and threaded evaluation and allocate nothing.
# CPU kernel vs plain loop, 1 vs N threads of propagate, GPU vs CPU and the double-double GPU gate need Tangent.propagate
# or a GPU: tangent-core suite and stage 04 (see the test plan).

# run n steps for every start state, serially or with Threads.@threads (whatever thread count the process has)
function sm_run_serial(sys, Z, n)
    out = similar(Z)
    for j in eachindex(Z)
        z = Z[j]
        for _ in 1:n
            z = forward(sys, z)
        end
        out[j] = z
    end
    return out
end

function sm_run_threaded(sys, Z, n)
    out = similar(Z)
    Threads.@threads for j in eachindex(Z)
        z = Z[j]
        for _ in 1:n
            z = forward(sys, z)
        end
        out[j] = z
    end
    return out
end

# bytes allocated by one call of each pure step after a warm-up call (inside a function, not at testset global scope);
# wrap is excluded: it belongs to the display-only wrapped mode
function sm_step_allocations(sys, z)
    forward(sys, z); inverse(sys, z); tangent(sys, z); inverse_tangent(sys, z)
    a1 = @allocated forward(sys, z)
    a2 = @allocated inverse(sys, z)
    a3 = @allocated tangent(sys, z)
    a4 = @allocated inverse_tangent(sys, z)
    return (a1, a2, a3, a4)
end

@testset "AC18 map-level reproducibility" begin
    # 10^5 seeded orbits at K in {1, 5}, n = n_H(K): serial loop and threaded loop agree bitwise (isequal), Float64
    rng = StableRNG(SM_SEED + 18001)
    lim = 4 * Float64(pi)
    Z = [SVector((2 * rand(rng) - 1) * lim, (2 * rand(rng) - 1) * lim) for _ in 1:100_000]
    for K in (1.0, 5.0)
        sys = StandardMap(K)
        n = sm_nH(K)
        @test isequal(sm_run_serial(sys, Z, n), sm_run_threaded(sys, Z, n))
    end
    # zero allocation of every pure step after warm-up, Float64 and Double64
    for T in (Float64, Double64)
        @test all(==(0), sm_step_allocations(StandardMap(T(5)), SVector(T(1.25), T(-0.5))))
    end
end

@testset "SC1 non-finite input" begin
    # scope: non-finite input gives non-finite output and is not an error (spec Scope, Domain)
    sys = StandardMap(1.0)
    for bad in (NaN, Inf, -Inf)
        zf = forward(sys, SVector(bad, 0.5))
        @test !any(isfinite, zf)
        zp = forward(sys, SVector(0.5, bad))
        @test !any(isfinite, zp)
        @test !all(isfinite, tangent(sys, SVector(bad, 0.5)))
    end
end
