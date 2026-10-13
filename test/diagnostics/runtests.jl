# Diagnostics suite (RED, stage 02). Spec: workspaces/module-build/stages/01-spec/output/diagnostics-spec.md
# API contract: M2 Batch A api-batch-a.md. Helpers are prefixed dg_. Criterion 7 and the GPU part of 11 are stage 04.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const dg_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const dg_DET_TOL = 4.0                    # _config/numerics-defaults.toml key `det_residual_tol`
const dg_PAIR_TOL = 16.0                  # _config/numerics-defaults.toml key `pairing_residual_tol`
const dg_ORTH_TOL = 1.0e-6                # absolute orthogonality tolerance for tests that do not probe it
const dg_TYPES = (Float64, Double64)      # Float32 dropped (D24)

# ---------------------------------------------------------------- helpers (no Chaoslab names at load time)

const dg_FAULTS_LOADED = Ref(false)
function dg_load_faults()                 # defines dg_FaultMap <: AbstractMap and dg_Renorm
    dg_FAULTS_LOADED[] || (include(joinpath(@__DIR__, "faultmaps.jl")); dg_FAULTS_LOADED[] = true)
    return nothing
end

# Run f; return the NumericalHealthError it throws, or nothing. Other exceptions propagate.
function dg_catch(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError || rethrow()
        return e
    end
end

# Powers of two strictly bracketing a threshold: raise side [2 theta, 4 theta), quiet side (theta/4, theta/2].
dg_above(theta) = [exp2(j) for j in -130:20 if 2theta <= exp2(j) < 4theta]
dg_below(theta) = [exp2(j) for j in -130:20 if theta / 4 < exp2(j) <= theta / 2]

# Reference orthonormalisation (modified Gram-Schmidt) of a 2 x 2 matrix, generic in T.
function dg_mgs(W::SMatrix{2,2,T}) where {T}
    q1 = W[:, 1] / sqrt(sum(abs2, W[:, 1]))
    v = W[:, 2] - dot(q1, W[:, 2]) * q1
    return hcat(q1, v / sqrt(sum(abs2, v)))
end

# ---------------------------------------------------------------- AC1

# AC1: integer symplectic maps have exact arithmetic, so eps_det = eps_symp = 0 at every step (wrap keeps the cat orbit finite).
function dg_ac1(::Type{T}) where {T}
    mats = (shear = (1, 0, 1, 1), cat = (2, 1, 1, 1), rot4 = (0, -1, 1, 0), rot6 = (1, 1, -1, 0))  # column-major
    return [dg_ac1_one(T, name, ent) for (name, ent) in pairs(mats)]
end
function dg_ac1_one(::Type{T}, name, ent) where {T}
    begin
        sys = LinearSymplecticMap(SMatrix{2,2,T}(ent...))
        mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
        st = reducer_init(mon, T, Val(2), Val(2))
        W = SMatrix{2,2,T}(1, 0, 0, 1)
        z = SVector{2,T}(1, 2)
        maxd = zero(T)
        maxs = zero(T)
        for j in 1:1000
            M = tangent(sys, z)
            maxd = max(maxd, eps_det(M))
            maxs = max(maxs, eps_symp(M))
            check_residual("det", eps_det(M), det_threshold(dg_DET_TOL, T, M); step = j, orbit = 1)
            st, W = reducer_update(mon, st, j, z, W)
            z = wrap(forward(sys, z))
        end
        (; name, maxd, maxs, fin = reducer_final(mon, st))
    end
end

# ---------------------------------------------------------------- AC2

# AC2: rotation by 0.7311, M <= 1 so theta_det <= 8 eps(T); eps_orth after reference re-orthonormalisation <= 8 eps(T).
function dg_ac2(::Type{T}, n) where {T}
    sys = LinearSymplecticMap(rotation_block(T, T(0.7311)))
    z = SVector{2,T}(1, 0)
    Q = SMatrix{2,2,T}(1, 0, 0, 1)
    maxdet = zero(T)
    maxorth = zero(T)
    viol = 0
    for _ in 1:n
        M = tangent(sys, z)
        d = eps_det(M)
        maxdet = max(maxdet, d)
        d <= det_threshold(dg_DET_TOL, T, M) || (viol += 1)
        z = forward(sys, z)
        Q = dg_mgs(M * Q)
        maxorth = max(maxorth, eps_orth(Q))
    end
    return (; maxdet, maxorth, viol)
end

# ---------------------------------------------------------------- AC3

# AC3: cat map, reference Householder QR per step, log|R_ii| accumulated with CompensatedSum, lambda_i(n) = sum / n.
# Q orthogonal => singular values of A Q equal those of A => G = sigma_1/sigma_2 = mu^2 at every window.
function dg_ac3(n)
    A = [2.0 1.0; 1.0 1.0]
    sv = svdvals(A)
    G = sv[1] / sv[2]
    Mmax = maximum(abs, A)
    theta = sum_threshold(dg_PAIR_TOL, Float64, G, 1, Mmax)
    Q = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]          # arbitrary orthonormal start (any basis works, I2)
    a1 = CompensatedSum{Float64}()
    a2 = CompensatedSum{Float64}()
    viol = 0
    pairgap = 0.0
    lam = SVector(0.0, 0.0)
    for k in 1:n
        F = qr(A * Q)
        Q = Matrix(F.Q)
        a1 = add(a1, log(abs(F.R[1, 1])))
        a2 = add(a2, log(abs(F.R[2, 2])))
        lam = SVector(total(a1) / k, total(a2) / k)
        eps_sum(lam) <= theta || (viol += 1)
        pairgap = max(pairgap, abs(eps_pair(lam) - eps_sum(lam)))
    end
    return (; G, theta, viol, pairgap, lam)
end

# ---------------------------------------------------------------- AC4

# AC4: pairing index i <-> 2N+1-i; the planted offsets d1, d2 are recovered, adjacent pairing would give ~1.27.
function dg_ac4(::Type{T}) where {T}
    mu = (3 + sqrt(T(5))) / 2
    L = log(mu)
    d = T(2)^-20
    lam1 = SVector{2,T}(L + d, -L)
    a = T(96) / 100
    b = T(31) / 100
    d1 = T(2)^-20
    d2 = T(2)^-18
    lam2 = SVector{4,T}(a, b, -b + d2, -a + d1)
    return (; L, d, d1, d2, lam1, lam2,
        pair1 = eps_pair(lam1), sum1 = eps_sum(lam1), pair2 = eps_pair(lam2), sum2 = eps_sum(lam2))
end

# ---------------------------------------------------------------- AC5 helpers

# AC5a: A_delta = [[2 + delta, 1], [1, 1]] has det exactly 1 + delta; Renorm keeps W finite over long runs.
function dg_ac5a(::Type{T}) where {T}
    A0 = SMatrix{2,2,T}(2, 1, 1, 1)
    theta = Float64(det_threshold(dg_DET_TOL, T, A0))
    x0 = SVector{2,T}(0.5, 0.25)
    W0 = SMatrix{2,1,T}(1, 0)
    run(delta, N) = begin
        sys = dg_FaultMap(SMatrix{2,2,T}(2 + T(delta), 1, 1, 1))
        mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
        dg_catch(() -> propagate_orbit(sys, x0, W0, N, 1, (dg_Renorm(), mon); orbit = 7))
    end
    above = [(delta, run(delta, 1)) for delta in dg_above(theta)]
    below = [(delta, run(delta, 10_000)) for delta in dg_below(theta)]
    return (; theta, above, below)
end

# AC5b: basis with one vector scaled by (1 + delta): eps_orth = 2 delta + delta^2; run tolerance theta_o is absolute.
function dg_ac5b(::Type{T}, theta_o) where {T}
    sys = dg_FaultMap(SMatrix{2,2,T}(1, 0, 0, 1))
    mon = ResidualMonitor(sys, dg_DET_TOL, theta_o)
    x = SVector{2,T}(0.5, 0.25)
    run(delta) = begin
        Q = SMatrix{2,2,T}(1 + T(delta), 0, 0, 1)
        st = reducer_init(mon, T, Val(2), Val(2))
        (eps_orth(Q), dg_catch(() -> reducer_update(mon, st, 5, x, Q)))
    end
    # residual 2 delta in [2 theta_o, 4 theta_o) means delta in [theta_o, 2 theta_o); quiet side delta = theta_o / 4
    above = [(delta, run(delta)...) for delta in dg_above(theta_o) .|> (r -> r / 2)]
    below = [(delta, run(delta)...) for delta in dg_below(theta_o) .|> (r -> r / 2)]
    return (; above, below)
end

# AC5c: delta added to lambda_1 (N = 1) gives eps_sum = eps_pair = delta; sum_threshold with G_max = mu^2.
function dg_ac5c(::Type{T}) where {T}
    mu = (3 + sqrt(T(5))) / 2
    L = log(mu)
    G = mu^2
    Mmax = T(2)
    theta = sum_threshold(dg_PAIR_TOL, T, G, 1, Mmax)
    mon = ResidualMonitor(StandardMap(T(2)), dg_DET_TOL, dg_ORTH_TOL)
    run(delta) = begin
        lam = SVector{2,T}(L + T(delta), -L)
        st = reducer_init(mon, T, Val(2), Val(2))
        e_mon = dg_catch(() -> monitor_exponents(st, lam, G, 1; pairing_tol = dg_PAIR_TOL, tau = 1,
            M_max = Mmax, step = 3, orbit = 9))
        e_pair = dg_catch(() -> check_residual("pairing", eps_pair(lam), theta; step = 3, orbit = 9))
        (; delta, es = eps_sum(lam), ep = eps_pair(lam), e_mon, e_pair)
    end
    thetaf = Float64(theta)
    return (; L, above = [run(d) for d in dg_above(thetaf)], below = [run(d) for d in dg_below(thetaf)])
end

# AC5d: 4 x 4 block map, one block det 1 + delta; threshold_symp = tol eps N^2 (1 + M)^2. Plus a det-1 non-symplectic map via the monitor.
function dg_ac5d(::Type{T}) where {T}
    B1 = SMatrix{2,2,T}(2, 1, 1, 1)
    A0 = block_diagonal(B1, B1)
    theta = symp_threshold(dg_DET_TOL, T, A0)
    run(delta) = begin
        A = block_diagonal(B1, SMatrix{2,2,T}(2 + T(delta), 1, 1, 1))
        (; delta, es = eps_symp(A),
            err = dg_catch(() -> check_residual("symplecticity", eps_symp(A), symp_threshold(dg_DET_TOL, T, A);
                step = 1, orbit = 3)))
    end
    thetaf = Float64(theta)
    D = SMatrix{4,4,T}(2, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)   # det 1, (q1,p1) entry of D' Omega D is 2
    mon = ResidualMonitor(dg_FaultMap(D), dg_DET_TOL, dg_ORTH_TOL)
    st = reducer_init(mon, T, Val(4), Val(1))
    e_mon = dg_catch(() -> reducer_update(mon, st, 4, SVector{4,T}(0.1, 0.2, 0.3, 0.4), SMatrix{4,1,T}(1, 0, 0, 0)))
    return (; thetaf, above = [run(d) for d in dg_above(thetaf)], below = [run(d) for d in dg_below(thetaf)], e_mon)
end

# AC5e: NaN or Inf raises with the right name for any threshold, including an infinite one.
function dg_ac5e(::Type{T}) where {T}
    x = SVector{2,T}(0.5, 0.25)
    W1 = SMatrix{2,1,T}(1, 0)
    ok = dg_FaultMap(SMatrix{2,2,T}(1, 0, 0, 1))
    cases = Tuple{String,String,Any}[]
    for what in ("det", "symplecticity", "orthogonality", "sum", "pairing"), v in (T(NaN), T(Inf))
        push!(cases, ("direct $what $v", what,
            dg_catch(() -> check_residual(what, v, T(Inf); step = 1, orbit = 1))))
    end
    for bad in (T(NaN), T(Inf))          # non-finite entry in DT
        mon = ResidualMonitor(dg_FaultMap(SMatrix{2,2,T}(bad, 1, 1, 1)), Inf, Inf)
        st = reducer_init(mon, T, Val(2), Val(1))
        push!(cases, ("DT $bad", "det", dg_catch(() -> reducer_update(mon, st, 1, x, W1))))
    end
    for bad in (T(NaN), T(Inf))          # non-finite entry in Q
        mon = ResidualMonitor(ok, Inf, Inf)
        st = reducer_init(mon, T, Val(2), Val(2))
        Q = SMatrix{2,2,T}(bad, 0, 0, 1)
        push!(cases, ("Q $bad", "orthogonality", dg_catch(() -> reducer_update(mon, st, 1, x, Q))))
    end
    for bad in (T(NaN), T(Inf))          # non-finite exponent
        mon = ResidualMonitor(ok, Inf, Inf)
        st = reducer_init(mon, T, Val(2), Val(2))
        lam = SVector{2,T}(bad, -0.5)
        push!(cases, ("lambda $bad", "sum", dg_catch(() ->
            monitor_exponents(st, lam, T(6.8), 1; pairing_tol = Inf, tau = 1, M_max = T(2), step = 1, orbit = 1))))
    end
    never_clipped = (isfinite(eps_orth(SMatrix{2,2,T}(T(NaN), 0, 0, 1))),
        isfinite(eps_sum(SVector{2,T}(T(NaN), 0))), isfinite(eps_pair(SVector{2,T}(T(Inf), 0))),
        isfinite(eps_det(SMatrix{2,2,T}(T(NaN), 1, 1, 1))))
    return (; cases, never_clipped)
end

# ---------------------------------------------------------------- AC6

# AC6: B1 bound eps_det <= 1.5 eps(T) M for the standard map along real orbits; starts uniform on [0, 2 pi)^2.
function dg_ac6(::Type{T}, K, rng, nstarts, nsteps) where {T}
    sys = StandardMap(T(K))
    ratio = 0.0
    viol = 0
    for _ in 1:nstarts
        z = SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng)))
        for _ in 1:nsteps
            M = tangent(sys, z)
            e = eps_det(M)
            ratio = max(ratio, Float64(e / (eps(T) * maximum(abs, M))))
            e <= det_threshold(dg_DET_TOL, T, M) || (viol += 1)
            z = forward(sys, z)
        end
    end
    return (; ratio, viol)
end

# ---------------------------------------------------------------- AC8

# AC8: eps_det of stored binary64 entries vs the same entries evaluated at 256 bits; they differ by one rounding <= u |det| < eps.
function dg_ac8()
    sys = StandardMap(5.0)
    rng = StableRNG(dg_SEED + 8)
    z = SVector(2pi * rand(rng), 2pi * rand(rng))
    worst = setprecision(BigFloat, 256) do
        w = 0.0
        for _ in 1:1000
            M = tangent(sys, z)
            Mb = BigFloat.(M)
            ref = abs(Mb[1, 1] * Mb[2, 2] - Mb[1, 2] * Mb[2, 1] - 1)
            w = max(w, abs(Float64(ref) - eps_det(M)))
            z = forward(sys, z)
        end
        w
    end
    return worst
end

# ---------------------------------------------------------------- AC9

# AC9: zero bytes after warm-up, measured inside functions on propagate_orbit and on the single update calls.
function dg_alloc_orbit(sys, mon, x0, W0)
    propagate_orbit(sys, x0, W0, 100, 1, (mon,))
    return @allocated propagate_orbit(sys, x0, W0, 100, 1, (mon,))
end
function dg_alloc_update(mon, st, x, W)
    reducer_update(mon, st, 1, x, W)
    return @allocated begin
        r = reducer_update(mon, st, 2, x, W)
        r
    end
end
function dg_alloc_exponents(mon, st, lam, G, M)
    monitor_exponents(st, lam, G, 1; pairing_tol = dg_PAIR_TOL, tau = 1, M_max = M, step = 1, orbit = 1)
    return @allocated begin
        r = monitor_exponents(st, lam, G, 2; pairing_tol = dg_PAIR_TOL, tau = 1, M_max = M, step = 2, orbit = 1)
        r
    end
end
function dg_ac9(::Type{T}) where {T}
    sys = StandardMap(T(2))
    mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
    x0 = SVector{2,T}(1, 2)
    st2 = reducer_init(mon, T, Val(2), Val(2))
    lam = SVector{2,T}(0.5, -0.5)
    return (orbit = dg_alloc_orbit(sys, mon, x0, SMatrix{2,1,T}(1, 0)),
        update = dg_alloc_update(mon, st2, x0, SMatrix{2,2,T}(1, 0, 0, 1)),
        exponents = dg_alloc_exponents(mon, st2, lam, T(2), T(2)))
end

# ---------------------------------------------------------------- AC10 / AC11

# Serial reference: the first error over orbits 1..n in id order (what a thread-count independent report must equal).
function dg_serial_first_error(sys, x0s, W0s, N, tau, reducers)
    for j in eachindex(x0s)
        e = dg_catch(() -> propagate_orbit(sys, x0s[j], W0s[j], N, tau, reducers; orbit = j))
        e === nothing || return e
    end
    return nothing
end

# AC10 series: propagate versus serial map over propagate_orbit, bitwise (x, W, reducer finals) for 64 orbits.
function dg_ac10_series()
    rng = StableRNG(dg_SEED + 10)
    sys = StandardMap(2.0)
    mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
    x0s = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:64]
    W0s = [(a = 2pi * rand(rng); SMatrix{2,1,Float64}(cos(a), sin(a))) for _ in 1:64]
    N = 100
    res = propagate(sys, x0s, N; k = 1, tau = 1, seed = dg_SEED + 10, reducers = (mon,), W0 = W0s)
    ser = map(1:64) do j
        x, W, st = propagate_orbit(sys, x0s[j], W0s[j], N, 1, (mon,); orbit = j)
        (x, W, (reducer_final(mon, st[1]),))
    end
    same = all(j -> res.x[j] === ser[j][1] && res.W[j] === ser[j][2] && res.finals[j] === ser[j][3], 1:64)
    return (; same, nontrivial = maximum(f[1].det_max for f in res.finals) > 0, nthreads = Threads.nthreads())
end

# AC10 faults: orbit i starts at (i, 0); the tangent has det 1 + 2^-20 only at the target orbits (first coordinate equal).
function dg_ac10_fault(targets)
    delta = 2.0^-20
    healthy = SMatrix{2,2,Float64}(1, 0, 0, 1)
    sys = dg_FaultMap(healthy, SMatrix{2,2,Float64}(1 + delta, 0, 0, 1), targets[1], targets[2])
    mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
    x0s = [SVector(Float64(i), 0.0) for i in 1:64]
    W0s = [SMatrix{2,1,Float64}(1, 0) for _ in 1:64]
    e_par = dg_catch(() -> propagate(sys, x0s, 10; k = 1, tau = 1, seed = dg_SEED + 10, reducers = (mon,), W0 = W0s))
    e_ser = dg_serial_first_error(sys, x0s, W0s, 10, 1, (mon,))
    return (; e_par, e_ser)
end

# AC11 item 1 on the CPU backend: integer maps through the kernel, residual maxima exactly zero.
function dg_ac11_item1(::Type{T}) where {T}
    mats = (shear = ((1, 0, 1, 1), 1000), cat = ((2, 1, 1, 1), 1000),
        rot4 = ((0, -1, 1, 0), 1000), rot6 = ((1, 1, -1, 0), 1000))
    return [dg_ac11_one(T, name, ent, N) for (name, (ent, N)) in pairs(mats)]
end
function dg_ac11_one(::Type{T}, name, ent, N) where {T}
    begin
        sys = LinearSymplecticMap(SMatrix{2,2,T}(ent...))
        mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
        # zero start is a fixed point of every linear map (lift cannot overflow); Renorm keeps W finite
        res = propagate(sys, [SVector{2,T}(0, 0)], N; k = 1, tau = 1, seed = dg_SEED + 11,
            reducers = (dg_Renorm(), mon), backend = :cpu, W0 = [SMatrix{2,1,T}(1, 0)])
        (; name, fin = res.finals[1][2])
    end
end

# AC11 item 5 on the CPU backend: det fault through propagate, both sides of the threshold.
function dg_ac11_item5()
    A0 = SMatrix{2,2,Float64}(2, 1, 1, 1)
    theta = Float64(det_threshold(dg_DET_TOL, Float64, A0))
    run(delta, N) = begin
        sys = dg_FaultMap(SMatrix{2,2,Float64}(2 + delta, 1, 1, 1))
        mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
        dg_catch(() -> propagate(sys, [SVector(0.5, 0.25)], N; k = 1, tau = 1, seed = dg_SEED + 11,
            reducers = (dg_Renorm(), mon), backend = :cpu, W0 = [SMatrix{2,1,Float64}(1, 0)]))
    end
    return (above = [(d, run(d, 1)) for d in dg_above(theta)], below = [(d, run(d, 1000)) for d in dg_below(theta)])
end

# AC11 item 6 on the CPU backend: standard-map panel through the kernel, 50 seeded starts per K, 10^4 steps.
function dg_ac11_item6()
    rng = StableRNG(dg_SEED + 11)
    out = []
    for K in (0.5, 0.9716, 2.0, 5.0, 10.0)
        sys = StandardMap(K)
        mon = ResidualMonitor(sys, dg_DET_TOL, dg_ORTH_TOL)
        x0s = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:50]
        W0s = [(a = 2pi * rand(rng); SMatrix{2,1,Float64}(cos(a), sin(a))) for _ in 1:50]
        err = nothing
        res = nothing
        err = dg_catch(() -> (res = propagate(sys, x0s, 10_000; k = 1, tau = 1, seed = dg_SEED + 11,
            reducers = (dg_Renorm(), mon), backend = :cpu, W0 = W0s)))
        dmax = err === nothing ? maximum(f[2].det_max for f in res.finals) : NaN
        push!(out, (; K, err, dmax))
    end
    return out
end

# ================================================================ testsets

@testset "diagnostics" begin

    @testset "AC1 integer maps exact" begin   # all arithmetic exact: eps_det = eps_symp = 0, no health error
        for T in dg_TYPES, r in dg_ac1(T)
            @test r.maxd == 0
            @test r.maxs == 0
            @test r.fin.det_max == 0
            @test r.fin.symp_max == 0
            @test r.fin.orth_max == 0
        end
    end

    @testset "AC2 rotation 0.7311" begin   # theta_det <= 8 eps (M <= 1); eps_orth <= 8 eps after reference re-orthonormalisation
        for T in dg_TYPES
            r = dg_ac2(T, 1_000_000)
            @test r.viol == 0
            @test r.maxdet <= 8 * eps(T)
            @test r.maxorth <= 8 * eps(T)
        end
    end

    @testset "AC3 cat map reference QR" begin   # eps_sum(n) <= theta_sum at every n, G_max = mu^2 = 6.854, pairing = sum for N = 1
        r = dg_ac3(100_000)
        @test r.G ≈ ((3 + sqrt(5)) / 2)^2 rtol = 1e-12
        @test r.viol == 0
        @test r.pairgap == 0
        @test abs(r.lam[1] - log((3 + sqrt(5)) / 2)) <= 1e-3   # sanity of the reference estimator (finite-time term ~ 1/n)
    end

    @testset "AC4 pairing index mapping" begin   # pairs (i, 2N+1-i); recovered d1, d2 within eps(T) x scale; adjacent pairing gives ~1.27
        for T in dg_TYPES
            r = dg_ac4(T)
            @test abs(r.pair1 - r.d) <= eps(T) * (abs(r.L) + r.d)
            @test abs(r.sum1 - r.d) <= eps(T) * (abs(r.L) + r.d)
            tol = eps(T) * 2 * maximum(abs, r.lam2)
            @test abs(r.pair2 - r.d2) <= tol
            @test abs(r.sum2 - (r.d1 + r.d2)) <= tol
            @test r.pair2 < 0.5
        end
    end

    @testset "AC5 injected faults" begin
        @testset "det" begin   # det = 1 + delta exactly; raise for delta in [2 theta, 4 theta), quiet for (theta/4, theta/2]
            dg_load_faults()
            for T in dg_TYPES
                r = Base.invokelatest(dg_ac5a, T)
                @test !isempty(r.above) && !isempty(r.below)
                for (delta, e) in r.above
                    @test e isa NumericalHealthError
                    @test e.what == "det" && e.value == delta && e.orbit == 7 && e.step == 1
                end
                for (_, e) in r.below
                    @test e === nothing
                end
            end
        end

        @testset "orthogonality" begin   # eps_orth = 2 delta + delta^2 for one vector scaled by (1 + delta), absolute run tolerance
            for (T, theta_o) in ((Float64, 2.0^-30), (Double64, 2.0^-60))
                dg_load_faults()
                r = Base.invokelatest(dg_ac5b, T, theta_o)
                @test !isempty(r.above) && !isempty(r.below)
                for (delta, eo, e) in r.above
                    @test abs(eo - (2 * delta + delta^2)) <= eps(T)
                    @test e isa NumericalHealthError
                    @test e.what == "orthogonality" && abs(e.value - (2 * delta + delta^2)) <= eps(T) && e.step == 5
                end
                for (delta, eo, e) in r.below
                    @test abs(eo - (2 * delta + delta^2)) <= eps(T)
                    @test e === nothing
                end
            end
        end

        @testset "sum and pairing" begin   # delta on lambda_1 (N = 1): eps_sum = eps_pair = delta; two-sided against theta_sum
            for T in dg_TYPES
                r = dg_ac5c(T)
                tol = eps(T) * (abs(r.L) + 1)
                @test !isempty(r.above) && !isempty(r.below)
                for c in r.above
                    @test abs(c.es - c.delta) <= tol && abs(c.ep - c.delta) <= tol
                    @test c.e_mon isa NumericalHealthError && c.e_mon.what == "sum"   # sum precedes pairing in the fixed order
                    @test c.e_mon.step == 3 && c.e_mon.orbit == 9 && abs(c.e_mon.value - c.delta) <= tol
                    @test c.e_pair isa NumericalHealthError && c.e_pair.what == "pairing" && abs(c.e_pair.value - c.delta) <= tol
                end
                for c in r.below
                    @test abs(c.es - c.delta) <= tol && abs(c.ep - c.delta) <= tol
                    @test c.e_mon === nothing && c.e_pair === nothing
                end
            end
        end

        @testset "symplecticity" begin   # one 2 x 2 block with det 1 + delta; threshold_symp = tol eps N^2 (1 + M)^2, N = 2
            dg_load_faults()
            for T in dg_TYPES
                r = Base.invokelatest(dg_ac5d, T)
                tol = 16 * eps(T)
                @test !isempty(r.above) && !isempty(r.below)
                for c in r.above
                    @test abs(c.es - c.delta) <= tol
                    @test c.err isa NumericalHealthError && c.err.what == "symplecticity" && abs(c.err.value - c.delta) <= tol
                end
                for c in r.below
                    @test abs(c.es - c.delta) <= tol
                    @test c.err === nothing
                end
                # diag(2, 1/2, 1, 1) has det 1 but D' Omega D - Omega has entry 1: the monitor names symplecticity
                @test r.e_mon isa NumericalHealthError && r.e_mon.what == "symplecticity" && abs(r.e_mon.value - 1) <= tol
            end
        end

        @testset "non-finite" begin   # NaN or Inf raises with the right name even for an infinite threshold, never clipped
            dg_load_faults()
            for T in dg_TYPES
                r = Base.invokelatest(dg_ac5e, T)
                for (_, what, e) in r.cases
                    @test e isa NumericalHealthError && e.what == what && e.threshold == Inf
                end
                @test !any(r.never_clipped)
            end
        end
    end

    @testset "AC6 standard-map calibration" begin   # B1: max eps_det / (eps M) <= 1.5; 50 starts x 10^4 steps per K, panel K from benchmarks.md golden-mean row (0.9716)
        for (T, offset) in ((Float64, 600), (Double64, 700)), (i, K) in enumerate((0.5, 0.9716, 2.0, 5.0, 10.0))
            r = dg_ac6(T, K, StableRNG(dg_SEED + offset + i), 50, 10_000)
            @test r.ratio <= 1.5
            @test r.viol == 0
        end
    end

    @testset "AC8 256-bit recomputation" begin   # same stored entries, final subtraction rounding only: |diff| <= u |det| < eps
        @test dg_ac8() <= eps(Float64)
    end

    @testset "AC9 zero allocation" begin   # @allocated == 0 after warm-up, binary64 and double-double
        for T in dg_TYPES
            r = dg_ac9(T)
            @test r.orbit == 0
            @test r.update == 0
            @test r.exponents == 0
        end
    end

    @testset "AC10 thread invariance" begin   # propagate == serial map over propagate_orbit, bitwise; lowest faulty orbit id reported
        r = dg_ac10_series()
        @test r.same
        @test r.nontrivial
        dg_load_faults()
        one = Base.invokelatest(dg_ac10_fault, (37.0, NaN))
        @test one.e_par isa NumericalHealthError && one.e_par.orbit == 37 && one.e_par.what == "det"
        @test one.e_ser isa NumericalHealthError && one.e_ser.orbit == 37
        @test one.e_par.step == one.e_ser.step && one.e_par.value == one.e_ser.value
        two = Base.invokelatest(dg_ac10_fault, (41.0, 17.0))
        @test two.e_par isa NumericalHealthError && two.e_par.orbit == 17
        @test two.e_ser isa NumericalHealthError && two.e_ser.orbit == 17
    end

    @testset "AC11 CPU backend items 1, 5, 6" begin   # kernel (backend :cpu) gives the verdicts of AC1, AC5 (det) and AC6; GPU part is stage 04
        dg_load_faults()
        for T in dg_TYPES, c in Base.invokelatest(dg_ac11_item1, T)
            @test c.fin.det_max == 0 && c.fin.symp_max == 0
        end
        r5 = Base.invokelatest(dg_ac11_item5)
        @test !isempty(r5.above) && !isempty(r5.below)
        for (delta, e) in r5.above
            @test e isa NumericalHealthError && e.what == "det" && e.value == delta && e.orbit == 1
        end
        for (_, e) in r5.below
            @test e === nothing
        end
        for c in Base.invokelatest(dg_ac11_item6)
            @test c.err === nothing
            @test c.dmax <= 1.5 * eps(Float64) * (1 + c.K)   # B1 with M <= 1 + K
        end
    end

end
