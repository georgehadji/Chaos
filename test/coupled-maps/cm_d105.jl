# Amendment D105 (rework of the chain state): the Y4 systems, lifted and wrapped, accept a fixed-size state (SVector of length 2N lifted, 4N wrapped),
# so that Y4 runs through the tangent stream like Y3. Assumed new names (test plan, "Rework D105"): the banded operator type `ChainTangent` (exported from
# Systems), `forward`, `inverse`, `tangent`, `inverse_tangent` methods on SVector, `torus_state(sys, ::SVector)`, `lift(sys, ::SVector)` for the wrapped chain.
# D105d (GPU against CPU for Y4) is stage 04 and has no test here. Reference values are the vector-state results and the banded frame actions of the
# existing suite (`tangent_action!`, `inverse_tangent_action!`), compared bitwise (isequal).

const CM105_KE = ((3.0, 0.2), (3.0, 0.6102), (1.0, -0.5), (0.5, 1.0), (5.0, 0.6102), (0.0, 1.0))   # (K, eps) of the Y4 panel, as CM_Y4_CONFIGS

cm105_sv(v::AbstractVector) = SVector{length(v)}(v)
cm105_sm(W::AbstractMatrix) = SMatrix{size(W, 1), size(W, 2)}(W)

# wrapped state: reduced coordinates in [0, 2 pi), small integer counts (layout q, p, q-counts, p-counts), element type T
function cm105_wstate(::Type{T}, rng, d::Integer) where {T}
    return T[i <= d ? T(2 * Float64(pi) * rand(rng)) : T(rand(rng, -3:3)) for i in 1:(2d)]
end

cm105_frame(::Type{T}, rng, d::Integer) where {T} = reshape(T[T(v) for v in cm_draw(rng, d * 4, 1.0)], d, 4)

# Number of failed comparisons of the fixed-size results against the vector-state results, for one case, lifted and wrapped, at one point each.
function cm105_mismatches(c::CMCase, ::Type{T}, rng) where {T}
    d = cm_d(c)
    bad = 0
    W = cm105_frame(T, rng, d)
    Ws = cm105_sm(W)
    out = zeros(T, d, 4)
    sys, wsys = cm_sys(T, c), cm_wsys(T, c)
    for (sy, chain, v) in ((sys, sys, cm_state(T, c, cm_draw(rng, d, 4 * Float64(pi)))), (wsys, wsys.chain, cm105_wstate(T, rng, d)))
        s = cm105_sv(v)
        for f in (forward, inverse)
            a = f(sy, s)
            bad += !(a isa SVector && length(a) == length(v) && isequal(Vector(a), f(sy, v)))
        end
        for (tn, act!) in ((tangent, tangent_action!), (inverse_tangent, inverse_tangent_action!))
            M = tn(sy, s)
            bad += !(M isa ChainTangent && M isa AbstractMatrix && size(M) == (d, d) && isequal(Matrix(M), tn(sy, v)))
            P = M * Ws
            act!(out, chain, v, W)
            bad += !(P isa SMatrix && size(P) == (d, 4) && isequal(Matrix(P), out))
        end
    end
    zl = [T(x) for x in cm_draw(rng, d, 4 * Float64(pi))]
    t = torus_state(wsys, cm105_sv(zl))
    bad += !(t isa SVector && isequal(Vector(t), torus_state(wsys, zl)))
    sw = cm105_wstate(T, rng, d)
    l = lift(wsys, cm105_sv(sw))
    bad += !(l isa SVector && isequal(Vector(l), lift(wsys, sw)))
    return bad
end

@testset "D105a fixed-size state equals vector state" begin
    # N in {2, 3, 5, 10, 100}, open and ring (ring needs N >= 3), Float64 and Double64: forward, inverse, entries of tangent and inverse tangent (Matrix of the
    # ChainTangent, bitwise), the frame products ChainTangent * SMatrix (bitwise the banded actions), torus_state and lift of the wrapped chain, lifted and wrapped.
    # N <= 10: the six (K, eps) of the Y4 panel, 2 seeded points each; N = 100 (compile time): two parameter sets, 1 point each. One testset per (T, N, boundary).
    rng = StableRNG(CM_SEED + 105001)
    for T in (Float64, Double64), N in (2, 3, 5, 10, 100), ring in (false, true)
        ring && N < 3 && continue
        @testset "T=$T N=$N ring=$ring" begin
            ke = N == 100 ? CM105_KE[2:3] : CM105_KE
            bad = 0
            for (K, e) in ke, _ in 1:(N == 100 ? 1 : 2)
                bad += cm105_mismatches(cm_y4(N, K, e; ring), T, rng)
            end
            @test bad == 0
        end
    end
end

# plain serial loop of the stream body: tangent at the old point, then forward (src/Tangent/orbit.jl)
function cm105_serial(sys, x, W, n::Integer)
    for _ in 1:n
        W = tangent(sys, x) * W
        x = forward(sys, x)
    end
    return x, W
end

@testset "D105b tangent stream equals plain loop" begin
    # propagate on the CPU backend, Y4 N in {2, 5} at (K, eps) = (3, 0.2), lifted and wrapped, 2000 seeded orbits, k = 4, n = n_H steps with one handoff (tau = n, no reducers),
    # explicit seeded initial frames W0: final states and deviation frames bitwise equal to the plain serial loop of forward and ChainTangent products
    rng = StableRNG(CM_SEED + 105002)
    for N in (2, 5), wrapped in (false, true)
        @testset "N=$N wrapped=$wrapped" begin
            c = cm_y4(N, 3.0, 0.2)
            sys = wrapped ? cm_wsys(Float64, c) : cm_sys(Float64, c)
            n, d = cm_nH(c), 2N
            X0 = [cm105_sv(wrapped ? cm105_wstate(Float64, rng, d) : cm_draw(rng, d, 4 * Float64(pi))) for _ in 1:2000]
            W0 = [cm105_sm(cm105_frame(Float64, rng, d)) for _ in 1:2000]
            res = propagate(sys, X0, n; k = 4, tau = n, seed = CM_SEED, W0)
            ref = [cm105_serial(sys, X0[j], W0[j], n) for j in eachindex(X0)]
            # the loop itself must stay in fixed-size types (without this the dense fallback also matches bitwise)
            @test tangent(sys, X0[1]) isa ChainTangent
            @test all(r -> r[1] isa SVector && r[2] isa SMatrix, ref)
            @test all(j -> isequal(res.x[j], ref[j][1]), eachindex(X0))
            @test all(j -> isequal(res.W[j], ref[j][2]), eachindex(X0))
        end
    end
end

# bytes allocated by one step and one frame product, each, after a warm-up call
function cm105_alloc(sys, s, W)
    forward(sys, s); inverse(sys, s); tangent(sys, s) * W; inverse_tangent(sys, s) * W
    return (@allocated(forward(sys, s)), @allocated(inverse(sys, s)), @allocated(tangent(sys, s) * W), @allocated(inverse_tangent(sys, s) * W))
end

@testset "D105c plain data and allocation" begin
    # isbits system values (lifted and wrapped, open and ring, N up to 100, both types); 0 bytes for forward, inverse and the frame products of tangent and
    # inverse tangent at a fixed-size state, N in {2, 10, 100}, lifted and wrapped, Float64 and Double64
    rng = StableRNG(CM_SEED + 105003)
    for T in (Float64, Double64), N in (2, 3, 10, 100), ring in (false, true)
        ring && N < 3 && continue
        c = cm_y4(N, 3.0, 0.2; ring)
        @test isbitstype(typeof(cm_sys(T, c)))
        @test isbitstype(typeof(cm_wsys(T, c)))
    end
    for T in (Float64, Double64), N in (2, 10, 100)
        c = cm_y4(N, 3.0, 0.2)
        d = 2N
        Ws = cm105_sm(cm105_frame(T, rng, d))
        @test all(==(0), cm105_alloc(cm_sys(T, c), cm105_sv(cm_state(T, c, cm_draw(rng, d, 4.0))), Ws))
        @test all(==(0), cm105_alloc(cm_wsys(T, c), cm105_sv(cm105_wstate(T, rng, d)), Ws))
    end
end
