# Helpers for the weighted-birkhoff suite (prefix wb_). Definitions only: nothing here calls a weighted_birkhoff, birkhoff_index, WeightedBirkhoff, CosObservable, SinObservable, LiftedObservable,
# observable_name or check_run (new method) name at load time, so a missing symbol errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/weighted-birkhoff-spec.md (D83). API proposal: workspaces/module-build/stages/02-tests/output/weighted-birkhoff-test-plan.md.
# The reference code of this file (weights from the definition in 256 bits, the window sums of WB14, the index of WB15, the exact rotation signals, the Poisson form WB17, the BigFloat standard-map orbit)
# is plain-array code written from the spec formulas; it never calls a function of the module under test to produce an expected value (Oracle, package suite).
# Orbits come from Systems (WrappedStandardMap, LinearSymplecticMap: verified modules), the regular class from the lyap-benettin folder (L1 S7), included once behind a guard.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats
@isdefined(lb_regular) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_helpers.jl"))

const WB_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const WB_BASE = WB_SEED + 39000           # this slug uses 20261012 + 39000 .. 39999 = 20300012 ..; testset AC<n> uses 50 n + sub (sub < 50); records R<n> use n = 20 + n; the API uses n = 0
wb_seed(ac::Integer, sub::Integer = 0) = WB_BASE + 50 * ac + sub
const WB_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs and the pure function
const WB_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
const WB_MULT = 32                        # floor multiplier of WB15 (spec P2, proposal)
# Package suite sizes (D38 class, a few minutes) or the sizes of the spec: verify/weighted-birkhoff.jl sets WB_FULL_RUN = true before including runtests.jl. Every reduction is marked "D38" at its use.
const WB_FULL = isdefined(Main, :WB_FULL_RUN) && Main.WB_FULL_RUN
const WB_CACHE = Dict{Any, Any}()

wb_eps(::Type{T}) where {T} = Float64(effective_eps(T))
wb_cp(p::Integer) = 4.0^p                 # c_p of WB18 (proposal, P9)

# ------------------------------------------------------------------ BigFloat conversions

wb_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, WB_BITS)
wb_big(x::AbstractFloat) = BigFloat(x)
wb_big(x::Integer) = BigFloat(x)
wb_dec(s::AbstractString) = setprecision(BigFloat, WB_BITS) do
    parse(BigFloat, s)                    # decimals are read as decimals
end
wb_round(::Type{Float64}, x::BigFloat) = Float64(x)
wb_round(::Type{BigFloat}, x::BigFloat) = x
function wb_round(::Type{Double64}, x::BigFloat)
    hi = Float64(x)
    return Double64(hi) + Double64(Float64(x - hi))
end
wb_pi() = BigFloat(pi)
wb_err(value, ref) = setprecision(BigFloat, WB_BITS) do
    Float64(abs(wb_big(value) - ref))
end
# 2 pi_T as realized by lift (decision D28) in T, widened exactly.
wb_tp(::Type{T}) where {T} = setprecision(BigFloat, WB_BITS) do
    wb_big(lift(SVector{4, T}(0, 0, 1, 0))[1])
end

# ------------------------------------------------------------------ the reference: weights, window sums (WB14) and the index (WB15), plain arrays in BigFloat(256)

# w_jj = w^[p](min(jj, N - jj)/N), jj = 0..N-1, from the definition exp(-[t (1 - t)]^-p) in 256 bits (w_0 = 0).
function wb_wref(N::Int, p::Int)
    return get!(WB_CACHE, (:w, N, p)) do
        setprecision(BigFloat, WB_BITS) do
            [jj == 0 ? zero(BigFloat) : (t = BigFloat(min(jj, N - jj)) / N; exp(-1 / (t * (1 - t))^p)) for jj in 0:(N - 1)]
        end
    end
end
# The window averages of one series of BigFloat samples fb (at least 2N entries): (avg_1, avg_2, B, S) with S the largest |f| of the 2N samples; `w` overrides the weights (recorded lines).
function wb_ref_avgs(fb::AbstractVector{BigFloat}, N::Int, p::Int; w::Union{Nothing, AbstractVector} = nothing)
    return setprecision(BigFloat, WB_BITS) do
        ww = w === nothing ? wb_wref(N, p) : BigFloat[wb_big(x) for x in w]
        B = sum(ww)
        a1 = sum(ww[j] * fb[j] for j in 1:N) / B
        a2 = sum(ww[j] * fb[N + j] for j in 1:N) / B
        (avg1 = a1, avg2 = a2, B = B, S = maximum(abs, fb[1:(2N)]))
    end
end
# WB15 in BigFloat for epsT (a Float64, a power of two). reldig is NaN where undefined (m <= phi), the digits are NaN for S = 0.
function wb_ref_index(a1::BigFloat, a2::BigFloat, S::BigFloat, epsT::Real; mult::Real = WB_MULT)
    return setprecision(BigFloat, WB_BITS) do
        delta = a1 - a2
        phi = mult * BigFloat(epsT) * S
        S == 0 && return (delta = delta, phi = phi, absdig = BigFloat(NaN), reldig = BigFloat(NaN), dig = BigFloat(NaN), m = zero(BigFloat), below_floor = true, zero = true)
        m = (abs(a1) + abs(a2)) / 2
        mx = max(abs(delta), phi)
        absdig = -log10(mx)
        reldig = m > phi ? -log10(mx / m) : BigFloat(NaN)
        dig = isnan(reldig) ? absdig : max(absdig, reldig)
        (delta = delta, phi = phi, absdig = absdig, reldig = reldig, dig = dig, m = m, below_floor = abs(delta) <= phi, zero = false)
    end
end
# The module-side weight vector of WB14 in T (unnormalised, mirrored index, built from the existing verified birkhoff_weight): used as the "weight vector equal to the module's weights" of criterion 6.
wb_wvec(::Type{T}, N::Int, p::Int) where {T} = T[birkhoff_weight(T(min(jj, N - jj)) / T(N), p) for jj in 0:(N - 1)]

# ------------------------------------------------------------------ exact signals (BigFloat, generated by exact recurrences, rounded once)

# (cos, sin) of 2 pi k (theta0 + j alpha), j = 0..n-1, as BigFloat vectors.
function wb_rot_big(theta0::BigFloat, alpha::BigFloat, k::Int, n::Int)
    return setprecision(BigFloat, WB_BITS) do
        tp = 2 * wb_pi()
        z = complex(cos(tp * k * theta0), sin(tp * k * theta0))
        w = complex(cos(tp * k * alpha), sin(tp * k * alpha))
        c = Vector{BigFloat}(undef, n)
        s = Vector{BigFloat}(undef, n)
        for j in 1:n
            c[j], s[j] = real(z), imag(z)
            z *= w
        end
        (c, s)
    end
end
wb_round_vec(::Type{T}, v::AbstractVector{BigFloat}) where {T} = T[wb_round(T, x) for x in v]

# ------------------------------------------------------------------ the Poisson form WB17 (criterion 12(a)): W(s) = integral_0^1 w(t) exp(2 pi i s t) dt by the trapezoid rule on nq nodes
# (the aliasing error is W(n q +- s), below 1e-40 for nq = 1000 and |s| <= 500), the average of cos 2 pi (theta0 + j alpha) over the first window as Re[e^(2 pi i theta0) sum_m W(N (alpha - m))] / sum_m W(N m).
function wb_Wint(s::BigFloat, nq::Int)
    return setprecision(BigFloat, 128) do
        wt = get!(WB_CACHE, (:wq, nq)) do
            [(t = BigFloat(i) / nq; exp(-1 / (t * (1 - t)))) for i in 1:(nq - 1)]
        end
        tp = 2 * BigFloat(pi)
        re = zero(BigFloat)
        im_ = zero(BigFloat)
        for i in 1:(nq - 1)
            a = tp * s * i / nq
            re += wt[i] * cos(a)
            im_ += wt[i] * sin(a)
        end
        complex(re, im_) / nq
    end
end
function wb_poisson_avg(N::Int, theta0::BigFloat, alpha::BigFloat; nq::Int = 1000, smax::Int = 500)
    return setprecision(BigFloat, 128) do
        ms = -(smax ÷ N + 2):(smax ÷ N + 3)
        num = sum(wb_Wint(N * (alpha - m), nq) for m in ms if abs(N * (alpha - m)) <= smax)
        den = sum(wb_Wint(BigFloat(N * m), nq) for m in ms if abs(N * m) <= smax)
        tp = 2 * BigFloat(pi)
        real(complex(cos(tp * theta0), sin(tp * theta0)) * num) / real(den)
    end
end

# ------------------------------------------------------------------ orbits (Systems: verified modules)

wb_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) : label === :five ? exact(T, 5) :
    label === :nine ? exact(T, 0.9) : label === :zero ? exact(T, 0) : error("unknown K $label")
wb_Kf(label::Symbol) = label === :two_pi ? 2pi : label === :half ? 0.5 : label === :ten ? 10.0 : label === :five ? 5.0 : label === :nine ? 0.9 : label === :zero ? 0.0 : error("unknown K $label")
wb_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]
wb_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))

# The states of an orbit at the handoffs 1..n (tau = 1): the samples that the stream reducers see. Wrapped states are (x, p, m, n).
function wb_states(sys, z0, n::Integer)
    out = Vector{typeof(z0)}(undef, n)
    z = z0
    for i in 1:n
        z = forward(sys, z)
        out[i] = z
    end
    return out
end
# Samples at the handoffs tau, 2 tau, ..., n tau (tau map steps per sample).
function wb_states_tau(sys, z0, n::Integer, tau::Integer)
    out = Vector{typeof(z0)}(undef, n)
    z = z0
    for i in 1:n
        for _ in 1:tau
            z = forward(sys, z)
        end
        out[i] = z
    end
    return out
end

# The regular class of L1 S7 (lb_regular) extended past its 50 starts: the same candidate stream, the same screening (X(1e5) <= 5e-4, Float64, tau = 1, batches of 60), the first m regular candidates; the first 50 are those of lb_regular.
function wb_regular(m::Int)
    m <= 50 && return lb_regular()
    return get!(WB_CACHE, (:regular, m)) do
        seed = LB_SEED + LB_OFF + 70
        sys = StandardMap(0.5)
        rng = StableRNG(seed)
        cands = SVector{2, Float64}[]
        ok = Int[]
        while length(ok) < m
            lo = length(cands) + 1
            append!(cands, lb_starts(rng, 60))
            W0 = [initial_deviations(seed, c, Val(2), Val(1), Float64) for c in lo:length(cands)]
            res = lb_run(sys, cands[lo:end], 100_000; k = 1, tau = 1, seed = seed, W0 = W0)
            append!(ok, (lo - 1) .+ findall(<=(5.0e-4), [f[1].lambda[1] for f in res.finals]))
        end
        keep = ok[1:m]
        (ids = keep, starts = cands[keep], rejected = last(keep) - m, seed = seed)
    end
end

# Regular-class and panel starts of the standard map, widened exactly into T with their wrapped state. Labels: :half (the regular class of L1 S7), :two_pi, :five, :ten, :nine, :zero.
function wb_case(label::Symbol, ::Type{T}, m::Int, seedsub::Int) where {T}
    return get!(WB_CACHE, (:case, label, T, m, seedsub)) do
        sys = WrappedStandardMap(wb_K(T, label))
        if label === :half
            r = wb_regular(m)
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = wb_seed(0, seedsub)
            ids = collect(1:m)
            starts = wb_starts(StableRNG(fseed), m)
        end
        (sys = sys, x0 = [torus_state(wb_cast(T, s)) for s in starts], starts = starts, ids = ids, fseed = fseed, label = label)
    end
end

# The five observables of the standard map in BigFloat (columns: cos q, sin q, cos p, sin p, lifted p) at the realized wrapped states, widened exactly; the lifted p uses 2 pi as realized in T.
function wb_fref_std(sts::AbstractVector{<:SVector{4, T}}) where {T}
    return setprecision(BigFloat, WB_BITS) do
        tp = wb_tp(T)
        out = Matrix{BigFloat}(undef, length(sts), 5)
        for (i, s) in enumerate(sts)
            q, p = wb_big(s[1]), wb_big(s[2])
            out[i, 1], out[i, 2], out[i, 3], out[i, 4] = cos(q), sin(q), cos(p), sin(p)
            out[i, 5] = p + tp * wb_big(s[4])
        end
        out
    end
end
# The same five observables evaluated by the module's kinds in T (the stream reducer's input).
wb_series(obs::Tuple, sts) = begin
    T = eltype(first(sts))
    T[T(o(x)) for x in sts, o in obs]
end

# The independent orbit of criterion 12(b): the standard map on lifted coordinates in BigFloat(256) from the widened start, (q, p) after steps 1..n, K a Float64 (exactly representable).
function wb_orbit_big(Kf::Float64, z0::SVector{2, Float64}, n::Int)
    return setprecision(BigFloat, WB_BITS) do
        K = BigFloat(Kf)
        q, p = BigFloat(z0[1]), BigFloat(z0[2])
        out = Vector{Tuple{BigFloat, BigFloat}}(undef, n)
        for i in 1:n
            p = p + K * sin(q)
            q = q + p
            out[i] = (q, p)
        end
        out
    end
end

# Y5 rotation blocks (a linear symplectic map of the state (q, p)), uniform starts in [-1, 1]^2: used by the health tests.
wb_rot_map(::Type{T}, theta) where {T} = LinearSymplecticMap(rotation_block(T, theta))

# ------------------------------------------------------------------ the stream (reducer protocol) and a reader recording the orbit

struct wb_Trace end
Chaoslab.reducer_init(::wb_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::wb_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::wb_Trace, st) = (xs = st,)

# A reducer that only normalises the columns (the replacer of the chains whose other members are readers, so that a long run of an expanding map does not overflow: K_std = 2 pi gives a growth
# of about 2 per step and K_std = 10 of about 5, so a column overflows within 1000 steps in Float64).
struct wb_Norm end
Chaoslab.reducer_init(::wb_Norm, ::Type, ::Val, ::Val) = 0
function Chaoslab.reducer_update(::wb_Norm, st, n::Int, x, W::SMatrix{d, k}) where {d, k}
    return st, hcat(ntuple(j -> W[:, j] / safe_norm(W[:, j]), Val(k))...)
end
Chaoslab.reducer_final(::wb_Norm, st) = (;)

# One orbit through the per-orbit body with a chain; returns (x, W, finals) with finals the tuple of reducer_final(r, state, ws).
function wb_run(rs::Tuple, sys, x0, W0, nsteps::Integer, tau::Integer; orbit::Integer = 1)
    T = eltype(x0)
    d = size(W0, 1)
    k = size(W0, 2)
    wslen = sum(r -> Chaoslab.reducer_workspace(r, T, Val(d), Val(k)), rs; init = 0)
    ws = zeros(UInt32, wslen)
    x, W, st = propagate_orbit(sys, x0, W0, nsteps, tau, rs; orbit = orbit, workspace = ws)
    return (x, W, map((r, s) -> Chaoslab.reducer_final(r, s, ws), rs, st))
end
wb_frame(seed, j, ::Type{T}, d::Int = 2, k::Int = 1) where {T} = initial_deviations(seed, j, Val(d), Val(k), T)

# ------------------------------------------------------------------ callables for the tests (any isbits callable of the orbit point is an observable)

struct wb_Const{C}                        # the constant c in the element type of the orbit point
    c::C
end
(o::wb_Const)(x) = oftype(x[1], o.c)
struct wb_Aff{F}                          # a f(x) + b, a and b Float64 (exact powers of two in the covariance tests)
    f::F
    a::Float64
    b::Float64
end
(o::wb_Aff)(x) = oftype(x[1], o.a) * oftype(x[1], o.f(x)) + oftype(x[1], o.b)
struct wb_Coord                           # the coordinate itself: the wrapped angle (a discontinuous observable) or any component
    i::Int
end
(o::wb_Coord)(x) = x[o.i]
struct wb_StdKick                         # f(x) = lifted p + K sin q of the wrapped standard map (criterion 7): the lifted increment of q, a caller function
    K::Float64
end
(o::wb_StdKick)(x) = (x[2] + 2 * Chaoslab.Systems._pi(x[1]) * x[4]) + oftype(x[1], o.K) * sin(x[1])
struct wb_BadAt                           # NaN as soon as the first coordinate exceeds c (a finite orbit with a non-finite observable)
    c::Float64
end
(o::wb_BadAt)(x) = x[1] > o.c ? oftype(x[1], NaN) : x[1]
struct wb_Tuple                           # not a real number: a "value" that is not of the required kind
end
(o::wb_Tuple)(x) = (x[1], x[1])

# ------------------------------------------------------------------ small utilities

function wb_error(f)
    try
        f()
        return nothing
    catch e
        return e
    end
end
function wb_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end
# Bitwise equality of nested results (NaN equals NaN, -0.0 differs from 0.0).
wb_same(a, b) = isequal(a, b)
# Failure bookkeeping: a list of short tags, printed (first five) when a testset fails; true when empty.
wb_report(bad, name) = (isempty(bad) || @info "$name: $(length(bad)) failures, first five" first = first(bad, 5); isempty(bad))
wb_vec(v) = convert(Vector{typeof(first(v))}, v)
# Quantile of a vector (type 7).
wb_quant(v, q) = (s = sort(filter(!isnan, collect(v))); isempty(s) ? NaN : s[clamp(round(Int, 1 + q * (length(s) - 1)), 1, length(s))])

# ------------------------------------------------------------------ observable sets of the standard map (functions: they name the new kinds only when called)

# The five observables of the criteria: cos q, sin q, cos p, sin p of the wrapped coordinates and the lifted p (state (x, p, m, n) of the WrappedStandardMap).
wb_std5() = (CosObservable(1), SinObservable(1), CosObservable(2), SinObservable(2), LiftedObservable(2, 4))
# The four of criteria 6 and 8: cos q, sin q, lifted p, cos 2 q.
wb_std4() = (CosObservable(1), SinObservable(1), LiftedObservable(2, 4), CosObservable(1, 2))
# The index of every entry of a result with matrices avg1, avg2, S: matrices of the fields of birkhoff_index (a helper of the tests, called only inside testsets).
wb_idx(res) = (ix = birkhoff_index.(res.avg1, res.avg2, res.S); (delta = map(t -> t.delta, ix), absdig = map(t -> t.absdig, ix), reldig = map(t -> t.reldig, ix), dig = map(t -> t.dig, ix),
                                                                  below_floor = map(t -> t.below_floor, ix), zero = map(t -> t.zero, ix)))
