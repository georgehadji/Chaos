# Helpers for the naff-fma suite (prefix nf_). Definitions only: nothing here calls a window_weights, window_frequency, frequency_map, check_frequency_map, circle_distance, frequency_distance,
# frequency_index, birkhoff_weight, birkhoff_weights, phase_increments, weighted_phase_frequency, WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable or check_run (new method) name at load
# time, so a missing symbol errors the testsets that use it and nothing else.
# Spec: workspaces/module-build/stages/01-spec/output/naff-fma-spec.md (D80, with Amendment D80). API proposal: workspaces/module-build/stages/02-tests/output/naff-fma-test-plan.md.
# The reference code of this file (cosine windows, the sampled scalar product and its derivatives, the BigFloat stationary point, the coarse search, the weighted phase advance, the closed
# forms, the signals) is plain-array code written from the spec formulas (NF2, NF6, NF9 to NF20, Eq. (57)); it never calls a function of the module under test to produce an expected value.
# nf_chi / nf_S / nf_stat_ref / nf_coarse are the independent implementation of the Oracle section (a direct-integration reference of the windowed Fourier integral).
# Orbits come from Systems (WrappedStandardMap, LinearSymplecticMap: verified modules), the regular class from the lyap-benettin folder (L1 S7), included once behind a guard.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats
@isdefined(lb_regular) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_helpers.jl"))

const NF_SEED = 20261012                  # base seed, _config/numerics-defaults.toml key `seed`
const NF_BASE = NF_SEED + 38000           # this slug uses 20261012 + 38000 .. 38999; testset AC<n> uses 50 n + sub (sub < 50); D80a/b/c use n = 20, 21, 22
nf_seed(ac::Integer, sub::Integer = 0) = NF_BASE + 50 * ac + sub
const NF_TYPES = (Float64, Double64)      # Float32 is out (D24); BigFloat only for oracle runs (and the pure E2 function, criterion 8(c))
const NF_BITS = 256                       # numerics-defaults.toml key `bigfloat_bits`
# Package suite sizes (a few minutes, D38) or the sizes of the spec: verify/naff-fma.jl sets NF_FULL_RUN = true before including runtests.jl. Every reduction is marked "D38" at its use.
const NF_FULL = isdefined(Main, :NF_FULL_RUN) && Main.NF_FULL_RUN
const NF_CACHE = Dict{Any, Any}()
# Relative tolerance of a coarse-grid tie: 64 eps of Float64 for both element types. The coarse transform is a Float64 transform (it only has to find the bin), so it cannot resolve the spec's 8 eps_T in
# Double64; and in Float64 the two grid values of an exact real cosine at +-0.2 differ by 14 eps (3.1e-15) in a direct summation with exactly reduced phases (plan, Items for the owner).
const NF_TIE_TOL = 64 * eps(Float64)

nf_eps(::Type{T}) where {T} = Float64(effective_eps(T))
nf_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
nf_big(x::AbstractFloat) = BigFloat(x)
nf_big(x::Integer) = BigFloat(x)
nf_big(z::Complex) = complex(nf_big(real(z)), nf_big(imag(z)))
nf_dec(s::AbstractString) = setprecision(BigFloat, NF_BITS) do
    parse(BigFloat, s)                    # decimals are read as decimals
end
nf_round(::Type{Float64}, x::BigFloat) = Float64(x)
nf_round(::Type{BigFloat}, x::BigFloat) = x
function nf_round(::Type{Double64}, x::BigFloat)
    hi = Float64(x)
    return Double64(hi) + Double64(Float64(x - hi))
end
nf_round(::Type{T}, z::Complex{BigFloat}) where {T} = complex(nf_round(T, real(z)), nf_round(T, imag(z)))
nf_series(::Type{T}, zb::AbstractVector{<:Complex{BigFloat}}) where {T} = Complex{T}[nf_round(T, z) for z in zb]
nf_cis(x) = complex(cos(x), sin(x))
nf_pi() = BigFloat(pi)

# Reduction to [-1/2, 1/2), the circle distance of NF9 and the Euclidean vector distance, in whatever (Big)Float the arguments have.
nf_red(x) = x - floor(x + 1 // 2)
nf_cd(a, b) = (d = a - b; abs(d - round(d)))
function nf_vdist(nu1, nu2)                                    # NF9 over planes
    return sqrt(sum(abs2, [nf_cd(a, b) for (a, b) in zip(nu1, nu2)]))
end
# Absolute error of a value of the element type against a BigFloat reference, as a Float64.
nf_err(value, ref) = setprecision(BigFloat, NF_BITS) do
    Float64(abs(nf_big(value) - ref))
end
# Circle distance between a module frequency (any element type) and a BigFloat truth, as a Float64.
nf_cderr(value, ref) = setprecision(BigFloat, NF_BITS) do
    Float64(nf_cd(nf_big(value), ref))
end

# ------------------------------------------------------------------ cosine windows (NF20) and the closed form of their transform (Eq. (57))

# chi_p(t_j) of NF20 in BigFloat: (2^p (p!)^2 / (2p)!) (1 + cos(pi t))^p, t_j = -1 + (2 j + 1)/N (:half, positions (j + 1/2)/N) or -1 + 2 j/N (:edge, positions j/N), j = 0..N-1.
function nf_chi(p::Int, N::Int, pos::Symbol)
    return get!(NF_CACHE, (:chi, p, N, pos)) do
        setprecision(BigFloat, NF_BITS) do
            c = BigFloat(2)^p * BigFloat(factorial(big(p)))^2 / BigFloat(factorial(big(2p)))
            pb = nf_pi()
            [c * (1 + cos(pb * (pos === :half ? -1 + BigFloat(2j + 1) / N : -1 + BigFloat(2j) / N)))^p for j in 0:(N - 1)]
        end
    end
end
# The transform of chi_p in closed form, Eq. (57): (-1)^p pi^(2p) (p!)^2 sin x / (x (x^2 - pi^2) ... (x^2 - p^2 pi^2)); BigFloat.
function nf_phi_closed(p::Int, x::BigFloat)
    return setprecision(BigFloat, NF_BITS) do
        pb = nf_pi()
        den = x
        for k in 1:p
            den *= x^2 - (k * pb)^2
        end
        (-1)^p * pb^(2p) * BigFloat(factorial(big(p)))^2 * sin(x) / den
    end
end
# The sampled transform (1/N) sum_j chi_p(t_j) exp(i x t_j) of given weights (any iterable of numbers), t_j as in nf_chi.
function nf_phi_sampled(w, x::BigFloat, pos::Symbol)
    return setprecision(BigFloat, NF_BITS) do
        N = length(w)
        s = zero(BigFloat)
        for j in 0:(N - 1)
            t = pos === :half ? -1 + BigFloat(2j + 1) / N : -1 + BigFloat(2j) / N
            s += nf_big(w[j + 1]) * cos(x * t)             # the weights are even about the centre: the sine part cancels
        end
        s / N
    end
end

# ------------------------------------------------------------------ the sampled scalar product and its derivatives (NF2, NF17), direct summation in BigFloat

# S_m(nu) = sum_j chi_j z_j (-2 pi i j)^m exp(-2 pi i nu j), m = 0, 1, 2, for the BigFloat series zb and weights chi (0-based j = index - 1). The exponentials are powers of one BigFloat
# exponential (the error after 10^4 products is 1e-73, far below the BigFloat(256) epsilon scale used for the comparisons).
function nf_S(zb::AbstractVector{Complex{BigFloat}}, chi::AbstractVector{BigFloat}, nu::BigFloat)
    return setprecision(BigFloat, NF_BITS) do
        tp = 2 * nf_pi()
        w = nf_cis(-tp * nu)
        e = one(Complex{BigFloat})
        s0 = zero(Complex{BigFloat})
        s1, s2 = s0, s0
        for j in 0:(length(zb) - 1)
            c = chi[j + 1] * zb[j + 1] * e
            d = complex(zero(BigFloat), -tp * j)
            c1 = c * d
            s0 += c
            s1 += c1
            s2 += c1 * d
            e *= w
        end
        (s0, s1, s2)
    end
end
# f = |S_0|^2, f' = 2 Re(conj(S_0) S_1), f'' = 2 (|S_1|^2 + Re(conj(S_0) S_2)); the amplitude is |S_0| / sum chi (NF18).
nf_fp(S) = 2 * real(conj(S[1]) * S[2])
nf_fpp(S) = 2 * (abs2(S[2]) + real(conj(S[1]) * S[3]))
nf_amp(zb, chi, nu) = setprecision(BigFloat, NF_BITS) do
    abs(nf_S(zb, chi, nu)[1]) / sum(chi)
end

# |S_0| of Float64 data at a Float64 frequency (the coarse search of the reference), blocks of 64 terms with a fresh exponential each.
function nf_abs64(zf::Vector{ComplexF64}, chif::Vector{Float64}, nu::Float64)
    n = length(zf)
    tp = 2 * pi
    w = cis(-tp * nu)
    s = 0.0 + 0.0im
    j = 0
    while j < n
        jend = min(j + 63, n - 1)
        e = cis(-tp * mod(nu * j, 1.0))
        for jj in j:jend
            s += chif[jj + 1] * zf[jj + 1] * e
            e *= w
        end
        j = jend + 1
    end
    return abs(s)
end
# f' = 2 Re(conj(S_0) S_1) in Float64 (analytic, used for the sign scan only).
function nf_f1_64(zf::Vector{ComplexF64}, chif::Vector{Float64}, nu::Float64)
    n = length(zf)
    tp = 2 * pi
    w = cis(-tp * nu)
    s0 = 0.0 + 0.0im
    s1 = 0.0 + 0.0im
    j = 0
    while j < n
        jend = min(j + 63, n - 1)
        e = cis(-tp * mod(nu * j, 1.0))
        for jj in j:jend
            c = chif[jj + 1] * zf[jj + 1] * e
            s0 += c
            s1 += c * complex(0.0, -tp * jj)
            e *= w
        end
        j = jend + 1
    end
    return 2 * real(conj(s0) * s1)
end

# The coarse maximum of |phi| on the grid nu_m = -1/2 + m/M, M = pad N, m = 0..M-1: the 0-based index m of the largest value (the lowest index wins an exact tie). Direct evaluation for N <= 1024; above that a
# two-stage search (every pad/2-th node, then the nodes around the 8 largest local maxima), which finds the same global maximum.
function nf_coarse(zf::Vector{ComplexF64}, chif::Vector{Float64}, pad::Int)
    N = length(zf)
    M = pad * N
    grid(m) = -0.5 + mod(m, M) / M
    if N <= 1024
        best, bm = -1.0, 0
        for m in 0:(M - 1)
            v = nf_abs64(zf, chif, grid(m))
            if v > best
                best, bm = v, m
            end
        end
        return bm
    end
    s = max(pad ÷ 2, 1)
    nodes = collect(0:s:(M - 1))
    vals = [nf_abs64(zf, chif, grid(m)) for m in nodes]
    L = length(nodes)
    loc = [i for i in 1:L if vals[i] >= vals[mod1(i - 1, L)] && vals[i] >= vals[mod1(i + 1, L)]]
    sort!(loc; by = i -> -vals[i])
    best, bm = -1.0, 0
    for i in loc[1:min(8, length(loc))], m in (nodes[i] - s):(nodes[i] + s)
        v = nf_abs64(zf, chif, grid(m))
        if v > best
            best, bm = v, mod(m, M)
        end
    end
    return bm
end

# The reference answer for one window of a T series: own coarse index, the BigFloat stationary point of |phi|^2 in [nu_c - 1/M, nu_c + 1/M] started from the coarse index `bin` (the module's, or the
# reference's own), whether it is unique with a negative second derivative in the bracket (a sign scan of f' at 33 nodes in Float64, then Newton in BigFloat to 1e-60), the amplitude there.
function nf_stat_ref(zT::AbstractVector{<:Complex}, p::Int, pos::Symbol, pad::Int; bin::Union{Nothing, Int} = nothing)
    N = length(zT)
    M = pad * N
    chi = nf_chi(p, N, pos)
    zb = Complex{BigFloat}[nf_big(z) for z in zT]
    zf = ComplexF64.(Float64.(real.(zb)), Float64.(imag.(zb)))
    chif = Float64.(chi)
    refbin = nf_coarse(zf, chif, pad)
    use = bin === nothing ? refbin : bin
    return setprecision(BigFloat, NF_BITS) do
        nuc = -BigFloat(1) / 2 + BigFloat(use) / M
        h = BigFloat(1) / M
        nodes = [Float64(nuc - h + 2h * i / 32) for i in 0:32]
        g = [nf_f1_64(zf, chif, nu) for nu in nodes]
        changes = [(g[i] > 0 && g[i + 1] < 0) ? 1 : ((g[i] < 0 && g[i + 1] > 0) ? -1 : 0) for i in 1:32]
        nchg = count(!=(0), changes)
        unique_max = nchg == 1 && sum(changes) == 1
        nu = nothing
        amp = nothing
        if unique_max
            i = findfirst(==(1), changes)
            x = BigFloat(nodes[i]) + (BigFloat(nodes[i + 1]) - BigFloat(nodes[i])) * g[i] / (g[i] - g[i + 1])
            for _ in 1:12
                S = nf_S(zb, chi, x)
                dx = nf_fp(S) / nf_fpp(S)
                x -= dx
                abs(dx) < BigFloat(10)^-60 && break
            end
            S = nf_S(zb, chi, x)
            unique_max = nf_fpp(S) < 0 && abs(x - nuc) <= 2h
            if unique_max
                nu = x
                amp = abs(S[1]) / sum(chi)
            end
        end
        (refbin = refbin, nu = nu, amplitude = amp, unique = unique_max, sign_changes = nchg)
    end
end

# Compare one module window result `w` (fields nu, bin, amplitude) with the reference of nf_stat_ref on the series zT. Returns (status, err) with status in
# :ok, :tie (the module bin differs from the reference bin but its value is not below the reference value by more than NF_TIE_TOL relative), :nonunique (the reference has no unique stationary point:
# counted, not failed), :binfail, :nufail.
function nf_compare_window(::Type{T}, zT, w, p::Int, pos::Symbol, pad::Int; tol::Real = 16) where {T}
    r = nf_stat_ref(zT, p, pos, pad; bin = Int(w.bin))
    epsT = nf_eps(T)
    N = length(zT)
    M = pad * N
    tie = false
    if Int(w.bin) != r.refbin
        chi = nf_chi(p, N, pos)
        zb = Complex{BigFloat}[nf_big(z) for z in zT]
        tie = setprecision(BigFloat, NF_BITS) do
            vm = nf_amp(zb, chi, -BigFloat(1) / 2 + BigFloat(Int(w.bin)) / M)
            vr = nf_amp(zb, chi, -BigFloat(1) / 2 + BigFloat(r.refbin) / M)
            vm >= vr * (1 - BigFloat(NF_TIE_TOL))
        end
        tie || return (:binfail, NaN)
    end
    r.unique || return (:nonunique, NaN)
    e = nf_cderr(w.nu, r.nu)
    a = setprecision(BigFloat, NF_BITS) do
        Float64(abs(nf_big(w.amplitude) - r.amplitude) / r.amplitude)
    end
    return (e <= tol * epsT && a <= 64 * epsT ? (tie ? :tie : :ok) : :nufail, e)
end

# ------------------------------------------------------------------ signals (BigFloat, generated by exact recurrences, rounded once)

# z_j = a exp(i (2 pi nu j + phi0)), j = 0..n-1.
function nf_line_big(nu::BigFloat, a::BigFloat, phi0::BigFloat, n::Int)
    return setprecision(BigFloat, NF_BITS) do
        w = nf_cis(2 * nf_pi() * nu)
        z = a * nf_cis(phi0)
        out = Vector{Complex{BigFloat}}(undef, n)
        for j in 1:n
            out[j] = z
            z *= w
        end
        out
    end
end
# Sum of lines with their own frequencies, amplitudes and phases: specs = [(nu, a, phi0), ...].
nf_lines_big(specs, n::Int) = setprecision(BigFloat, NF_BITS) do
    sum(nf_line_big(nu, a, phi0, n) for (nu, a, phi0) in specs)
end
# z_j = exp(2 pi i (nu0 j + alpha j^2 / 2)), j = 0..n-1 (NF11): the increments of the phase are 2 pi (nu0 + alpha (j + 1/2)).
function nf_chirp_big(nu0::BigFloat, alpha::BigFloat, n::Int)
    return setprecision(BigFloat, NF_BITS) do
        tp = 2 * nf_pi()
        c = nf_cis(tp * (nu0 + alpha / 2))
        r = nf_cis(tp * alpha)
        z = one(Complex{BigFloat})
        out = Vector{Complex{BigFloat}}(undef, n)
        for j in 1:n
            out[j] = z
            z *= c
            c *= r
        end
        out
    end
end
# Window means removed in exact arithmetic (per window of length N); the caller rounds the result to T.
function nf_demean_big(zb::AbstractVector{Complex{BigFloat}}, N::Int)
    return setprecision(BigFloat, NF_BITS) do
        out = copy(zb)
        for r in (1:N, (N + 1):(2N))
            m = sum(zb[r]) / N
            out[r] .= zb[r] .- m
        end
        out
    end
end

# ------------------------------------------------------------------ weighted Birkhoff weights and the weighted phase advance (NF4, NF6), BigFloat

nf_wb(t::BigFloat, p::Int = 1) = (t <= 0 || t >= 1) ? zero(BigFloat) : exp(-1 / (t^p * (1 - t)^p))
# The (unnormalised) weights w(j/N), j = 0..N-1, and the normalised ones.
nf_wb_raw(N::Int, p::Int = 1) = setprecision(BigFloat, NF_BITS) do
    [nf_wb(BigFloat(j) / N, p) for j in 0:(N - 1)]
end
nf_wb_hat(N::Int, p::Int = 1) = (w = nf_wb_raw(N, p); setprecision(BigFloat, NF_BITS) do
    w ./ sum(w)
end)
# NF6 for increments (radians, BigFloat) of 2N entries: (nu1, nu2) = sum_j w_hat_j Delta_j / (2 pi) per window, the second window restarting the weights.
function nf_E2_ref(incr::AbstractVector{BigFloat}, N::Int; p::Int = 1, weights = nothing, bits::Int = NF_BITS)
    return setprecision(BigFloat, bits) do
        what = weights
        if what === nothing
            w = [nf_wb(BigFloat(j) / N, p) for j in 0:(N - 1)]
            what = w ./ sum(w)
        end
        tp = 2 * nf_pi()
        (sum(what[j] * incr[j] for j in 1:N) / tp, sum(what[j] * incr[N + j] for j in 1:N) / tp)
    end
end
# Increments (BigFloat, principal argument of z_{j+1}/z_j) of a complex series given in T.
nf_incr_planar(zT::AbstractVector{<:Complex}) = setprecision(BigFloat, NF_BITS) do
    zb = [nf_big(z) for z in zT]
    [(r = zb[j + 1] * conj(zb[j]); atan(imag(r), real(r))) for j in 1:(length(zb) - 1)]
end
# The bound of criterion 10(d): B_E2 = 8 eps (1 + max_j |Delta_j| / (2 pi)).
nf_BE2(::Type{T}, incr) where {T} = 8 * nf_eps(T) * (1 + Float64(maximum(abs, incr)) / (2 * pi))

# NF15 (Das2017 Example 1): rho = sqrt(2) - 1, observed lifted angle in cycles phi_n = n rho + alpha cos(2 pi theta_n) + beta sin(2 pi theta_n), theta_n = n rho mod 1, alpha = 0.1, beta = 0.2.
# Returns the BigFloat increments in RADIANS, 2 pi (phi_{n+1} - phi_n), n = 0..nn-1, and rho.
function nf_nf15_big(nn::Int)
    return get!(NF_CACHE, (:nf15, nn)) do
        setprecision(BigFloat, NF_BITS) do
            rho = sqrt(BigFloat(2)) - 1
            tp = 2 * nf_pi()
            a, b = nf_dec("0.1"), nf_dec("0.2")
            phi(n) = (th = n * rho - floor(n * rho); n * rho + a * cos(tp * th) + b * sin(tp * th))
            ph = [phi(n) for n in 0:nn]
            ([tp * (ph[n + 1] - ph[n]) for n in 1:nn], rho)
        end
    end
end

# ------------------------------------------------------------------ orbits (Systems: verified modules) and the observables

nf_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) :
    label === :five ? exact(T, 5) : label === :zero ? exact(T, 0) : error("unknown K $label")
nf_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]
nf_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))

# The states of an orbit at the handoffs 1..n (tau = 1): the samples that the series feed and the stream reducers see. Wrapped states are (x, p, m, n).
function nf_states(sys, z0, n::Integer)
    out = Vector{typeof(z0)}(undef, n)
    z = z0
    for i in 1:n
        z = forward(sys, z)
        out[i] = z
    end
    return out
end
# Samples at the handoffs tau, 2 tau, ..., n tau of the same orbit (tau map steps per sample).
function nf_states_tau(sys, z0, n::Integer, tau::Integer)
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
nf_obs_cis(x) = nf_cis(x[1])                         # exp(i q) from the wrapped coordinate (P4)
nf_obs_planar(x) = complex(x[1], x[2])               # q + i p
nf_obs_plane(k, d) = x -> complex(x[k], x[k + d])    # plane k of a direct sum of blocks ordered (q_1..q_N, p_1..p_N)
nf_obs(f, states) = [f(s) for s in states]
# The lifted increments of the wrapped coordinate q from the states: (x' - x) + tp (m' - m) with tp = lift of one winding (the period of the wrapped map in its own arithmetic).
function nf_lifted_incr(states::AbstractVector{<:SVector{4, T}}) where {T}
    tp = lift(SVector{4, T}(0, 0, 1, 0))[1]
    return T[(states[j + 1][1] - states[j][1]) + tp * (states[j + 1][3] - states[j][3]) for j in 1:(length(states) - 1)]
end
# The same increments in BigFloat from the exactly widened states (period tp as realized in T).
nf_lifted_incr_big(states::AbstractVector{<:SVector{4, T}}) where {T} = setprecision(BigFloat, NF_BITS) do
    tp = nf_big(lift(SVector{4, T}(0, 0, 1, 0))[1])
    [(nf_big(states[j + 1][1]) - nf_big(states[j][1])) + tp * (nf_big(states[j + 1][3]) - nf_big(states[j][3])) for j in 1:(length(states) - 1)]
end

# Y5 rotation blocks: the map of a direct sum of blocks (angles in thetas), ordered (q_1..q_N, p_1..p_N); a start is uniform in [-1, 1]^(2N).
nf_rot_map(::Type{T}, thetas...) where {T} = LinearSymplecticMap(length(thetas) == 1 ? rotation_block(T, thetas[1]) : block_diagonal(map(th -> rotation_block(T, th), thetas)...))
nf_rot_start(::Type{T}, rng, d::Int) where {T} = SVector{d, T}(ntuple(_ -> T(2 * rand(rng) - 1), d))
# Frequency of a plane with block angle theta: -theta / (2 pi) reduced (NF14), BigFloat from the Float64 value of theta.
nf_rot_nu(theta::Real) = setprecision(BigFloat, NF_BITS) do
    nf_red(-BigFloat(theta) / (2 * nf_pi()))
end

# Regular-class and panel starts of the standard map, widened exactly into T with their wrapped state.
function nf_case(label::Symbol, ::Type{T}, m::Int, seedsub::Int) where {T}
    return get!(NF_CACHE, (:case, label, T, m, seedsub)) do
        sys = WrappedStandardMap(nf_K(T, label))
        if label === :half
            r = lb_regular()
            ids, starts, fseed = r.ids[1:m], r.starts[1:m], r.seed
        else
            fseed = nf_seed(0, seedsub)
            ids = collect(1:m)
            starts = nf_starts(StableRNG(fseed), m)
        end
        (sys = sys, x0 = [torus_state(nf_cast(T, s)) for s in starts], ids = ids, fseed = fseed, label = label)
    end
end

# ------------------------------------------------------------------ the stream (reducer protocol) and a reader recording the orbit

# Read-only reducer recording the orbit point at every handoff (CPU only: Vector state).
struct nf_Trace end
Chaoslab.reducer_init(::nf_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::nf_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::nf_Trace, st) = (xs = st,)

# One orbit through the per-orbit body with a chain; returns (x, W, finals) with finals the tuple of reducer_final(r, state, ws). The workspace is allocated here when a reducer asks for one.
function nf_run(rs::Tuple, sys, x0, W0, nsteps::Integer, tau::Integer; orbit::Integer = 1)
    T = eltype(x0)
    d = size(W0, 1)
    k = size(W0, 2)
    wslen = sum(r -> Chaoslab.reducer_workspace(r, T, Val(d), Val(k)), rs; init = 0)
    ws = zeros(UInt32, wslen)
    x, W, st = propagate_orbit(sys, x0, W0, nsteps, tau, rs; orbit = orbit, workspace = ws)
    return (x, W, map((r, s) -> Chaoslab.reducer_final(r, s, ws), rs, st))
end
nf_frame(seed, j, ::Type{T}, d::Int = 2, k::Int = 2) where {T} = initial_deviations(seed, j, Val(d), Val(k), T)

# ------------------------------------------------------------------ small utilities

function nf_error(f)
    try
        f()
        return nothing
    catch e
        return e
    end
end
function nf_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end
nf_mean(v) = sum(v) / length(v)
# Bitwise equality of nested results (NaN equals NaN, -0.0 differs from 0.0).
nf_same(a, b) = isequal(a, b)
# The fields of a window result that make up its identity (for bitwise comparisons of whole results).
nf_wkey(w) = (w.nu, w.amplitude, w.mirror_amplitude, w.second_ratio, w.bin, w.tie, w.near_zero, w.mirror, w.refined, w.status)

# A reducer that only normalises the columns (the replacer of the chains whose only other members are readers, so that a long run of an expanding map does not overflow).
struct nf_Norm end
Chaoslab.reducer_init(::nf_Norm, ::Type, ::Val, ::Val) = 0
function Chaoslab.reducer_update(::nf_Norm, st, n::Int, x, W::SMatrix{d, k}) where {d, k}
    return st, hcat(ntuple(j -> W[:, j] / safe_norm(W[:, j]), Val(k))...)
end
Chaoslab.reducer_final(::nf_Norm, st) = (;)

# Failure bookkeeping: a list of short tags, printed (first five) when a testset fails; true when empty.
nf_report(bad, name) = (isempty(bad) || @info "$name: $(length(bad)) failures, first five" first = first(bad, 5); isempty(bad))

# A vector of the concrete state type from a Vector{Any} of states (the trace reader's state).
nf_vec(v) = convert(Vector{typeof(first(v))}, v)
