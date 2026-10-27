# Windowed Fourier frequency of one window of a complex series (spec naff-fma: E1, NF2, NF16 to NF20, Amendments D80 and D81): the
# cosine window family, the coarse maximum of |phi| on the grid nu_m = -1/2 + m/(pad N), the stationary point of |phi|^2 in the bracket.

const NF_TIE_TOL = 64 * eps(Float64)         # D81 item 7: ties of the coarse search, for every element type
const NF_CANDIDATES = 64                     # local maxima of the transform whose neighbours are evaluated directly
const NF_CAND_REL = 1.5e-2                   # transform padding of at least 8 N loses under 0.5 percent of a Hann peak (derived), so 3x margin

_nf_pi(::Type{T}) where {T} = Systems._pi(one(T))

function _check_order_positions(order, positions)
    (order isa Integer && 0 <= order <= 5) || throw(ArgumentError("order must be an integer in 0:5, got $order"))
    positions in (:half, :edge) || throw(ArgumentError("positions must be :half or :edge, got $(repr(positions))"))
    return Int(order)
end

"""
    window_weights(T, N; order = 1, positions = :half) -> Vector{T}

Weights `chi_p(t_j) = (4^p (p!)^2 / (2p)!) sin(y_j)^(2p)` of the cosine window of order `p = order` (Laskar2003 Eq. (57), NF20; `p = 1` is
the Hanning window) in the Laskar normalisation, so the sum is `N`. `positions = :half` samples at `(j + 1/2)/N`, `:edge` at `j/N`
(`j = 0..N-1`). The index is mirrored (`min(j, N - 1 - j)` or `min(j, N - j)`) and the sine form is used, so the weights are bitwise
symmetric and keep full relative accuracy next to the ends (a `1 + cos` form would lose about `N^2 eps` there, D81 item 6).
`ArgumentError` naming `order` (outside 0:5) or `positions`.
"""
function window_weights(::Type{T}, N::Integer; order::Integer = 1, positions::Symbol = :half) where {T <: AbstractFloat}
    p = _check_order_positions(order, positions)
    N >= 1 || throw(ArgumentError("window length must be positive, got $N"))
    half = positions === :half
    c = T(4^p * factorial(p)^2) / T(factorial(2p))
    pit = _nf_pi(T)
    w = Vector{T}(undef, N)
    for j in 0:(N - 1)
        jj = half ? min(j, N - 1 - j) : min(j, N - j)
        y = pit * (T(jj) + (half ? one(T) / 2 : zero(T))) / N
        w[j + 1] = c * sin(y)^(2p)
    end
    return w
end

"""
    WindowFrequency{T}

Result of [`window_frequency`](@ref): `nu` (cycles per sample, in `[-1/2, 1/2)`), `amplitude = |phi(nu)| / sum(chi)`, `mirror_amplitude`
(`|phi(-nu)| / sum(chi)`, 0 when `-nu` is within `2(p+1)/N` of `nu` on the circle), `second_ratio` (largest local maximum of the transform
outside `nu` +- `2(p+1)/N` over the peak), `bin` (0-based index `m` of the coarse maximum), `tie` (another grid value within 64 eps of Float64
of the largest), `near_zero` (`abs(nu) < (p+1)/N`), `mirror`, `refined` (false if the stationary point was not found in the bracket and the
coarse value is kept) and `status` (`:ok`, `:degenerate`, `:failed`; `nu`, `amplitude`, `mirror_amplitude`, `second_ratio` are `NaN` unless `:ok`).
"""
struct WindowFrequency{T}
    nu::T
    amplitude::T
    mirror_amplitude::T
    second_ratio::T
    bin::Int
    tie::Bool
    near_zero::Bool
    mirror::Bool
    refined::Bool
    status::Symbol
end

_undefined_window(::Type{T}, status::Symbol) where {T} = WindowFrequency{T}(T(NaN), T(NaN), T(NaN), T(NaN), 0, false, false, false, false, status)

# Fractional part of nu * j, reduced to about [-1/2, 1/2], from the exact product of the operands (no rounding of nu * j for large j).
_frac(nu::AbstractFloat, j::Integer) = (x = nu * j; x - round(x))
function _frac(nu::Float64, j::Integer)
    x = nu * j
    return (x - round(x)) + fma(nu, Float64(j), -x)
end
function _frac(nu::Double64, j::Integer)
    jf = Float64(j)
    p = nu.hi * jf
    q = nu.lo * jf
    return ((Double64(p - round(p)) + fma(nu.hi, jf, -p)) + q) + fma(nu.lo, jf, -q)
end

# S_m = sum_j chi_j z_j (-2 pi i j)^m exp(-2 pi i nu j), m = 0, 1, 2, with Neumaier sums (D81 item 9).
function _sums(z::AbstractVector{Complex{T}}, chi::AbstractVector{T}, nu::T) where {T}
    tp = 2 * _nf_pi(T)
    a0r = a0i = a1r = a1i = a2r = a2i = CompensatedSum{T}()
    for j in 0:(length(z) - 1)
        s, c = sincos(tp * _frac(nu, j))
        t = chi[j + 1] * z[j + 1] * Complex(c, -s)
        d = Complex(zero(T), -tp * j)
        t1 = t * d
        t2 = t1 * d
        a0r, a0i = add(a0r, real(t)), add(a0i, imag(t))
        a1r, a1i = add(a1r, real(t1)), add(a1i, imag(t1))
        a2r, a2i = add(a2r, real(t2)), add(a2i, imag(t2))
    end
    return Complex(total(a0r), total(a0i)), Complex(total(a1r), total(a1i)), Complex(total(a2r), total(a2i))
end

# |phi| of the Float64 copies at nu, unnormalised: blocks of 32 samples, each restarted from the exactly reduced phase, the sum Neumaier-compensated
# (a plain sum has a rounding noise above the tie tolerance from N of about 8000, so the tie of exactly symmetric spectra would be lost).
function _abs_phi64(zf::Vector{ComplexF64}, chif::Vector{Float64}, nu::Float64)
    tp, n = 2 * pi, length(zf)
    w = cis(-tp * nu)
    sr = si = CompensatedSum{Float64}()
    j = 0
    while j < n
        jend = min(j + 31, n - 1)
        e = cis(-tp * _frac(nu, j))
        for jj in j:jend
            t = chif[jj + 1] * zf[jj + 1] * e
            sr, si = add(sr, real(t)), add(si, imag(t))
            e *= w
        end
        j = jend + 1
    end
    return hypot(total(sr), total(si))
end

# Coarse maximum of |phi| on the grid of M = pad N points: the transform of the zero-padded window (length the next power of two of at least
# max(M, 8 N), so the bin spacing is at most 1/M and the scalloping loss of a peak is small) names the local maxima near the top, whose neighbouring grid points are evaluated by direct sums; the largest
# decides, the lowest index among the values within NF_TIE_TOL wins and sets the tie flag (D81 item 7). Returns (bin, tie, spectrum), or `nothing`
# if the transform is not finite.
function _coarse(zf::Vector{ComplexF64}, chif::Vector{Float64}, pad::Int)
    N = length(zf)
    M = pad * N
    L = nextpow(2, max(M, 8N))
    a = zeros(ComplexF64, L)
    for j in 1:N
        a[j] = chif[j] * zf[j]
    end
    spec = abs.(_fft!(a))
    top = maximum(spec)
    (isfinite(top) && top > 0) || return nothing
    cand = [i for i in 1:L if spec[i] >= top * (1 - NF_CAND_REL) && spec[i] >= spec[mod1(i - 1, L)] && spec[i] >= spec[mod1(i + 1, L)]]
    sort!(cand; by = i -> -spec[i])
    grid = Int[]
    vals = Float64[]
    for i in first(cand, NF_CANDIDATES)
        nuk = (i - 1) / L
        nuk >= 0.5 && (nuk -= 1)
        m0 = round(Int, (nuk + 0.5) * M)
        for m in (m0 - 2):(m0 + 2)
            mm = mod(m, M)
            mm in grid && continue
            push!(grid, mm)
            push!(vals, _abs_phi64(zf, chif, -0.5 + mm / M))
        end
    end
    best = maximum(vals)
    near = [k for k in eachindex(vals) if vals[k] >= best * (1 - NF_TIE_TOL)]
    return minimum(grid[near]), length(near) > 1, spec
end

# Largest local maximum of the transform outside the circle distance `excl` of nu (NF2 step 5: a measure of a second strong line).
function _second_peak(spec::Vector{Float64}, nu::Float64, excl::Float64)
    L = length(spec)
    sec = 0.0
    for i in 1:L
        s = spec[i]
        (s > sec && s >= spec[mod1(i - 1, L)] && s >= spec[mod1(i + 1, L)]) || continue
        nui = (i - 1) / L
        nui >= 0.5 && (nui -= 1)
        d = abs(nui - nu)
        abs(d - round(d)) <= excl && continue
        sec = s
    end
    return sec
end

# Stationary point of |phi|^2 (Newton on its derivative, in the element type) in the bracket `nuc` +- h of the spec (NF2); iterates may wander to
# twice the half-width on the way, the converged point must lie inside it. (nuc, false) if the second derivative is not negative, the iterate
# leaves twice the bracket, the point is outside the bracket or it does not converge.
function _stationary(z, chi, nuc::T, h::T) where {T}
    nu = nuc
    extra = -1
    for _ in 1:20
        S0, S1, S2 = _sums(z, chi, nu)
        fp = 2 * real(conj(S0) * S1)
        fpp = 2 * (abs2(S1) + real(conj(S0) * S2))
        fpp < 0 || return nuc, false
        dx = fp / fpp
        nu -= dx
        abs(nu - nuc) > 2h && return nuc, false
        extra = extra < 0 ? (abs(dx) < 1e-9 ? 2 : -1) : extra - 1
        extra == 0 && return (abs(nu - nuc) <= h ? (nu, true) : (nuc, false))
    end
    return nuc, false
end

function _demean(z::AbstractVector{Complex{T}}) where {T}
    sr = si = CompensatedSum{T}()
    for x in z
        sr, si = add(sr, real(x)), add(si, imag(x))
    end
    m = Complex(total(sr) / length(z), total(si) / length(z))
    return Complex{T}[x - m for x in z]
end

function _sum_weights(chi::Vector{T}) where {T}
    s = CompensatedSum{T}()
    for w in chi
        s = add(s, w)
    end
    return total(s)
end

_reduce_circle(nu::T) where {T} = nu >= one(T) / 2 ? nu - 1 : (nu < -one(T) / 2 ? nu + 1 : nu)

# Result for a clean, non-degenerate window `z` (any mean already removed).
function _window_result(z::Vector{Complex{T}}, p::Int, positions::Symbol, pad::Int, mirror_ratio::Real) where {T}
    N = length(z)
    M = pad * N
    chi = window_weights(T, N; order = p, positions = positions)
    sumchi = _sum_weights(chi)
    zf = ComplexF64[ComplexF64(Float64(real(x)), Float64(imag(x))) for x in z]
    coarse = _coarse(zf, Float64.(chi), pad)
    coarse === nothing && return _undefined_window(T, :failed)
    bin, tie, spec = coarse
    nuc = -T(1) / 2 + T(bin) / T(M)
    nu, refined = _stationary(z, chi, nuc, one(T) / M)
    nu = _reduce_circle(nu)
    amp = abs(first(_sums(z, chi, nu))) / sumchi
    excl = 2 * (p + 1) / N
    dcirc = abs(2 * nu - round(2 * nu))
    mamp = dcirc < excl ? zero(T) : abs(first(_sums(z, chi, -nu))) / sumchi
    ratio = T(_second_peak(spec, Float64(nu), excl) / (Float64(amp) * Float64(sumchi)))
    return WindowFrequency{T}(nu, amp, mamp, ratio, bin, tie, abs(nu) < (p + 1) / N, mamp > T(mirror_ratio) * amp, refined, :ok)
end

# The window scaled by an exact power of two so that its largest component lies in [1, 2) (frequency, ties and flags are scale-free; the
# squares of the derivative sums cannot overflow or underflow); returns the scaled copy and the exponent `e` (the amplitudes scale by 2^e).
function _rescale(z::AbstractVector{Complex{T}}) where {T}
    m = zero(T)
    for x in z
        m = max(m, abs(real(x)), abs(imag(x)))
    end
    e = iszero(m) ? 0 : exponent(Float64(m))
    e == 0 && return Vector{Complex{T}}(z), 0
    return Complex{T}[Complex(ldexp(real(x), -e), ldexp(imag(x), -e)) for x in z], e
end

function _unscale(w::WindowFrequency{T}, e::Int) where {T}
    return WindowFrequency{T}(w.nu, ldexp(w.amplitude, e), ldexp(w.mirror_amplitude, e), w.second_ratio, w.bin, w.tie, w.near_zero, w.mirror, w.refined, w.status)
end

"""
    window_frequency(z; order = 1, positions = :half, pad = 8, mean_removal = false, mirror_ratio = 0.9) -> WindowFrequency

Fundamental frequency of one window `z::AbstractVector{Complex{T}}` of `N` samples (spec NF2, NF17): optional subtraction of the window
mean (compensated), the cosine window of order `order` ([`window_weights`](@ref)), the coarse maximum of `|phi|` on the grid
`-1/2 + m/(pad N)` (a Float64 transform locates the candidates, direct sums decide; the lowest index within 64 eps of Float64 of the
largest wins and sets `tie`), and the stationary point of `|phi|^2` in `[nu_c - 1/M, nu_c + 1/M]` by Newton's method in the element
type with compensated sums and exactly reduced phases. A window that is constant after the optional mean removal and either has
mean removal on or the constant 0 is `:degenerate` (D81 item 8); a non-finite sample gives `:failed`. Nothing is thrown for the data.
`ArgumentError` naming `window` (fewer than `8 (order + 1)` samples), `order`, `positions`, `pad` (below 2), `mirror_ratio`, before any work.
The window is first scaled by an exact power of two, so samples of any finite magnitude are handled and nothing overflows; a non-finite transform
(not reachable for finite data) gives `:failed`.
"""
function window_frequency(z::AbstractVector{Complex{T}}; order::Integer = 1, positions::Symbol = :half, pad::Integer = 8,
                          mean_removal::Bool = false, mirror_ratio::Real = 0.9) where {T <: AbstractFloat}
    p = _check_order_positions(order, positions)
    N = length(z)
    N >= 8 * (p + 1) || throw(ArgumentError("window of $N samples is below the minimum window $(8 * (p + 1)) for order $p"))
    pad >= 2 || throw(ArgumentError("pad must be at least 2 (the grid must be finer than the bin), got $pad"))
    0 < mirror_ratio <= 1 || throw(ArgumentError("mirror_ratio must be in (0, 1], got $mirror_ratio"))
    all(isfinite, z) || return _undefined_window(T, :failed)
    zs, e = _rescale(z)
    zz = mean_removal ? _demean(zs) : zs
    (mean_removal ? all(==(zz[1]), zz) : all(iszero, zz)) && return _undefined_window(T, :degenerate)
    w = _window_result(zz, p, positions, Int(pad), mirror_ratio)
    return e == 0 ? w : _unscale(w, e)
end
