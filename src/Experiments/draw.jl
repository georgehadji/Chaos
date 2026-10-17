# Streams, seeds and draws (E12 to E16, E19b).

"""
    cell_stream_seed(run_seed, i, class) -> UInt64

Seed of the stream of panel cell (K index `i`, `class`) (E15): the first 8 bytes, big-endian, of SHA-256 over the UTF-8 text
`chaoslab-panel-1|<run_seed>|<i>|<class>`. Independent of threads and devices.
"""
function cell_stream_seed(run_seed::Integer, i::Integer, class)
    h = sha256("chaoslab-panel-1|$(run_seed)|$(i)|$(class)")
    return foldl((acc, b) -> (acc << 8) | UInt64(b), h[1:8]; init = UInt64(0))
end

"Generator of panel cell (`i`, `class`): a `StableRNG` keyed by [`cell_stream_seed`](@ref)."
cell_rng(run_seed::Integer, i::Integer, class) = StableRNG(cell_stream_seed(run_seed, i, class))

"Per-orbit seed of orbit `id` in run `run_seed` (E16): tangent-core's rule, the key of the generator behind `initial_deviations`."
orbit_seed(run_seed::Integer, id::Integer) = Tangent._mix(Tangent._mix(run_seed % UInt64) + id % UInt64)

"""
    draw_candidates(rng, cands, n_c) -> Vector{Int}

Partial Fisher-Yates draw (E14) of `min(n_c, M)` distinct entries of the ascending candidate list `cands`, one uniform number per draw.
"""
function draw_candidates(rng, cands, n_c::Integer)
    c = collect(Int, cands)
    M = length(c)
    out = Int[]
    for t in 1:min(n_c, M)
        j = t + floor(Int, rand(rng) * (M - t + 1))
        c[t], c[j] = c[j], c[t]
        push!(out, c[t])
    end
    return out
end

"X start (E12): `r = 10^l`, `l` uniform on `[lo, hi]`, angle `phi = 2 pi u` from the hyperbolic fixed point (0, 0). Returns `(x, p, r, phi)`."
function place_X(rng, lo, hi)
    u1, u2 = rand(rng), rand(rng)
    r = 10.0^(lo + (hi - lo) * u1)
    phi = 2pi * u2
    return (x = r * cos(phi), p = r * sin(phi), r = r, phi = phi)
end

"""
    place_R(rng, pts, delta_R)

R start at K > 0 (E13): a uniformly chosen point of `pts` displaced by `delta_R` in a uniform direction, the sum formed at the
ambient BigFloat precision and rounded once. Returns `(index, phi, x, p, xb, pb)` with the unrounded sum in `xb`, `pb`.
"""
function place_R(rng, pts, delta_R)
    u1, u2 = rand(rng), rand(rng)
    index = 1 + floor(Int, u1 * length(pts))
    phi = 2pi * u2
    xb = pts[index][1] + BigFloat(delta_R) * cos(BigFloat(phi))
    pb = pts[index][2] + BigFloat(delta_R) * sin(BigFloat(phi))
    return (index = index, phi = phi, x = Float64(xb), p = Float64(pb), xb = xb, pb = pb)
end

"R start at K = 0 (E13): a reduced fraction `a / b`, `b <= q_max`, in the order `(b, a)`; `p = 2 pi a / b`, `x = 2 pi u`. Returns `(a, b, x, p)`."
function place_R0(rng, q_max::Integer)
    fr = [(a, b) for b in 1:q_max for a in 0:(b - 1) if gcd(a, b) == 1]
    u1, u2 = rand(rng), rand(rng)
    a, b = fr[1 + floor(Int, u1 * length(fr))]
    p = setprecision(() -> Float64(2 * BigFloat(pi) * a / b), BigFloat, 256)
    return (a = a, b = b, x = 2pi * u2, p = p)
end

"""
    y5_control(mu) -> (sys, x0, exponent)

Hyperbolic linear symplectic control of dimension 2 with eigenvalues `mu` and `1/mu` (the block of `standard-map-spec.md` item 10),
started at (0, 0); the exponent is `ln mu`.
"""
function y5_control(mu::Real)
    A, ex = setprecision(BigFloat, 256) do
        lm = log(BigFloat(mu))
        SMatrix{2, 2, Float64}(Float64.(hyperbolic_block(BigFloat, lm))), Float64(lm)
    end
    return (sys = LinearSymplecticMap(A), x0 = SVector(0.0, 0.0), exponent = ex)
end
