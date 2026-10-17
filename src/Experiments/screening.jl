# Screening (E11): L5 at Float64 with k = 1 and the first deviation vector (cos alpha, sin alpha). With one vector the
# Householder reflection gives R_11 = ||W||, so the loop carries the vector, takes its norm at each handoff and keeps the
# compensated cumulative sum S(n) of ln R_11 (Skokos2010 Eq. (82), p.47). The hot loop allocates nothing.

"Buffers of one screening run: the block exponents and the three values the predicates read."
mutable struct ScreenBuffer
    blocks::Vector{Float64}
    lambda::Float64
    S_tenth::Float64
    S_N::Float64
end

"`ScreenBuffer` for `N` steps in blocks of `w` steps."
screen_buffer(N::Integer, w::Integer) = ScreenBuffer(zeros(N ÷ w), 0.0, 0.0, 0.0)

"Bytes held per grid cell during screening: 2 state and 4 tangent numbers (48), a compensated accumulator (16), `N / w` block exponents."
screen_bytes(N::Integer, w::Integer) = 48 + 16 + 8 * (N ÷ w)

"""
    screen_point!(buf, K, x, p, N, w, alpha, tau) -> buf

Run `N` steps of the lifted standard map from `(x, p)` with handoff interval `tau` and fill `buf`: `lambda = S(N) / N`,
`S_tenth = S(N / 10)` (at the last handoff not beyond it), `S_N = S(N)` and the block exponents `(S((i + 1) w) - S(i w)) / w`.
`N` must be a multiple of `w` and `w` of `tau`.
"""
function screen_point!(buf::ScreenBuffer, K::Float64, x::Float64, p::Float64, N::Integer, w::Integer, alpha::Float64, tau::Integer)
    sys = StandardMap(K)
    z = SVector(x, p)
    v = SVector(cos(alpha), sin(alpha))
    acc = CompensatedSum{Float64}()
    h10, hw = (N ÷ 10) ÷ tau, w ÷ tau
    prev, s10, blk = 0.0, 0.0, 1
    for h in 1:(N ÷ tau)
        for _ in 1:tau
            v = tangent(sys, z) * v
            z = forward(sys, z)
        end
        r = safe_norm(v)
        acc = add(acc, log(r))
        v = v / r
        h == h10 && (s10 = total(acc))
        if h % hw == 0
            s = total(acc)
            buf.blocks[blk] = (s - prev) / w
            prev, blk = s, blk + 1
        end
    end
    buf.S_N = total(acc)
    buf.S_tenth = s10
    buf.lambda = buf.S_N / N
    return buf
end

"""
    screen_point(K, x, p; N_scr, w, alpha, tau = 1) -> (lambda, blocks, S_tenth, S_N)

Allocating form of [`screen_point!`](@ref).
"""
function screen_point(K::Real, x::Real, p::Real; N_scr::Integer, w::Integer, alpha::Real, tau::Integer = 1)
    buf = screen_buffer(N_scr, w)
    screen_point!(buf, Float64(K), Float64(x), Float64(p), N_scr, w, Float64(alpha), tau)
    return (lambda = buf.lambda, blocks = buf.blocks, S_tenth = buf.S_tenth, S_N = buf.S_N)
end
