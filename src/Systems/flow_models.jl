# The Hamiltonian flows of the integrator-catalog tests (catalog systems.md, section B): F1 harmonic oscillator, F4 pendulum, F2 planar Kepler,
# T2 Henon-Heiles (test-only here, Hairer2008 Sect. 2) and the nonseparable T1 (Tao2016 Sec. II A). Parameters are built in the state's element type.

@inline _sq(x) = x * x

"""
    HarmonicOscillator(w)

F1: `H = sum_j (p_j^2 + w_j^2 q_j^2)/2` with frequencies `w` (a vector in the element type); exact flow a rotation in each plane.
"""
struct HarmonicOscillator{N, T} <: AbstractSeparableFlow
    w::SVector{N, T}
end
HarmonicOscillator(w::AbstractVector{T}) where {T} = HarmonicOscillator{length(w), T}(SVector{length(w), T}(w))

"F4: `H = p^2/2 - eps cos q`, the pendulum (`eps = 1` the standard one, small `eps` the near-integrable form of LaskarRobutel2001)."
struct Pendulum{T} <: AbstractSeparableFlow
    eps::T
end

"F2: planar Kepler problem `H = |p|^2/2 - 1/|q|` (GM = 1), state (q1, q2, p1, p2)."
struct Kepler <: AbstractSeparableFlow end

"T2: Henon-Heiles `H = |p|^2/2 + |q|^2/2 + q1^2 q2 - q2^3/3`, state (q1, q2, p1, p2)."
struct HenonHeilesFlow <: AbstractSeparableFlow end

"T1: the nonseparable `H = (Q^2 + 1)(P^2 + 1)/2` of Tao2016 Sec. II A, state (Q, P)."
struct TaoNonseparable <: AbstractFlow end

dimension(::HarmonicOscillator{N}) where {N} = 2N
dimension(::Pendulum) = 2
dimension(::Kepler) = 4
dimension(::HenonHeilesFlow) = 4
dimension(::TaoNonseparable) = 2

# ---- F1
function hamiltonian(f::HarmonicOscillator{N}, z::SVector) where {N}
    return sum(ntuple(j -> (_sq(z[N + j]) + _sq(f.w[j]) * _sq(z[j])) / 2, Val(N)))
end
force(f::HarmonicOscillator, q::SVector) = (f.w .* f.w) .* q
force_jacobian(f::HarmonicOscillator{N, T}, ::SVector) where {N, T} = SMatrix{N, N, T}(ntuple(l -> (l - 1) % N == (l - 1) ÷ N ? f.w[(l - 1) % N + 1]^2 : zero(T), Val(N * N)))
corrector_force(f::HarmonicOscillator, z::SVector) = 2 .* (f.w .^ 4) .* _q(z)
corrector_hessian(f::HarmonicOscillator{N, T}, ::SVector) where {N, T} = SMatrix{N, N, T}(ntuple(l -> (l - 1) % N == (l - 1) ÷ N ? 2 * f.w[(l - 1) % N + 1]^4 : zero(T), Val(N * N)))

function exact_flow(f::HarmonicOscillator{N}, z::SVector{D}, t) where {N, D}
    return SVector{D}(ntuple(i -> begin
        j = i <= N ? i : i - N
        w = f.w[j]; c = cos(w * t); s = sin(w * t)
        i <= N ? z[j] * c + z[N + j] / w * s : -z[j] * w * s + z[N + j] * c
    end, Val(D)))
end
function exact_flow_matrix(f::HarmonicOscillator{N, T}, t) where {N, T}
    S = promote_type(T, typeof(t))
    entry(i, j) = begin
        a = i <= N ? i : i - N; b = j <= N ? j : j - N
        a == b || return zero(S)
        w = f.w[a]; c = cos(w * t); s = sin(w * t)
        i <= N ? (j <= N ? c : s / w) : (j <= N ? -w * s : c)
    end
    return SMatrix{2N, 2N, S}(ntuple(l -> entry((l - 1) % (2N) + 1, (l - 1) ÷ (2N) + 1), Val(4N * N)))
end

# ---- F4
hamiltonian(f::Pendulum, z::SVector) = _sq(z[2]) / 2 - f.eps * _cos(z[1])
force(f::Pendulum, q::SVector) = SVector(f.eps * _sin(q[1]))
force_jacobian(f::Pendulum, q::SVector) = SMatrix{1, 1}(f.eps * _cos(q[1]))
perturbation_scale(f::Pendulum) = f.eps
corrector_force(::Pendulum, z::SVector) = SVector(_sin(2 * z[1]))
corrector_hessian(::Pendulum, z::SVector) = SMatrix{1, 1}(2 * _cos(2 * z[1]))

# ---- F2
hamiltonian(::Kepler, z::SVector) = (_sq(z[3]) + _sq(z[4])) / 2 - inv(sqrt(_sq(z[1]) + _sq(z[2])))
function force(::Kepler, q::SVector)
    r2 = _sq(q[1]) + _sq(q[2])
    return q / (r2 * sqrt(r2))
end
function force_jacobian(::Kepler, q::SVector{2, T}) where {T}
    r2 = _sq(q[1]) + _sq(q[2]); r3 = r2 * sqrt(r2); r5 = r3 * r2
    return SMatrix{2, 2, T}(inv(r3) - 3 * _sq(q[1]) / r5, -3 * q[1] * q[2] / r5, -3 * q[1] * q[2] / r5, inv(r3) - 3 * _sq(q[2]) / r5)
end
function corrector_force(::Kepler, z::SVector)   # W = r^-4
    q = _q(z); r2 = _sq(q[1]) + _sq(q[2])
    return -4 * q / (r2 * r2 * r2)
end
function corrector_hessian(::Kepler, z::SVector{4, T}) where {T}
    q = _q(z); r2 = _sq(q[1]) + _sq(q[2]); r6 = r2 * r2 * r2; r8 = r6 * r2
    return SMatrix{2, 2, T}(-4 / r6 + 24 * _sq(q[1]) / r8, 24 * q[1] * q[2] / r8, 24 * q[1] * q[2] / r8, -4 / r6 + 24 * _sq(q[2]) / r8)
end
angular_momentum(::Kepler, z::SVector) = z[1] * z[4] - z[2] * z[3]
function runge_lenz(f::Kepler, z::SVector)
    r = sqrt(_sq(z[1]) + _sq(z[2])); L = angular_momentum(f, z)
    return SVector(z[4] * L - z[1] / r, -z[3] * L - z[2] / r)
end

# ---- T2
hamiltonian(::HenonHeilesFlow, z::SVector) = (_sq(z[3]) + _sq(z[4])) / 2 + (_sq(z[1]) + _sq(z[2])) / 2 + _sq(z[1]) * z[2] - z[2]^3 / 3
force(::HenonHeilesFlow, q::SVector) = SVector(q[1] + 2 * q[1] * q[2], q[2] + _sq(q[1]) - _sq(q[2]))
force_jacobian(::HenonHeilesFlow, q::SVector{2, T}) where {T} = SMatrix{2, 2, T}(1 + 2 * q[2], 2 * q[1], 2 * q[1], 1 - 2 * q[2])
corrector_force(f::HenonHeilesFlow, z::SVector) = 2 * (force_jacobian(f, _q(z)) * force(f, _q(z)))
function corrector_hessian(f::HenonHeilesFlow, z::SVector{4, T}) where {T}   # W = |g|^2: 2 (H H + sum_k g_k dH/dq_k); third derivatives of V: V112 = 2, V222 = -2
    q = _q(z); g = force(f, q); Hs = force_jacobian(f, q)
    T1 = SMatrix{2, 2, T}(0, 2, 2, 0); T2 = SMatrix{2, 2, T}(2, 0, 0, -2)   # dH/dq1, dH/dq2
    return 2 * (Hs * Hs + g[1] * T1 + g[2] * T2)
end

# ---- T1 (nonseparable)
hamiltonian(::TaoNonseparable, z::SVector) = (_sq(z[1]) + 1) * (_sq(z[2]) + 1) / 2
hamiltonian_gradient(::TaoNonseparable, z::SVector) = SVector(z[1] * (_sq(z[2]) + 1), z[2] * (_sq(z[1]) + 1))
hamiltonian_hessian(::TaoNonseparable, z::SVector{2, T}) where {T} = SMatrix{2, 2, T}(_sq(z[2]) + 1, 2 * z[1] * z[2], 2 * z[1] * z[2], _sq(z[1]) + 1)
vector_field(f::TaoNonseparable, z::SVector) = (g = hamiltonian_gradient(f, z); SVector(g[2], -g[1]))
function vector_field_jacobian(::TaoNonseparable, z::SVector{2, T}) where {T}
    return SMatrix{2, 2, T}(2 * z[1] * z[2], -(_sq(z[2]) + 1), _sq(z[1]) + 1, -2 * z[1] * z[2])
end
