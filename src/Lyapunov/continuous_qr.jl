# Continuous QR method for Hamiltonian flows (catalog L8): the deviation frame Q obeys an ODE (the penalty form of the
# continuous Gram-Schmidt flow), integrated by Heun (order 2) or RK4 (order 4) on the orbit grid of the delivered one-step
# map; the exponents are the compensated time integrals of the diagonal growth rates (Q' A Q)_ii. Spec: continuous-qr C1 to C21.
# The module owns a per-orbit loop (not a reducer of the tangent stream); the orbit is advanced by `forward(map, z)` alone.

"""
    ContinuousQR(map, flow; order = 4, penalty_factor = nothing, orthogonality_tol = CQ_EMAX, form = :auto)

Continuous QR estimator. `map` is the one-step map of a symplectic integrator (`forward` and a `tau` field), `flow` the Hamiltonian
flow (`vector_field_jacobian`). `order` is 2 (Heun, one orbit step per Q-step) or 4 (RK4, two orbit steps per Q-step); `form` is `:auto`
(form 84 for `k = d`, form 92 for `k < d`), `:f84` (`k = d` only) or `:f92`; the penalty is `nu = c_nu ||A||_F` with `c_nu` defaulting to 1
for `k = d` and 2 for `k < d`. A Q-step whose `max |Q'Q - I|` exceeds `orthogonality_tol` throws `NumericalHealthError("orthogonality", ...)`.
"""
struct ContinuousQR{M, F}
    map::M
    flow::F
    order::Int
    penalty_factor::Union{Nothing, Float64}
    orthogonality_tol::Float64
    form::Symbol
end

"Default orthogonality guard `E_max` of [`ContinuousQR`](@ref) (spec P5 as amended at stage 02, test plan U6)."
const CQ_EMAX = 1.0e-2

function ContinuousQR(map, flow; order::Integer = 4, penalty_factor = nothing, orthogonality_tol::Real = CQ_EMAX, form::Symbol = :auto)
    order in (2, 4) || throw(ArgumentError("order must be 2 or 4, got $order"))
    form in (:auto, :f84, :f92) || throw(ArgumentError("form must be :auto, :f84 or :f92, got $form"))
    return ContinuousQR(map, flow, Int(order), penalty_factor === nothing ? nothing : Float64(penalty_factor), Float64(orthogonality_tol), form)
end

_cqr_form(cq::ContinuousQR, ::Val{d}, ::Val{k}) where {d, k} = cq.form === :auto ? (k == d ? :f84 : :f92) : cq.form
_cqr_cnu(cq::ContinuousQR, ::Val{d}, ::Val{k}) where {d, k} = cq.penalty_factor === nothing ? (k == d ? 1.0 : 2.0) : cq.penalty_factor

# Argument checks, all before any work.
function _check_cqr(cq::ContinuousQR, ::Val{d}, ::Val{k}, N::Integer) where {d, k}
    1 <= k <= d || throw(ArgumentError("k must be in 1:$d, got $k"))
    cq.form === :f84 && k != d && throw(ArgumentError("form :f84 needs k = d (got k = $k, d = $d)"))
    m = cq.order ÷ 2
    (N >= 1 && N % m == 0) || throw(ArgumentError("N must be a positive multiple of $m for order $(cq.order), got $N"))
    return nothing
end

"`cqr_growth(Q, A)`: the growth rates `u_i = (Q' A Q)_ii` of the columns of `Q`."
function cqr_growth(Q::SMatrix{d, k}, A) where {d, k}
    AQ = A * Q
    return SVector{k}(ntuple(i -> dot(Q[:, i], AQ[:, i]), Val(k)))
end

# S = strict lower triangle of U minus its transpose (skew generator that keeps Q upper-triangular in R).
_skew_lower(U::SMatrix{k, k, T}) where {k, T} =
    SMatrix{k, k, T}(ntuple(l -> (i = (l - 1) % k + 1; j = (l - 1) ÷ k + 1; i > j ? U[i, j] : i < j ? -U[j, i] : zero(T)), Val(k * k)))

"""
    cqr_rhs(Q, A, nu; form = :auto)

Right-hand side `G(Q, A)` of the penalty form of the continuous QR flow: form 84 (`k = d`) is `Q S - nu Q E`, form 92 (`k <= d`) is
`A Q - Q U + Q S (I + E) - nu Q E`, with `U = Q' A Q`, `S` the skew part built from the strict lower triangle of `U`, `E = Q'Q - I`. Pure.
"""
function cqr_rhs(Q::SMatrix{d, k, T}, A, nu; form::Symbol = :auto) where {d, k, T}
    form === :auto && (form = k == d ? :f84 : :f92)
    form === :f84 && k != d && throw(ArgumentError("form :f84 needs k = d"))
    AQ = A * Q
    U = Q' * AQ
    S = _skew_lower(U)
    E = Q' * Q - one(SMatrix{k, k, T})
    pen = nu * (Q * E)
    return form === :f84 ? Q * S - pen : AQ - Q * U + (Q * S) * (one(SMatrix{k, k, T}) + E) - pen
end

_emax(Q::SMatrix{d, k, T}) where {d, k, T} = maximum(abs, Q' * Q - one(SMatrix{k, k, T}))

# One stage: right-hand side and growth rates at (Q, A) with penalty nu = c_nu ||A||_F.
@inline _stage(Q::SMatrix, A, cnu, form::Symbol) = (cqr_rhs(Q, A, cnu * norm(A); form = form), cqr_growth(Q, A))

# Jacobian at an orbit point, with the health check of the deviation matrix (step convention of test plan U4).
@inline function _next_jacobian(cq::ContinuousQR, z, orbit, n, backend)
    A = vector_field_jacobian(cq.flow, z)
    all(isfinite, A) || throw(NumericalHealthError("deviation matrix", NaN, 0.0, n, orbit, backend))
    return A
end

# Q-step of order 2 (Heun) over one orbit step: 2 right-hand sides, 1 new Jacobian. Returns (z, A, Q, increment).
@inline function _step_heun(cq::ContinuousQR, z, Q, A, p, h, orbit, n, backend)
    z1 = forward(cq.map, z)
    A1 = _next_jacobian(cq, z1, orbit, n + 1, backend)
    K1, u1 = _stage(Q, A, p.cnu, p.form)
    K2, u2 = _stage(Q + h * K1, A1, p.cnu, p.form)
    return z1, A1, Q + (h / 2) * (K1 + K2), (h / 2) * (u1 + u2)
end

# Q-step of order 4 (RK4) over two orbit steps with nodes A, A1, A1, A2: 4 right-hand sides, 2 new Jacobians.
@inline function _step_rk4(cq::ContinuousQR, z, Q, A, p, h, orbit, n, backend)
    z1 = forward(cq.map, z)
    A1 = _next_jacobian(cq, z1, orbit, n + 1, backend)
    z2 = forward(cq.map, z1)
    A2 = _next_jacobian(cq, z2, orbit, n + 2, backend)
    D = 2h
    K1, u1 = _stage(Q, A, p.cnu, p.form)
    K2, u2 = _stage(Q + h * K1, A1, p.cnu, p.form)
    K3, u3 = _stage(Q + h * K2, A1, p.cnu, p.form)
    K4, u4 = _stage(Q + D * K3, A2, p.cnu, p.form)
    return z2, A2, Q + (D / 6) * (K1 + 2K2 + 2K3 + K4), (D / 6) * (u1 + 2u2 + 2u3 + u4)
end

@inline _addall(acc::NTuple{k}, inc) where {k} = ntuple(i -> add(acc[i], inc[i]), Val(k))
@inline _totals(acc::Tuple) = SVector(map(total, acc))

# Series sink: `nothing` records nothing (the zero-allocation path); a vector records (n, l, E, Q) every `every` Q-steps.
@inline _record!(::Nothing, every, m, n, acc, E, Q) = nothing
@inline function _record!(series::Vector, every, m, n, acc, E, Q)
    (n ÷ m) % every == 0 && push!(series, (n = n, l = _totals(acc), E = E, Q = Q))
    return nothing
end

# The run: N orbit steps in Q-steps of m = order / 2 orbit steps. Returns (z, Q, final); N + 1 Jacobian calls in total.
function _cqr_run(cq::ContinuousQR, x0::SVector, Q0::SMatrix{d, k, T}, N::Int, orbit::Int, backend::Symbol, series, every::Int) where {d, k, T}
    _check_cqr(cq, Val(d), Val(k), N)
    all(isfinite, x0) || throw(NumericalHealthError("state", NaN, 0.0, 0, orbit, backend))
    all(isfinite, Q0) || throw(NumericalHealthError("deviation matrix", NaN, 0.0, 0, orbit, backend))
    p = (cnu = T(_cqr_cnu(cq, Val(d), Val(k))), form = _cqr_form(cq, Val(d), Val(k)))
    h = T(cq.map.tau); m = cq.order ÷ 2
    z = x0; Q = Q0; A = _next_jacobian(cq, x0, orbit, 0, backend)
    acc = ntuple(_ -> CompensatedSum{T}(), Val(k)); Emax = _emax(Q0); n = 0
    while n < N
        z, A, Q, inc = m == 1 ? _step_heun(cq, z, Q, A, p, h, orbit, n, backend) : _step_rk4(cq, z, Q, A, p, h, orbit, n, backend)
        n += m; acc = _addall(acc, inc)
        E = _emax(Q); Emax = max(Emax, E)
        (isfinite(E) && E <= cq.orthogonality_tol) || throw(NumericalHealthError("orthogonality", Float64(E), cq.orthogonality_tol, n, orbit, backend))
        _record!(series, every, m, n, acc, E, Q)
    end
    l = _totals(acc); t = N * h
    return z, Q, (n = N, h = cq.map.tau, t = t, l = l, lambda = l / t, E = _emax(Q), Emax = Emax, nrhs = 2N, njac = N + 1)
end

"""
    propagate_cqr_orbit(cq, x0, Q0, N; orbit = 1, backend = :cpu) -> (x, Q, final)

One orbit of `N` orbit steps (`N` even for order 4). `final` has `n, h, t, l, lambda, E, Emax, nrhs, njac` (`l` the integrated growth
`l_i`, `lambda = l / t`). No allocation after warm-up. Throws `NumericalHealthError` ("state", "deviation matrix", "orthogonality").
"""
propagate_cqr_orbit(cq::ContinuousQR, x0::SVector, Q0::SMatrix, N::Integer; orbit::Integer = 1, backend::Symbol = :cpu) =
    _cqr_run(cq, x0, Q0, Int(N), Int(orbit), backend, nothing, 1)

"""
    propagate_cqr_series(cq, x0, Q0, N; every, orbit = 1) -> (x, Q, final, series)

The same run as [`propagate_cqr_orbit`](@ref), also recording `(n, l, E, Q)` every `every` Q-steps; the series is a bitwise prefix of the run.
"""
function propagate_cqr_series(cq::ContinuousQR, x0::SVector, Q0::SMatrix, N::Integer; every::Integer, orbit::Integer = 1, backend::Symbol = :cpu)
    every >= 1 || throw(ArgumentError("every must be at least 1, got $every"))
    series = NamedTuple[]
    x, Q, fin = _cqr_run(cq, x0, Q0, Int(N), Int(orbit), backend, series, Int(every))
    return x, Q, fin, series
end

function _cqr_provenance(cq::ContinuousQR, ::Val{d}, ::Val{k}, ::Type{T}, choice, seed) where {d, k, T}
    return (backend = choice.backend, device = choice.device, backend_fallback_reason = choice.reason, precision = precision_label(T), seed = seed, k = k,
            cqr_h = cq.map.tau, cqr_integrator = method_name(cq.map.method), cqr_integrator_order = nominal_order(cq.map.method), cqr_order = cq.order,
            cqr_form = _cqr_form(cq, Val(d), Val(k)), cqr_penalty_factor = _cqr_cnu(cq, Val(d), Val(k)), cqr_orthogonality_tol = cq.orthogonality_tol,
            cqr_flow = sprint(show, cq.flow))
end

"""
    propagate_cqr(cq, xs, N; k, seed, backend = :cpu, chunk = nothing, Q0 = nothing) -> PropagationResult

Many orbits, one per element of `xs`, each from its own frame `Q0[j]` (default `initial_deviations(seed, j, Val(d), Val(k), T)`).
`finals[j]` equals the `final` of the orbit alone; `W` holds the final frames. CPU only: a `:gpu` request falls back to the CPU with its reason
in the provenance (the device path is stage 04). `chunk` is accepted for the common signature and does not change the result.
"""
function propagate_cqr(cq::ContinuousQR, xs::AbstractVector{<:SVector{s, T}}, N::Integer; k::Integer, seed::Integer, backend::Symbol = :cpu, chunk = nothing,
                       Q0 = nothing) where {s, T}
    d = dimension(cq.flow)
    _check_cqr(cq, Val(d), Val(Int(k)), N)
    choice = resolve_backend(backend, T; cpu_only_reducer = true)
    Qs = Q0 === nothing ? [initial_deviations(seed, j, Val(d), Val(Int(k)), T) for j in eachindex(xs)] : Q0
    runs = [_cqr_run(cq, xs[j], Qs[j], Int(N), j, :cpu, nothing, 1) for j in eachindex(xs)]
    return PropagationResult([r[1] for r in runs], [r[2] for r in runs], [r[3] for r in runs], _cqr_provenance(cq, Val(d), Val(Int(k)), T, choice, seed))
end

"`cqr_bytes(cq, Val(d), Val(k), T)`: working memory of one orbit, `sizeof(T) (d + 4dk + 2k + 3 (d/2)^2) + 16` bytes (spec C21)."
cqr_bytes(::ContinuousQR, ::Val{d}, ::Val{k}, ::Type{T}) where {d, k, T} = sizeof(T) * (d + 4d * k + 2k + 3 * (d ÷ 2)^2) + 16
