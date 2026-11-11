# OFLI reducer (spec ofli, O1 to O13, Algorithm 1 to 7): the FLI of the component of the deviation columns orthogonal to the flow.

"""
    OFLI(flow; replace = true)

Orthogonal Fast Lyapunov Indicator reducer (log form) for `propagate` on the one-step map of a Hamiltonian flow `flow`
(`1 <= k <= d` columns, `d = dimension(flow)`; any other `k` or `d` throws `ArgumentError`). At each handoff `n` it evaluates the exact
vector field `f` of `flow` at the numerical state (candidate A, spec P1), reads the arriving, unnormalised columns, and for column `i`
forms `g_i = S_i + log(beta_i)`, `beta_i` the scaled norm of the column minus its component along `f`, `S_i` the compensated sum of
`log(alpha_i)` before this handoff (the same operations as [`FLI`](@ref) and [`SALI`](@ref), so the sums are bitwise equal to theirs).
The running maxima are `OFLI_B = max_i G_i` (basis) and `OFLI_1 = G_1` (first column); they range over the handoffs only, the
`s = 0` term is [`ofli_initial`](@ref). The companion unprojected maxima `FLI_B`, `FLI_1` (with the `s = 0` term 0) come from the same sums.

`replace` is the chain position of the one-replacer rule: with `replace = true` (alone, or last) `W` is replaced by its column-normalised
matrix; with `replace = false` the columns pass on unchanged for a later replacer. The reducer cannot see the chain or the map: it holds
`flow`, and a flow different from the one the map integrates gives wrong numbers silently (spec P10; the run record carries its name).

`reducer_final` gives `(n, ofli, ofli1, inst, fli, fli1, sum, comp, stepB, colB, step1, sinB, sin1, sinmin, nminf)`: `inst = max_i g_i(n)`
the instantaneous value, the maximisers (map step, column) and `sin theta = beta / alpha` at them, the least `sin theta` seen, and `nminf`
the count of minus infinity terms (a zero orthogonal component). A field that is zero or not finite is a health failure ("vector field"),
a zero or non-finite column norm one ("deviation matrix", the column).
"""
struct OFLI{F}
    flow::F
    replace::Bool
end
OFLI(flow; replace::Bool = true) = OFLI(flow, replace)

# `code` 0 = ok, 1 = zero or non-finite column norm `value` of column `column`, 2 = field norm `value` zero or not finite, at handoff `step`.
struct OFLIState{T, k}
    acc::SVector{k, CompensatedSum{T}}
    g1::T
    gB::T
    inst::T
    f1::T
    fB::T
    sin1::T
    sinB::T
    sinmin::T
    value::T
    nminf::Int
    n::Int
    step::Int
    step1::Int
    stepB::Int
    colB::Int8
    code::Int8
    column::Int8
end

function reducer_init(r::OFLI, ::Type{T}, ::Val{d}, ::Val{k}) where {T, d, k}
    dimension(r.flow) == d || throw(ArgumentError("OFLI: the flow has dimension $(dimension(r.flow)), the state has $d"))
    1 <= k <= d || throw(ArgumentError("OFLI needs 1 <= k <= d, got k = $k, d = $d"))
    acc = SVector{k, CompensatedSum{T}}(ntuple(_ -> CompensatedSum{T}(), Val(k)))
    ninf = T(-Inf)
    return OFLIState{T, k}(acc, ninf, ninf, ninf, zero(T), zero(T), T(Inf), T(Inf), T(Inf), zero(T), 0, 0, 0, 0, 0, Int8(0), Int8(0), Int8(0))
end

# Scaled norms (alpha, beta) of a column w and of its component orthogonal to the unit vector u.
function ofli_column(u::SVector, w::SVector)
    return safe_norm(w), safe_norm(w - dot(u, w) * u)
end

function ofli_fail(st::OFLIState{T, k}, n::Int, code::Int8, column::Int8, value) where {T, k}
    return OFLIState{T, k}(st.acc, st.g1, st.gB, st.inst, st.f1, st.fB, st.sin1, st.sinB, st.sinmin, T(value), st.nminf, st.n, n,
                           st.step1, st.stepB, st.colB, code, column)
end

# Update of the maxima and sums at handoff n from the checked column norms `cols` (tuple of (alpha, beta)).
function ofli_advance(st::OFLIState{T, k}, n::Int, cols::Tuple) where {T, k}
    newacc = SVector{k, CompensatedSum{T}}(ntuple(i -> add(st.acc[i], log(cols[i][1])), Val(k)))
    g1, gB, inst, f1, fB = st.g1, st.gB, T(-Inf), st.f1, st.fB
    sin1, sinB, sinmin, step1, stepB, colB, nminf = st.sin1, st.sinB, st.sinmin, st.step1, st.stepB, st.colB, st.nminf
    for i in 1:k
        alpha, beta = cols[i]
        g = beta > 0 ? total(st.acc[i]) + log(beta) : T(-Inf)    # S_i before this handoff
        sn, S = beta / alpha, total(newacc[i])
        sinmin = min(sinmin, sn)
        nminf += beta > 0 ? 0 : 1
        inst = max(inst, g)
        g > gB && ((gB, stepB, colB, sinB) = (g, n, Int8(i), sn))             # strict increase: ties keep the earlier maximiser
        i == 1 && g > g1 && ((g1, step1, sin1) = (g, n, sn))
        i == 1 && S > f1 && (f1 = S)
        S > fB && (fB = S)
    end
    return OFLIState{T, k}(newacc, g1, gB, inst, f1, fB, sin1, sinB, sinmin, st.value, nminf, n, st.step, step1, stepB, colB, st.code, st.column)
end

function reducer_step(r::OFLI, st::OFLIState{T, k}, n::Int, x::SVector, W::SMatrix{d, k}) where {T, d, k}
    st.code == 0 || return st, W
    f = vector_field(r.flow, x)
    fn = safe_norm(f)
    (isfinite(fn) && fn > floatmin(T)) || return ofli_fail(st, n, Int8(2), Int8(0), fn), W
    u = f / fn
    cols = ntuple(i -> ofli_column(u, W[:, i]), Val(k))
    for i in 1:k
        a = cols[i][1]
        (isfinite(a) && a > 0) || return ofli_fail(st, n, Int8(1), Int8(i), a), W
    end
    W1 = r.replace ? SMatrix{d, k}(ntuple(idx -> W[idx] / cols[(idx - 1) ÷ d + 1][1], Val(d * k))) : W
    return ofli_advance(st, n, cols), W1
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::OFLI, st::OFLIState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::OFLI, st::OFLIState) = st.code != 0

function reducer_error(::OFLI, st::OFLIState, orbit, backend)
    value = Float64(st.value)
    what = st.code == 2 ? "vector field" : "deviation matrix"
    threshold = st.code == 2 ? Float64(floatmin(typeof(st.value))) : (isfinite(value) ? 0.0 : Inf)
    return NumericalHealthError(what, value, threshold, st.step, orbit, backend, Int(st.column))
end

function reducer_final(::OFLI, st::OFLIState{T, k}) where {T, k}
    return (n = st.n, ofli = st.gB, ofli1 = st.g1, inst = st.inst, fli = st.fB, fli1 = st.f1,
            sum = SVector{k, T}(ntuple(i -> st.acc[i].s, Val(k))), comp = SVector{k, T}(ntuple(i -> st.acc[i].c, Val(k))),
            stepB = st.stepB, colB = Int(st.colB), step1 = st.step1, sinB = st.sinB, sin1 = st.sin1, sinmin = st.sinmin, nminf = st.nminf)
end

reducer_provenance(r::OFLI) = (indicator = "ofli", ofli_flow = string(nameof(typeof(r.flow))), ofli_replace = r.replace)

"""
    ofli_initial(flow, x0, W0) -> SVector{k}

The `s = 0` terms `log(||P w_i||)` of the columns of `W0` (`P` the projection orthogonal to `vector_field(flow, x0)`), minus infinity for a
zero component (spec O2, Algorithm 1). The OFLI series is the running maximum of these and of the reducer's handoff terms.
"""
function ofli_initial(flow, x0::SVector, W0::SMatrix{d, k, T}) where {d, k, T}
    f = vector_field(flow, x0)
    u = f / safe_norm(f)
    return SVector{k, T}(ntuple(i -> (b = ofli_column(u, W0[:, i])[2]; b > 0 ? log(b) : T(-Inf)), Val(k)))
end

"""
    check_frame(OFLI(flow), W0, x0; orbit = 0) -> defect

Frame check of a supplied initial frame (spec Algorithm 1, O13): every column must be a unit vector, `abs(norm(w_i) - 1) <= 8 effective_eps(T)`,
else `NumericalHealthError("initial frame", ...)` at step 0 naming the column; for a single column (`k = 1`) the column must also not be parallel
to `vector_field(flow, x0)` (its sine against the field above `8 effective_eps(T)`), because a start along the flow stays along it. A basis is
accepted with a column parallel to the flow. Returns the larger norm defect.
"""
function check_frame(r::OFLI, W0::AbstractMatrix{T}, x0::AbstractVector; orbit::Integer = 0) where {T}
    threshold = 8 * effective_eps(T)
    defect = zero(T)
    for i in axes(W0, 2)
        d = abs(safe_norm(W0[:, i]) - 1)
        (isfinite(d) && d <= threshold) || throw(NumericalHealthError("initial frame", Float64(d), Float64(threshold), 0, orbit, :cpu, i))
        defect = max(defect, d)
    end
    if size(W0, 2) == 1
        f = vector_field(r.flow, SVector{length(x0), T}(x0))
        w = SVector{size(W0, 1), T}(W0[:, 1])
        s = safe_norm(w - dot(f, w) / dot(f, f) * f) / safe_norm(w)
        s > threshold || throw(NumericalHealthError("initial frame", Float64(s), Float64(threshold), 0, orbit, :cpu, 1))
    end
    return defect
end
