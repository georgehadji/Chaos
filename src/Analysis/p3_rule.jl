# Study p3-abstention (Amendment D93): the rule R_a and R_ab (E52), the search space Theta (E52), the cost model and the matched budgets (E53). Pure functions of stored
# verdicts and of the parameter table; nothing here reads a file. Verdicts are `:chaotic`, `:regular`, `:undecided`; a vector of ten verdicts is indexed by the budget
# index k = 0 to 9 (position k + 1). The output of a rule is `:chaotic_supported`, `:regular_supported` or `:unresolved`.

"""
    p3_stable(v, k, w) -> Bool

The verdict of one indicator at budget index `k` is stable at depth `w` iff it is decided and equals the verdicts at `k - 1, ..., k - w` (E51; `w = 0`: decided).
An `ArgumentError` naming the budget if `b_k / 2^w` is not an element of B (an error of the caller).
"""
function p3_stable(v::AbstractVector{Symbol}, k::Integer, w::Integer)
    k - w >= 0 || throw(ArgumentError("budget b_$k / 2^$w is not an element of B (the evaluation needs k - w >= 0)"))
    c = v[k + 1]
    c === :undecided && return false
    for j in 1:w
        v[k + 1 - j] === c || return false
    end
    return true
end

"""
    p3_evidence(V, k, theta) -> (nC, nR, nCany)

Evidence (a) of E52 for the verdict table `V` (indicator => ten verdicts) at budget index `k` and `theta.V`, `theta.w`: the chaotic and the regular votes that are
stable at depth `w`, and the chaotic votes stable or not.
"""
function p3_evidence(V, k::Integer, theta)
    nC = nR = nA = 0
    for X in theta.V
        v = V[X]
        c = v[k + 1]
        s = p3_stable(v, k, theta.w)
        c === :chaotic && (nA += 1; s && (nC += 1))
        c === :regular && s && (nR += 1)
    end
    return (nC = nC, nR = nR, nCany = nA)
end

"The output of the rule form R_a (E52): chaotic-supported iff n_C >= q_C, else regular-supported iff n_Cany <= v and n_R >= q_R, else unresolved."
function p3_rule_a(V, k::Integer, theta)
    e = p3_evidence(V, k, theta)
    e.nC >= theta.qC && return :chaotic_supported
    (e.nCany <= theta.v && e.nR >= theta.qR) && return :regular_supported
    return :unresolved
end

"""
    p3_rule_ab(V, Vd, k, theta) -> Symbol

R_ab: R_a, overridden to unresolved if for some X in `theta.V` the Float64 verdict and the Double64 verdict of `Vd` at `k` are both decided and differ. `Vd` must hold
every indicator of `theta.V`; otherwise an `ArgumentError` naming Double64 (never silently R_a).
"""
function p3_rule_ab(V, Vd, k::Integer, theta)
    (Vd !== nothing && all(X -> haskey(Vd, X), theta.V)) || throw(ArgumentError("R_ab is not defined without Double64 verdicts of every indicator of V"))
    out = p3_rule_a(V, k, theta)
    for X in theta.V
        a, b = V[X][k + 1], Vd[X][k + 1]
        a !== :undecided && b !== :undecided && a != b && return :unresolved
    end
    return out
end

# ---------------------------------------------------------------- the search space (E52)

p3_space_size(q) = sum(binomial(length(q.voter_bits), s) * s * (s - s ÷ 2) * length(q.w_range) * length(q.v_range) for s in q.min_voters:length(q.voter_bits))

"""
    p3_enumerate(q) -> Vector{NamedTuple}

The elements `(index, V, mask, qC, qR, w, v)` of Theta in the order of E52: by `|V|` ascending, the bitmask of `V` (bit i - 1 for the i-th voter) ascending, then `q_C`,
`q_R`, `w`, `v`; `V` lists its indicators in the bit order.
"""
function p3_enumerate(q)
    names = Symbol.(q.voter_bits)
    n = length(names)
    out = NamedTuple[]
    for sz in q.min_voters:n
        for m in sort!([m for m in 0:(2^n - 1) if count_ones(m) == sz])
            V = Symbol[names[i] for i in 1:n if (m >> (i - 1)) & 1 == 1]
            for qC in 1:sz, qR in (sz ÷ 2 + 1):sz, w in q.w_range, v in q.v_range
                push!(out, (index = length(out) + 1, V = V, mask = m, qC = qC, qR = qR, w = w, v = v))
            end
        end
    end
    return out
end

# ---------------------------------------------------------------- the cost model (E53); all costs and products are exact rationals of the decimal texts

_p3_spec(q, X) = q.spec_of[string(X)]
_p3_cost(q, spec) = p3_rat(getfield(q, Symbol("cost_", spec)))

"C_R(theta, b) = b times the sum of the costs of the distinct specs needed by `V` (a list of indicators), exact."
p3_rule_cost(q, V, b::Integer) = b * sum(_p3_cost(q, s) for s in unique(_p3_spec(q, X) for X in V))

"The cost cap of F1 (C4): the minimum over the voters of the cost of its spec times `N_c`."
p3_cost_cap(q) = minimum(_p3_cost(q, _p3_spec(q, X)) for X in Symbol.(q.voter_bits)) * q.N_c

"""
    p3_matched(q, V) -> (CR, X, le_read, le_unmatched)

`CR = C_R(theta, b_P)`; `X[ind] = (b, capped, round_down, round_down_none)` for the voters and `:LE`: the smallest element of B with `b cost >= CR` (the ceiling rule;
`b = nothing` when capped) and the largest with `b cost <= CR` (`nothing` with the flag when none). `le_read = min(b_LE, N_c)` (`N_c` if capped), flagged
`le_unmatched` iff its read budget times its cost is below `CR`.
"""
function p3_matched(q, V)
    CR = p3_rule_cost(q, V, q.b_primary)
    out = Dict{Symbol, NamedTuple}()
    for X in vcat(Symbol.(q.voter_bits), :LE)
        c = _p3_cost(q, _p3_spec(q, X))
        i = findfirst(b -> b * c >= CR, q.budgets)
        j = findlast(b -> b * c <= CR, q.budgets)
        out[X] = (b = i === nothing ? nothing : q.budgets[i], capped = i === nothing, round_down = j === nothing ? nothing : q.budgets[j], round_down_none = j === nothing)
    end
    le = out[:LE]
    read = le.capped ? q.N_c : min(le.b, q.N_c)
    return (CR = CR, X = out, le_read = read, le_unmatched = read * _p3_cost(q, _p3_spec(q, :LE)) < CR)
end

"""
    p3_d64_costs(q, V) -> (pair, float, ratio)

The Double64 costs (reported, E53): the pair cost `b_P` times the sum over the distinct specs of `cost + ratio cost`, the Float64 vote cost, and `float // pair`.
"""
function p3_d64_costs(q, V)
    specs = unique(_p3_spec(q, X) for X in V)
    ratio(s) = p3_rat(getfield(q, Symbol("ratio_d64_", s)))
    pair = q.b_primary * sum(_p3_cost(q, s) * (1 + ratio(s)) for s in specs)
    fl = p3_rule_cost(q, V, q.b_primary)
    return (pair = pair, float = fl, ratio = fl // pair)
end
