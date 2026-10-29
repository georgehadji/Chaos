# Study p2-indicators (Amendment D89): the items of the reported quantities T9, T11, T14 to T17 (E31, E35).

_vk(d) = join(["$X $(String(v))" for (X, v) in sort!(collect(d); by = kv -> something(findfirst(==(kv[1]), P2_I8), 99))], ", ")

p2_t9(f4) = _tbl(f4.columns, [r for r in f4.rows if r[4] == "pooled"])

function p2_t11(sens, p)
    rows = [[String(pair), k, p.budgets[k + 1], s.n, s.A, s.kappa, s.median_diff, s.max_diff] for ((pair, k), s) in sort!(collect(sens); by = kv -> (string(kv[1][1]), kv[1][2]))]
    return _tbl(["pair", "k", "budget", "n", "A", "kappa", "median |dD|", "largest |dD|"], rows)
end

function p2_t14(range, leaving, below, orb)
    rows = Vector{Any}[]
    for ((X, kidx, k), r) in sort!(collect(range); by = kv -> (findfirst(==(kv[1][1]), P2_ALL), kv[1][2], kv[1][3]))
        push!(rows, Any["range", String(X), kidx, k, nothing, r === nothing ? 0 : r.n, r === nothing ? nothing : r.lo, r === nothing ? nothing : r.hi, r === nothing ? nothing : r.median, nothing])
    end
    for ((id, X), v) in sort!(collect(leaving); by = kv -> (kv[1][1], findfirst(==(kv[1][2]), P2_ALL)))
        push!(rows, Any["leaving", String(X), orb[id].K, nothing, id, nothing, nothing, nothing, nothing, v])
    end
    for ((X, kidx), v) in sort!(collect(below); by = kv -> (findfirst(==(kv[1][1]), P2_ALL), kv[1][2]))
        push!(rows, Any["below_median", String(X), kidx, nothing, nothing, nothing, nothing, nothing, nothing, v])
    end
    return _tbl(["part", "indicator", "K", "k", "orbit", "n", "lo", "hi", "median", "value"], rows)
end

function p2_t15(byrho, p)
    rows = Vector{Any}[]
    for ((name, pair, k), d) in sort!(collect(byrho); by = kv -> (kv[1][1], string(kv[1][2]), kv[1][3])), ((vx, vy, rho), n) in sort!(collect(d); by = kv -> string(kv[1]))
        push!(rows, Any[name, "$(pair[1]) with $(pair[2])", k, p.budgets[k + 1], String(vx), String(vy), String(rho), n])
    end
    return _tbl(["set", "package", "k", "budget", "verdict first", "verdict second", "rho", "orbits"], rows)
end

function p2_t16(flag_counts, omissions, controls, k0)
    rows = Vector{Any}[]
    for (s, d) in sort!(collect(flag_counts); by = first), (f, n) in sort!(collect(d); by = first)
        push!(rows, Any["flag", s, f, n, ""])
    end
    for (s, ids) in sort!(collect(omissions); by = first)
        push!(rows, Any["omission", s, "orbits", join(ids, ", "), ""])
    end
    for c in controls
        push!(rows, Any["control", c.id, c.scope === nothing ? "no scope" : c.scope, c.deviation,
                        "exponent $(_fmt(c.exponent)); lambda_1(N_1) $(_fmt(c.lambda)); stored deviation $(_fmt(c.stored_deviation)); verdicts at the last budget: $(_vk(c.verdicts9))"])
    end
    for (W, d) in sort!(collect(k0); by = first)
        push!(rows, Any["k0line", W, "largest circle distance of nu1 to p0 / (2 pi)", d.fq1,
                        "pa1: $(_fmt(d.pa1)); orbits left out as non-finite: fq1 $(d.fq1_nonfinite), pa1 $(d.pa1_nonfinite)"])
    end
    return _tbl(["part", "key", "name", "value", "detail"], rows)
end

function p2_t17(cost)
    return _tbl(["spec", "measured ns per orbit-step", "budget line ns", "bytes per orbit", "transform seconds"],
                [[s, c.measured_ns, c.budget_line_ns, c.bytes_per_orbit, c.transform_seconds] for (s, c) in sort!(collect(cost); by = first)])
end
