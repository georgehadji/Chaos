# Study p3-abstention (Amendment D95): the store of the test run as `count-round` and `analyze` read it (E61 f, E62 c). The root attributes of every chunk file are read and checked;
# the datasets of the ref arm are read only for the rounds that a command needs (I39); a closure arm is never opened by `count-round`.

# The root attributes of every chunk file of the store, with the checks of E61 f: readable root, the study, the panel and preregistration hashes, no two chunks with the same cell, arm and item.
function p3t_scan(store, slug, panel_hash, prereg_hash, on_read)
    isdir(store) || _fail("store $store does not exist")
    out, seen = NamedTuple[], Dict{Tuple, String}()
    for f in sort!([f for f in readdir(store) if occursin(r"^chunk-\d{5}\.h5$", f)])
        path = joinpath(store, f)
        a = try
            Experiments.read_chunk_attrs(path)
        catch e
            _fail("chunk file $f: the root attributes cannot be read ($(first(sprint(showerror, e), 120)))")
        end
        on_read === nothing || on_read(:attrs, path)
        for k in ("study", "panel_sha256", "prereg_sha256", "chunk_key", "chunk_hash", "arms", "status", "complete", "method", "precision", "horizon", "tau")
            haskey(a, k) || _fail("chunk file $f: the root attribute $k is missing")
        end
        a["study"] == slug || _fail("chunk file $f belongs to study \"$(a["study"])\", not to study $slug")
        a["panel_sha256"] == panel_hash || _fail("panel_sha256 of chunk file $f differs from the SHA-256 of the panel file")
        a["prereg_sha256"] == prereg_hash || _fail("prereg_sha256 of chunk file $f differs from the preregistration header hash")
        key = (String(a["chunk_key"]), join(sort(String.(a["arms"])), ","))
        haskey(seen, key) && _fail("chunk files $(seen[key]) and $f have the same cell, arm and item ($(key[1]))")
        seen[key] = f
        push!(out, _p3t_chunk(f, path, a))
    end
    return out
end

# One chunk file as a record: arm ref or closure, the cell, the round, the K index and the class of a cell of the form of E74 (`k<K index>:<class>:r<round>`).
function _p3t_chunk(f, path, a)
    cell = something(match(r"cell=(\S+)", String(a["chunk_key"])), (captures = [""],)).captures[1]
    m = match(r"^k(\d+):([A-Za-z0-9]+):r(\d+)$", cell)
    arms = String.(a["arms"])
    "ref" in arms && m === nothing && _fail("chunk file $f: the cell \"$cell\" of a ref chunk is not of the form k<K index>:<class>:r<round>")
    return (file = f, path = path, ordinal = parse(Int, f[7:11]), key = String(a["chunk_key"]), hash = String(a["chunk_hash"]), arms = arms, ref = "ref" in arms, cell = String(cell),
            kindex = m === nothing ? 0 : parse(Int, m.captures[1]), class = m === nothing ? "" : String(m.captures[2]), round = m === nothing ? 0 : parse(Int, m.captures[3]),
            method = String(a["method"]), status = String(a["status"]), complete = a["complete"] === true, attrs = a)
end

"The panel file of the test run: the orbit tables, each with the round (D94 R23) and the flag `in_D`; an `ArgumentError` naming the orbit and the key otherwise (E61 e)."
function p3t_read_panel(panel, r_max)
    isfile(panel) || _fail("panel file $panel does not exist")
    orbits = TOML.parsefile(panel)["orbit"]
    for o in orbits
        haskey(o, "round") && _isint(o["round"]) || _fail("panel file $panel: orbit $(get(o, "id", "?")) has no round (D94 R23)")
        haskey(o, "in_D") && o["in_D"] isa Bool || _fail("panel file $panel: orbit $(get(o, "id", "?")) has no flag in_D (D94 R23)")
        for k in ("id", "k_index", "K", "x0_dec", "p0_dec")
            haskey(o, k) || _fail("panel file $panel: orbit $(get(o, "id", "?")) has no key $k")
        end
        1 <= o["round"] <= r_max || _fail("panel file $panel: orbit $(o["id"]) has the round $(o["round"]), outside 1 to r_max = $(Int(r_max))")
    end
    return orbits
end

# The cell of E74 of an orbit of the panel.
p3t_cell_of(o) = "k$(Int(o["k_index"])):$(o["class"]):r$(Int(o["round"]))"

# The cells expected for round r from the panel (E62 c, E74): one per (K index, class) of the orbits of that round; both L1 and L5 chunks are expected in each.
p3t_cells(panel_orbits, r) = sort!(unique([p3t_cell_of(o) for o in panel_orbits if Int(o["round"]) == r]))

# Whether the ref chunks of round r are all present and complete; returns the problems as texts (the missing cells are named).
function p3t_round_problems(chunks, panel_orbits, r)
    ref = [c for c in chunks if c.ref && c.round == r]
    probs = String[]
    for cell in p3t_cells(panel_orbits, r), m in ("L1", "L5")
        c = findfirst(x -> x.cell == cell && x.method == m, ref)
        c === nothing ? push!(probs, "cell $cell item $m is missing") : ref[c].complete || push!(probs, "cell $cell item $m (chunk file $(ref[c].file)) is not complete (complete is false)")
    end
    return probs
end

"The digest of E62 a: the SHA-256 of the chunk hashes of the ref chunks of the rounds up to `r`, ascending as byte strings, each followed by LF."
p3t_ref_digest(chunks, r) = bytes2hex(sha256(codeunits(join(sort!([c.hash for c in chunks if c.ref && c.complete && 1 <= c.round <= r]) .* "\n"))))

# ---------------------------------------------------------------- the ref state: one row per orbit with K > 0 of the rounds read

_p3t_xl(it) = it.mode == "wrapped" ? _p3_lifted_p2 : _p3_lifted_p1

# The row of one orbit from its ref items (E62 b, E64 a): excluded if an item or a recorded step is missing, the drift flag of L1 or L5, the label rho and lambda_ref, lambda_ns (exact rationals, `nothing` outside E_test).
function _p3t_ref_row(o, it1, it5, q, p2)
    id, N = Int(o["id"]), q.N_1
    base = (id = id, round = Int(o["round"]), kindex = Int(o["k_index"]), excluded = true, drift = false, rho = :missing, lref = nothing, lns = nothing)
    (it1 === nothing || it5 === nothing) && return base
    c0 = _p3_c0_inputs(it5, it1, N, _p3t_xl(it5))
    rn = _p3_step(it1, q.n_s)
    (c0.data === nothing || rn === nothing) && return base
    drift = "drift" in it1.flags || "drift" in it5.flags
    rho = drift ? :missing : p3_relabel(c0.data, p2).rho
    need = rho === :C || rho === :R      # only the orbits that stay in E_test need the exponents (E64 b, c)
    return merge(base, (excluded = false, drift = drift, rho = rho, lref = need ? _p3t_lam(it1.lam[_p3_step(it1, N)], id, "lambda_ref") : nothing, lns = need ? _p3t_lam(it1.lam[rn], id, "lambda_ns") : nothing))
end

# A stored exponent as an exact rational: the binary64 value itself (the owner's reading of "exact", M5); a non-finite value refuses naming the orbit.
function _p3t_lam(x, id, what)
    isfinite(x) || _fail("orbit $id: $what is $x, not a finite number")
    return Rational{BigInt}(x)
end

"The seven counts of E62 b after round `n` (a Dict by the names of the count lines) from the rows of the ref state and the targets."
function p3t_counts(rows, q, targets, n)
    rs = [r for r in rows if r.round <= n]
    ok = [r for r in rs if !r.excluded && !r.drift && r.rho === :C]
    fs, lw = p3_rat(q.f_s), p3_rat(q.lam_w2)
    ts(r) = r.lns <= fs * r.lref
    return Dict("orbits" => length(rs), "excluded" => count(r -> r.excluded, rs), "drift" => count(r -> !r.excluded && r.drift, rs), "t_s" => count(ts, ok),
                "w2" => count(r -> ts(r) || r.lref <= lw, ok), "target_t_s" => targets.t_s, "target_w2" => targets.w2)
end
