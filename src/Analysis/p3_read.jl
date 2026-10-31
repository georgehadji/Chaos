# Study p3-abstention (Amendment D93): the inputs of `calibrate` (E46, E48): the paths by option or by the study layout, the checks that need the files only, and the reading of the
# items of the p2 and P1 stores that E49 and E51 read (nothing else is read, so nothing else is checked).

# The inputs by option or by the layout of the study (workspaces/study/stages): the files of studies p1-le-benchmark and p2-indicators next to the preregistration.
function _p3_inputs(prereg, params, out, kw)
    stages = joinpath(dirname(abspath(prereg)), "..", "..")
    st(sub, name) = normpath(joinpath(stages, sub, name))
    pick(k, default) = get(kw, k, nothing) === nothing ? default : kw[k]
    p2_prereg = pick(:p2_prereg, joinpath(dirname(abspath(prereg)), P2_SLUG * "-preregistration.md"))
    p1_prereg = pick(:p1_prereg, joinpath(dirname(abspath(prereg)), P1_SLUG * "-preregistration.md"))
    for (k, _) in kw
        k in (:p2_store, :p2_manifest, :p2_prereg, :p2_params, :panel, :p1_store, :p1_manifest, :p1_prereg, :p1_params) || throw(ArgumentError("unknown keyword $k of calibrate"))
    end
    return (prereg = prereg, params = params, out = out === nothing ? dirname(abspath(prereg)) : out, p2_prereg = p2_prereg, p1_prereg = p1_prereg,
            p2_manifest = pick(:p2_manifest, st("03-run/output", P2_SLUG * "-manifest.json")), p2_params = pick(:p2_params, st("04-analyze/references", P2_SLUG * "-analysis-parameters.toml")),
            panel = pick(:panel, st("02-panel/output", P1_SLUG * "-panel.toml")), p1_manifest = pick(:p1_manifest, st("03-run/output", P1_SLUG * "-manifest.json")),
            p1_params = pick(:p1_params, st("04-analyze/references", P1_SLUG * "-analysis-parameters.toml")),
            p2_store = pick(:p2_store, isfile(p2_prereg) ? default_store(p2_prereg) : ""), p1_store = pick(:p1_store, isfile(p1_prereg) ? default_store(p1_prereg) : ""))
end

# The data tree of the repository: an output folder inside it refuses (I30).
function _p3_check_out(out)
    root = lowercase(replace(abspath(joinpath(pkgdir(parentmodule(@__MODULE__))::String, "data")), '\\' => '/'))
    o = lowercase(replace(abspath(out), '\\' => '/'))
    (o == root || startswith(o, root * "/")) && _fail("the output folder $out is inside the data folder; nothing is written into data/")
end

# E48 f: the P1 inputs (manifest complete, chunks ok with their hashes, the preregistration hash rule, the panel file) and the study of its manifest; the panel hash equals the p2 one.
function _p3_check_p1(i, p2man)
    c = check_inputs(i.p1_prereg, i.p1_manifest, i.panel, i.p1_params, i.p1_store)
    c.man["study"] == P1_SLUG || _fail("the P1 manifest study \"$(c.man["study"])\" is not \"$P1_SLUG\"")
    c.hashes.panel_sha256 == get(get(p2man, "panel", Dict()), "sha256_end", "") || _fail("the panel file hash of the P1 manifest differs from the panel hash of the p2 manifest")
    return c
end

# E49 f (base labelling of the P1 inputs, in memory, as D89 E23 a): the old class of every panel orbit that the base pipeline can label; an orbit whose items lack what the base
# pipeline needs has no old label (it is then missing for C0 too, or counted as such).
function _p3_p1_old(c, i)
    items = read_store(i.p1_manifest, i.p1_store)
    orbits = TOML.parsefile(i.panel)["orbit"]
    idx = index_items(items)
    ctrl = _controls(idx, orbits, c.p)
    ctrl.ok || _fail("the P1 reference run stops at a control defect: " * join(["control $(d.id) method $(d.method)" for d in ctrl.defects], "; "))
    flags, dirs, bf = _flag_map(items), _dirs_map(items, c.p), _bf_label(items)
    old = Dict{Int, Symbol}()
    for o in orbits
        o["control"] === true && continue
        r = try
            orbit_record(idx, o, c.p, flags, dirs, bf)
        catch e
            e isa Union{ArgumentError, BoundsError, KeyError} ? nothing : rethrow()    # an item that lacks a recorded step: no old label
        end
        r === nothing || (old[Int(o["id"])] = p2_rho(r.label.instability, r.label.structure))
    end
    return (items = items, idx = idx, old = old, orbits = orbits)
end

# The chunks of the p2 store that hold an item of E49 or E51: the spec, the precision, the horizon and the direction of the chunk attributes.
function _p3_wanted(p)
    d1 = p.directions[1]
    F, D = "Float64", "Double64"
    h = (L1 = p.N_1, L5 = p.N_1, tg = p.N_1, fq1 = p.N_c, bk1 = p.N_c, pa1 = p.N_c + 1)
    w = Set{Tuple}([("L1", F, h.L1, d1), ("L5", F, h.L5, d1), ("tg", F, h.tg, d1), ("fq1", F, h.fq1, nothing), ("bk1", F, h.bk1, nothing), ("pa1", F, h.pa1, nothing)])
    union!(w, [("L1", D, p.N_dir, d1), ("tg", D, p.N_dir, d1), ("fq1", D, p.N_c, nothing), ("bk1", D, p.N_c, nothing), ("pa1", D, p.N_c + 1, nothing)])
    return w
end

# The items of the wanted chunks, without the controls (no verdict of a control is read), grouped as `p2_check_items` expects.
function _p3_read_p2(man, store, p, panel_orbits)
    want = _p3_wanted(p)
    controls = Set(Int(o["id"]) for o in panel_orbits if o["control"] === true)
    items = NamedTuple[]
    for c in man["chunks"]
        path = chunk_path(store, c["ordinal"])
        a = Experiments.read_chunk_attrs(path)
        key = (String(a["method"]), String(a["precision"]), Int(a["horizon"]), haskey(a, "direction") ? String(a["direction"]) : nothing)
        key in want || continue
        append!(items, [it for it in _p2_chunk_items(read_chunk(path), c["ordinal"]) if !(it.orbit in controls)])
    end
    return items
end

# The verdict items of E51 for the orbits of Cal and of the gate set.
function _p3_verdict_requirements(p, ids)
    d1, F = p.directions[1], "Float64"
    return [("tg", F, p.N_1, d1, ids), ("fq1", F, p.N_c, nothing, ids), ("pa1", F, p.N_c + 1, nothing, ids), ("bk1", F, p.N_c, nothing, ids)]
end

"The SHA-256 of every input file (E55), by name, for the fingerprints and for the check that nothing changed during the run (E48 j)."
_p3_hashes(files) = Dict(f => open(io -> bytes2hex(sha256(io)), f) for f in files if isfile(f))

# The chunk files a manifest lists, in the store; none if the manifest cannot be read (the checks that follow refuse it).
function _p3_chunk_files(manifest, store)
    isfile(manifest) || return String[]
    try
        return String[chunk_path(store, c["ordinal"]) for c in read_json(manifest)["chunks"]]
    catch
        return String[]
    end
end
