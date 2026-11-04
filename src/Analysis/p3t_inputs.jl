# Study p3-abstention (Amendment D95): the input checks that both commands share (E59, E61 a to g) and the ones of `analyze` alone (E61 h, i). Nothing is written by any of them.

# A path inside the data tree of the repository refuses (I38); `what` names the option.
function p3t_check_not_data(path, what)
    root = lowercase(replace(abspath(joinpath(pkgdir(parentmodule(@__MODULE__))::String, "data")), '\\' => '/'))
    o = lowercase(replace(abspath(path), '\\' => '/'))
    (o == root || startswith(o, root * "/")) && _fail("$what $path is inside the data folder; nothing is written into data/")
end

# The SHA-256 of every input file by name (E61 g); the check at the end refuses a file that changed while the command ran.
p3t_recheck(raw) = for (f, h) in raw
    isfile(f) && bytes2hex(sha256(read(f))) == h || _fail("input file $f changed during the run")
end

"""
    p3t_context(; ...) -> NamedTuple

The checks a to g of E61 common to `count-round` and `analyze`, in the order of the rules: the slug, the hash rule, the 71-key table and the test-side table, the p2 manifest and table, the
panel file with its rounds, the root attributes of the store. Returns the tables, the panel, the chunk records and the hashes taken at the start.
"""
function p3t_context(; prereg, params, test_params, store, panel, p2_manifest, p2_prereg, p2_params, calibration_record, on_read, command, files = String[], p1_params = nothing, after_hash = nothing)
    slug = study_slug(prereg)
    slug == P3_SLUG || _fail("$command is the command of study $P3_SLUG; this preregistration is of study \"$slug\"")
    raw = _p3_hashes(vcat([prereg, params, test_params, calibration_record, p2_manifest, p2_prereg, p2_params, panel], files))
    ph = _p2_check_prereg(prereg)
    after_hash === nothing || after_hash()      # analyze: the round log is the first input after the hash rule
    qp = _p3_read(params, prereg, p2_params, p2_prereg, p2_manifest, p1_params, nothing)
    t = p3t_read_params(test_params; prereg, record = calibration_record, params)
    orbits = p3t_read_panel(panel, Int(qp.q.r_max))
    panel_hash = raw[panel]
    chunks = p3t_scan(store, slug, panel_hash, ph, on_read)
    merge!(raw, _p3_hashes([c.path for c in chunks]))
    return (slug = slug, q = qp.q, p2 = qp.p2, t = t, panel_orbits = orbits, panel_hash = panel_hash, prereg_hash = ph, chunks = chunks, raw = raw, targets = (t_s = t.target_t_s, w2 = t.target_w2))
end

# The inputs of `analyze` beyond the context: the round log is present, non-empty and closed (E61 h); returns its text.
function p3t_read_closed_log(path)
    isfile(path) || _fail("the round log $path is absent: analyze needs the closed round log of the test run")
    text = read(path, String)
    isvalid(text) || _fail("the round log $path is not valid UTF-8")
    isempty(text) && _fail("the round log $path is empty: analyze needs the closed round log of the test run")
    return text
end

# E76: the manifest status is `complete`, or `complete with failed chunks` exactly when a chunk entry has status failed.
function _p3t_check_status(man)
    st = String(get(man, "status", "absent"))
    nf = count(e -> get(e, "status", "ok") == "failed", man["chunks"])
    ok = st == "complete" ? nf == 0 : st == "complete with failed chunks" && nf > 0
    ok || _fail("manifest status is \"$st\" with $nf failed chunk entries: the status must be complete (no failed entry) or complete with failed chunks (at least one)")
end

# E61 i, E76, E81: the manifest of the test run: the study, the status, the six arm names of the run block, and every chunk file against its manifest entry.
function p3t_check_manifest(man, c, store)
    get(man, "study", nothing) == c.slug || _fail("manifest study \"$(get(man, "study", "absent"))\" is not the study $(c.slug)")
    _p3t_check_status(man)
    arms = get(man, "arms", Dict())
    for a in P3T_ARMS
        haskey(arms, a) || _fail("the manifest has no arm \"$a\" in its run block (the arms of the test run are $(join(P3T_ARMS, ", ")))")
    end
    byfile = Dict(x.file => x for x in c.chunks)
    for e in man["chunks"]
        f = basename(chunk_path("", _p3t_entry(e, "ordinal")))
        haskey(byfile, f) || _fail("chunk file $f is listed in the manifest but missing from the store")
        byfile[f].hash == _p3t_entry(e, "hash") || _fail("chunk hash of $f differs from the hash in the manifest")
        (byfile[f].status == get(e, "status", "ok")) || _fail("chunk $f has status \"$(byfile[f].status)\" in the file and \"$(get(e, "status", "ok"))\" in the manifest")
    end
    listed = Set(basename(chunk_path("", e["ordinal"])) for e in man["chunks"])
    for x in c.chunks
        x.file in listed || _fail("store file $(x.file) is not in the manifest")
    end
end

const P3T_ARMS = ("ref", "ind", "ind_pa1", "d64", "d64_pa1", "chk")

# A field of a manifest chunk entry that analyze reads; a missing field is a refusal naming it.
_p3t_entry(e, k) = haskey(e, k) ? e[k] : _fail("the manifest chunk entry with ordinal $(get(e, "ordinal", "absent")) has no field \"$k\"")

# E61 l: a chunk of a closure arm that is not complete and not a chunk that failed twice.
function p3t_check_closure_complete(c)
    for x in c.chunks
        (x.ref || x.complete || x.status == "failed") && continue
        _fail("the chunk file $(x.file) of a closure arm is not complete (complete is false) and is not listed as failed twice")
    end
end
