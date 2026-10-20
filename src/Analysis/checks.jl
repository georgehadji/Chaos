# Input checks before any output (E17) and the consistency checks of the stored facts (invariant I4). Every failure is an ArgumentError that
# names the cause; nothing has been written when it is raised.

_fail(msg) = throw(ArgumentError(msg))
_chunk_name(ordinal) = basename(chunk_path("", ordinal))

function _check_prereg(prereg, man)
    h = try
        check_prereg(prereg).header_sha256
    catch e
        e isa ArgumentError ? _fail("preregistration check failed: " * e.msg) : rethrow()
    end
    h == get(get(man, "prereg", Dict()), "header", "") || _fail("preregistration hash $h differs from the hash in the manifest")
    return h
end

function _check_panel(panel, man)
    isfile(panel) || _fail("panel file $panel does not exist")
    h = bytes2hex(sha256(read(panel)))
    pm = get(man, "panel", Dict())
    h == get(pm, "sha256_end", "") || _fail("panel file hash $h differs from the hash in the manifest")
    get(pm, "sha256_start", h) == h || _fail("the panel file hash at the start of the campaign differs from its hash at the end in the manifest")
    return h
end

function _check_manifest(man)
    get(man, "schema_version", nothing) == 2 || _fail("manifest has schema_version $(get(man, "schema_version", "absent")), not 2")
    get(man, "status", "") == "complete" || _fail("manifest status is \"$(get(man, "status", "absent"))\", not complete")
    haskey(man, "chunks") && haskey(man, "study") || _fail("manifest has no chunk list or no study name")
    for c in man["chunks"]
        c["status"] == "ok" && continue
        _fail("chunk $(c["ordinal"]) has status \"$(c["status"])\" in the manifest" *
              (c["status"] == "failed" ? "; chunks that failed twice are not analysed by this version (the study has none)" : ""))
    end
end

function _check_chunk(store, c, prereg_hash, panel_hash)
    path = chunk_path(store, c["ordinal"])
    isfile(path) || _fail("chunk file $(_chunk_name(c["ordinal"])) is missing from the store $store")
    a = Experiments.read_chunk_attrs(path)
    get(a, "pilot", false) === true && _fail("$(_chunk_name(c["ordinal"])) is a pilot file")
    get(a, "complete", false) === true || _fail("$(_chunk_name(c["ordinal"])) is not complete")
    get(a, "schema_version", 0) == 2 || _fail("$(_chunk_name(c["ordinal"])) has schema_version $(get(a, "schema_version", "absent")), not 2")
    get(a, "chunk_hash", "") == c["hash"] || _fail("chunk_hash of $(_chunk_name(c["ordinal"])) differs from the manifest")
    get(a, "prereg_sha256", prereg_hash) == prereg_hash || _fail("prereg_sha256 of $(_chunk_name(c["ordinal"])) differs from the preregistration hash")
    get(a, "panel_sha256", panel_hash) == panel_hash || _fail("panel_sha256 of $(_chunk_name(c["ordinal"])) differs from the panel file hash")
end

function _check_store(man, store, prereg_hash, panel_hash)
    isdir(store) || _fail("store $store does not exist")
    for c in man["chunks"]
        _check_chunk(store, c, prereg_hash, panel_hash)
    end
    listed = Set(_chunk_name(c["ordinal"]) for c in man["chunks"])
    for f in readdir(store)
        occursin(r"^chunk-\d{5}\.h5$", f) && !(f in listed) && _fail("store file $f is not in the manifest")
    end
end

"""
    check_inputs(prereg, manifest, panel, params, store) -> (man, params, hashes)

Read and check the four inputs and the store (E17). The `stale` area of the store is not looked at.
"""
function check_inputs(prereg, manifest, panel, params_file, store)
    p = read_params(params_file)
    isfile(manifest) || _fail("manifest $manifest does not exist")
    man = read_json(manifest)
    _check_manifest(man)
    ph = _check_prereg(prereg, man)
    pnl = _check_panel(panel, man)
    _check_store(man, store, ph, pnl)
    raw = Dict(f => bytes2hex(sha256(read(f))) for f in (prereg, manifest, panel, params_file))
    hashes = (prereg_header_sha256 = ph, panel_sha256 = pnl, manifest_sha256 = bytes2hex(sha256(read(manifest))),
              params_sha256 = bytes2hex(sha256(read(params_file))))
    return (man = man, p = p, hashes = hashes, panel_path = panel, raw = raw)
end

# Non-finite values in a used row, and the drift flag against the recomputation from the stored points (I4). Failed chunks are blanked, not used.
function check_items(items, man, p)
    failed = Set(c["ordinal"] for c in man["chunks"] if c["status"] == "failed")
    for it in items
        it.chunk in failed && continue
        "breakdown_tau" in it.flags && continue       # a result of the tau sweep (A5 item 2): the rows after the breakdown are blank (t = 0, NaN)
        all(isfinite, it.lam) && all(isfinite, it.state) || _fail("non-finite value in $(_chunk_name(it.chunk)), orbit $(it.orbit), method $(it.method)")
        it.t[end] == it.horizon || _fail("$(_chunk_name(it.chunk)), orbit $(it.orbit): last recorded step $(it.t[end]) is not the horizon $(it.horizon)")
        pp = it.winding === nothing ? view(it.state, :, 2) : it.state[:, 2] .+ 2pi .* it.winding[:, 2]      # the rule of the runner (wrapped mode adds the winding)
        drifted = any(abs(x) > p.drift_bound for x in pp)
        drifted == ("drift" in it.flags) || _fail("drift flag of orbit $(it.orbit) in $(_chunk_name(it.chunk)) differs from the recomputation from the stored points")
    end
end

"Refuse the run if an input file changed while it ran (I2): the SHA-256 of each file against the one taken at the start."
function recheck_inputs(raw, files...)
    for f in files
        bytes2hex(sha256(read(f))) == raw[f] || _fail("input file $f changed during the run")
    end
end
