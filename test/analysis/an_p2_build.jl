# Amendment D89, stage 02: assembly of the synthetic study "p2-indicators" (chunk files, manifest, inputs) and editing helpers. Definitions only.

const P2_SPEC_INFO = Dict("tg" => ("tangent", ["FLI", "MEGNO", "SALI"]), "pr" => ("pair", ["RLI", "MEGNO"]), "st" => ("stretching", ["stretching"]),
    "fq1" => ("frequency", ["frequency"]), "fq1_raw" => ("frequency", ["frequency"]), "fq1_qp" => ("frequency", ["frequency"]),
    "pa1" => ("phase_advance", ["phase_advance"]), "pa_planar" => ("phase_advance", ["phase_advance"]), "bk1" => ("birkhoff", ["birkhoff"]))

function p2_study(tag::AbstractString; variant::Symbol = :default, p = nothing, status = "complete", pilot = false)
    dir = an_path("p2-$tag")
    rm(dir; recursive = true, force = true)
    mkpath(dir)
    p1 = p2_p1_study(joinpath(dir, "p1"); cover = Dict(2 => Dict("L9" => 100.0)),
                     flags = (id, sp) -> (id == 13 && sp.method == "L5" && "main" in sp.arms) ? ["diagnostics exceeded"] : String[])
    panel = joinpath(dir, "$P2_SLUG-panel.toml")
    cp(p1.panel, panel)
    ph = an_sha256(read(panel))
    p = p === nothing ? p2_params(; panel_sha256 = ph) : merge(p, (; panel_sha256 = ph))
    text, h = p2_prereg_text(p)
    prereg = joinpath(dir, "$P2_SLUG-preregistration.md"); write(prereg, text)
    params = p2_write_params(joinpath(dir, "$P2_SLUG-analysis-parameters.toml"), p)
    design = p2_design(variant)
    h4 = p2_h4_tables(p)
    store = joinpath(dir, "store"); mkpath(store)
    specs = p2_specs(p)
    chunks = Dict{String, Any}[]
    total = 0
    for (ord, sp) in enumerate(specs)
        recs = AX.OrbitRecord[]
        for id in sp.orbits
            info = p2_info(id, h4)
            ex = Dict{String, Any}()
            if "lyap_p" in sp.arms && sp.role === :L1 && id in (17, 18, 20)
                ex["control_scope"] = id == 18 ? "info" : "bound"
                ex["control_deviation"] = abs(design.lref[id] - Float64(AN_ORBITS[id]["exponent"]))
            end
            push!(recs, p2_orbit_record(sp, info, design, p; flags = id in design.drift ? ["drift"] : String[], extra = ex))
        end
        total += length(recs)
        key = "method=$(sp.id) kind=$(sp.kind) precision=$(sp.precision) horizon=$(sp.horizon) tau=1 direction=$(sp.direction) mode=wrapped arm=$(sp.arms[1])"
        ch = "c"^60 * lpad(string(ord), 4, '0')
        root = Dict{String, Any}("schema_version" => 2, "study" => P2_SLUG, "chunk" => ord, "chunk_key" => key, "chunk_hash" => ch, "plan_sha256" => p.plan_sha256,
            "panel_sha256" => ph, "prereg_sha256" => h, "git_sha" => "a"^40, "code_version" => "an-code-1", "backend_request" => "cpu", "backend" => "cpu",
            "device" => "cpu", "backend_fallback_reason" => "", "precision" => sp.precision, "method" => sp.id, "horizon" => sp.horizon, "tau" => 1,
            "vectors" => sp.role in (:tg, :L5) ? 2 : 1, "coordinate_mode" => "wrapped", "grid_sha256" => "g"^64,
            "status" => "ok", "attempt" => 1, "setup_check" => "not applicable", "first_attempt_reasons" => String[], "failure_reasons" => String[],
            "pilot" => pilot, "arms" => sp.arms, "session" => 1, "control_stop" => false, "dirty" => false, "created_utc" => "2026-10-29T00:00:00Z")
        sp.role in (:L1, :L5, :tg, :st) && (root["direction"] = sp.direction)
        if sp.kind != ""
            kd, mem = P2_SPEC_INFO[sp.id]
            st_text = "indicator.id = \"$(sp.id)\"\nindicator.kind = \"$kd\"\n"
            merge!(root, Dict{String, Any}("kind" => kd, "spec" => st_text, "spec_sha256" => an_sha256(st_text), "chain" => [kd], "members" => mem))
            kd in ("frequency", "phase_advance", "birkhoff") && merge!(root, Dict{String, Any}("windows" => p.windows, "run_length" => sp.horizon))
            kd == "birkhoff" && (root["observable_names"] = P2_BKROWS)
        end
        AX.write_chunk(joinpath(store, "chunk-" * lpad(string(ord), 5, '0') * ".h5"), root, recs)
        ent = Dict{String, Any}("ordinal" => ord, "key" => key, "hash" => ch, "arms" => sp.arms, "status" => "ok", "attempt" => 1, "orbits" => length(recs), "backend" => "cpu",
            "steps" => sp.horizon * length(recs), "wall_seconds" => 1.0, "thread_seconds" => 2.0, "fallback_reason" => "", "first_attempt_reasons" => String[],
            "setup_check" => "not applicable", "control_stop" => false, "control_info" => Any[],
            "flags" => Dict{String, Any}("diagnostics exceeded" => Int[], "drift" => design.drift, "excluded" => Int[], "breakdown_tau" => Int[]))
        if sp.kind != ""
            ent["indicator_flags"] = Dict{String, Any}(k => v for (k, v) in p2_flag_counts(recs))
            sp.kind == "frequency" && (ent["transform_seconds"] = 1.5)
            sp.role in (:fq1, :bk1) && sp.arms[1] == "win_p" && (ent["not_applicable"] = Dict{String, Any}(sp.id => [20]))
        end
        push!(chunks, ent)
    end
    arms = Dict{String, Any}()
    for sp in specs, a in sp.arms
        e = get!(arms, a, Dict{String, Any}("status" => "complete", "steps" => 0, "cost_ns_per_orbit_step" => 100.0, "budget_ns" => 415.2, "specs" => Dict{String, Any}()))
        haskey(P2_COST, sp.id) && (e["specs"][sp.id] = Dict{String, Any}("budget_line_ns" => P2_COST[sp.id][1], "measured_ns_per_orbit_step" => P2_COST[sp.id][2]))
    end
    inds = Dict{String, Any}(id => Dict{String, Any}("kind" => P2_SPEC_INFO[id][1], "windows" => p.windows,
                                                    "bytes_per_orbit" => Dict("B" => 0, "M" => 80, "R" => 24, "S" => 9128, "total" => 9232 + 100 * i),
                                                    "budget" => Dict("Float64" => Dict("ns_cpu_thread" => P2_COST[id][1] / 2, "safety" => 2.0)))
                             for (i, id) in enumerate(sort(collect(keys(P2_SPEC_INFO)))))
    manifest = joinpath(dir, "$P2_SLUG-manifest.json")
    AX._write_json(manifest, Dict{String, Any}("study" => P2_SLUG, "schema_version" => 2, "status" => status, "plan_sha256" => p.plan_sha256,
        "plan_counts" => Dict("items" => total, "chunks" => length(chunks), "steps_before" => 1, "steps_after" => 1),
        "git" => Dict("sha" => "b"^40, "clean" => true), "sessions" => Any[Dict("clean" => true, "git_sha" => "b"^40, "index" => 1, "threads_used" => 4, "wall_seconds" => 10.0)],
        "prereg" => Dict("header" => h, "sha256_start" => h, "sha256_end" => h),
        "panel" => Dict("sha256_start" => ph, "sha256_end" => ph, "build_hash" => p.panel_build_sha256), "overrides" => Dict{String, Any}(), "control_stops" => Any[],
        "reasons" => String[], "chunks" => chunks, "arms" => arms, "indicators" => inds))
    return (; dir, store, prereg, panel, manifest, params, p, design, specs, chunks = length(chunks), items = total, study = P2_SLUG, prereg_hash = h, panel_hash = ph,
            out = joinpath(dir, "out"), h4, p1, variant)
end
const P2_STUDIES = Dict{String, Any}()
p2_cached(tag; kw...) = get!(() -> p2_study(tag; kw...), P2_STUDIES, tag)

# copy of a study (P2 and its P1 reference) into a new directory, with the paths rewritten
function p2_copy(st, tag::AbstractString)
    dir = an_path("p2copy-$tag")
    rm(dir; recursive = true, force = true)
    cp(st.dir, dir)
    rm(joinpath(dir, "out"); recursive = true, force = true)
    sub(s) = replace(s, st.dir => dir)
    p1 = merge(st.p1, (; dir = sub(st.p1.dir), store = sub(st.p1.store), prereg = sub(st.p1.prereg), panel = sub(st.p1.panel), manifest = sub(st.p1.manifest),
                       params = sub(st.p1.params)))
    return merge(st, (; dir, store = sub(st.store), prereg = sub(st.prereg), panel = sub(st.panel), manifest = sub(st.manifest), params = sub(st.params),
                      out = joinpath(dir, "out"), p1))
end

p2_run(st; params = st.params, out = st.out, store = st.store, kw...) = AN.analyze(; prereg = st.prereg, manifest = st.manifest, panel = st.panel, params, store, out,
    p1_prereg = st.p1.prereg, p1_manifest = st.p1.manifest, p1_params = st.p1.params, p1_store = st.p1.store, kw...)
p2_cli_args(st; extra = String[]) = vcat(["analyze", st.prereg, st.manifest, st.panel, st.params, "--store", st.store, "--out", st.out, "--p1-store", st.p1.store,
    "--p1-preregistration", st.p1.prereg, "--p1-manifest", st.p1.manifest, "--p1-parameters", st.p1.params], extra)
p2_chunk(st, i) = joinpath(st.store, "chunk-" * lpad(string(i), 5, '0') * ".h5")

# editing helpers: HDF5 attributes and values of a chunk file, the manifest
# the opened group or dataset is closed before the file, so that the library flushes the change to the disk (the byte-level checks of the read-only tests rely on it)
function p2_set_attr(file, path, name, val)
    AX.HDF5.h5open(file, "r+") do h
        g = path == "" ? h : h[path]
        AX.HDF5.delete_attribute(g, name)
        AX.HDF5.write_attribute(g, name, val)
        path == "" || close(g)
        flush(h)
    end
end
function p2_set_value(file, path, idx, val)
    AX.HDF5.h5open(file, "r+") do h
        d = h[path]
        d[idx...] = val
        close(d)
        flush(h)
    end
end
function p2_edit_manifest(st, f)
    m = AX.read_json(st.manifest)
    f(m)
    AX._write_json(st.manifest, m)
end
# a hash of every byte of every input and store file, also of the P1 reference (the read-only checks, I12); hash of a long byte vector samples a few elements only, the hash of a String covers all
p2_fingerprint(st) = [hash(String(read(f))) for f in vcat([st.manifest, st.prereg, st.panel, st.params, st.p1.manifest, st.p1.prereg, st.p1.params],
    sort(filter(f -> endswith(f, ".h5"), readdir(st.store; join = true))), sort(filter(f -> endswith(f, ".h5"), readdir(st.p1.store; join = true))))]

# the stored items of a chunk file read back by the real reader (not by the module): orbit id => OrbitRecord
p2_read(st, i) = Dict(o.id => o for o in AX.read_chunk(p2_chunk(st, i)).orbits)

# listing of the data tree (relative path, size, modification time of every file), the criterion of I12
function p2_datatree_hash()
    root = joinpath(@__DIR__, "..", "..", "data")
    isdir(root) || return an_sha256("no data tree")
    lines = String[]
    for (d, _, fs) in walkdir(root), f in fs
        p = joinpath(d, f)
        push!(lines, "$(relpath(p, root))|$(filesize(p))|$(mtime(p))")
    end
    return an_sha256(join(sort!(lines), "\n"))
end


p2_tcol(t, nm) = (i = findfirst(==(nm), t.columns); i === nothing ? error("table has no column $nm; columns: $(t.columns)") : i)
p2_rowsof(t, pairs...) = [r for r in t.rows if all(((nm, v),) -> r[p2_tcol(t, nm)] == v, pairs)]
