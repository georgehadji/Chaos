# Independent reader of the run store for the files of the indicator specs (Amendment D86, oracle O20, criterion D86s (V)). It extends the reader of the
# D53 oracle O13 (verify/experiments-runner-reader.jl). It uses HDF5.jl and SHA only: no Chaoslab code is loaded, and nothing of the module is
# assumed except the layout written in the spec (R15, E24, E25). Usage:
#   julia --project=verify verify/experiments-indicators-reader.jl <store dir> <out file>
# One line per attribute and dataset (fields separated by |), then the checks:
#   A | file | -       | name | type | value              root attributes
#   O | file | orbit   | name | type | value              orbit attributes
#   D | file | orbit   | path | eltype | file dims (C order) | julia dims | sha256 of the raw bytes
#   X | file | orbit   | path | name | value              dataset attributes (columns, rows, bits, range, threshold, ...)
#   P | file | member  | name | value                     attributes of provenance/<member>
#   Q | file | id      | eltype | file dims | julia dims | sha256 | orbits    pooled/<id>/counts
#   K | file | orbit   | path | sha256 of the cost entries without the two time values
#   C | file | PASS/FAIL/NOTE | message
using HDF5, SHA

const ROOT_ATTRS = ["schema_version", "study", "chunk", "chunk_key", "chunk_hash", "plan_sha256", "panel_sha256", "git_sha", "code_version",
                    "prereg_sha256", "config_sha256", "julia_version", "packages", "cpu", "threads", "backend_request", "backend", "device",
                    "backend_fallback_reason", "precision", "method", "horizon", "tau", "vectors", "coordinate_mode", "grid_sha256",
                    "status", "attempt", "setup_check", "first_attempt_reasons", "failure_reasons", "pilot", "arms", "renorm_interval_tau",
                    "deviation_vectors_k", "seed", "chunk_size", "gpu_chunks", "steps_per_launch", "launches_per_chunk", "orthonormalizer",
                    "implementation", "created_utc", "complete"]
const ORBIT_ATTRS = ["system", "params", "ic", "class_label", "seed", "item_hash", "design_class", "flags", "steps_recorded"]
const NEW_ROOT = ["kind", "spec", "spec_sha256", "chain", "members"]
const KINDS = ["tangent", "pair", "stretching", "frequency", "phase_advance", "birkhoff"]
const MEMBER_COLS = Dict("SALI" => ["SALI.sali", "SALI.L"], "FLI" => ["FLI.fli", "FLI.fli1"],
                         "MEGNO" => ["MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1", "MEGNO.Y2", "MEGNO.Ybar2", "MEGNO.Ybar20_2"], "RLI" => ["RLI.DL", "RLI.RLI"])
const STRETCH_SCALARS = ["handoffs", "sum_1", "sum_2", "comp_1", "comp_2", "sumsq_1", "sumsq_2", "compsq_1", "compsq_2", "mean_1", "mean_2", "std_1", "std_2",
                         "below_1", "below_2", "above_1", "above_2", "first_full_step"]
const FREQ_COLS = ["nu1", "nu2", "amp1", "amp2", "mirror1", "mirror2", "second1", "second2", "dist", "delta", "D"]
const FREQ_BITS = ["below_floor", "degenerate", "failed", "tie1", "tie2", "near_zero1", "near_zero2", "mirror1", "mirror2", "refined1", "refined2"]
const PA_COLS = ["nu1", "nu2", "delta", "D"]
const PA_BITS = ["below_floor", "undefined window"]
const PA_SCAL = ["handoffs", "near_pi", "complete"]
const BK_COLS = ["avg1", "avg2", "S", "B", "delta", "absdig", "reldig", "dig"]
const BK_BITS = ["below_floor", "zero", "reldig_undefined"]
const BK_SCAL = ["handoffs", "complete"]
const DIAG_NAMES = Set(["det", "symp", "orth", "sum", "pair"])

store, outfile = ARGS[1], ARGS[2]
out = open(outfile, "w")
nfail = Ref(0)
chk(file, ok, msg) = ok || (nfail[] += 1; println(out, "C|", file, "|FAIL|", msg))

attr_val(v) = replace(v isa AbstractArray ? join(string.(v), ";") : string(v), '\n' => "\\n")
rawsha(A) = bytes2hex(sha256(collect(reinterpret(UInt8, vec(A)))))
cdims(d) = Tuple(Int.(HDF5.API.h5s_get_simple_extent_dims(HDF5.dataspace(d))[1]))   # file (row-major) extents; the Julia shape is these reversed
dattrs(d) = Dict(nm => HDF5.read_attribute(d, nm) for nm in keys(HDF5.attributes(d)))

function dump_dataset(file, orbit, path, d)
    A = read(d)
    println(out, "D|", file, "|", orbit, "|", path, "|", eltype(A), "|", join(cdims(d), "x"), "|", join(size(A), "x"), "|", rawsha(A))
    at = dattrs(d)
    for nm in sort(collect(keys(at)))
        println(out, "X|", file, "|", orbit, "|", path, "|", nm, "|", attr_val(at[nm]))
    end
    return A, at
end

# value of `key = value` in the canonical spec text
spec_val(txt, key) = (m = match(Regex("(?m)^indicator\\." * key * " = (.*)\$"), txt); m === nothing ? nothing : String(m.captures[1]))

files = sort(readdir(store))
names_seen = Dict{String, Int}()     # every attribute and dataset name met, for the D86w scan
note(nm) = (names_seen[nm] = get(names_seen, nm, 0) + 1)
for fn in files
    if !occursin(r"^chunk-\d{5}\.h5$", fn)
        println(out, "C|", fn, "|NOTE|not a chunk file name (", fn, ")")
        endswith(fn, ".h5") && chk(fn, false, "a file with an .h5 name that is not chunk-NNNNN.h5")
        continue
    end
    path = joinpath(store, fn)
    h5open(path, "r") do f
        at = Dict(nm => HDF5.read_attribute(f, nm) for nm in keys(HDF5.attributes(f)))
        for nm in sort(collect(keys(at)))
            println(out, "A|", fn, "|-|", nm, "|", typeof(at[nm]), "|", attr_val(at[nm]))
            note("root@" * nm)
        end
        method, status = get(at, "method", ""), get(at, "status", "")
        for nm in ROOT_ATTRS
            chk(fn, haskey(at, nm), "root attribute $nm missing")
        end
        chk(fn, !haskey(at, "perturbation") || method == "L2", "perturbation present for method $method")
        get(at, "schema_version", 0) == 2 || chk(fn, false, "schema_version is $(get(at, "schema_version", nothing))")
        chk(fn, Bool(get(at, "complete", false)), "complete is not true in a file with a final name")
        chk(fn, !Bool(get(at, "pilot", true)), "pilot is true in a study file")
        chk(fn, status in ("ok", "failed"), "status $status")
        chk(fn, get(at, "attempt", 0) in (1, 2), "attempt $(get(at, "attempt", nothing))")
        chk(fn, get(at, "backend", "") in ("cpu", "gpu"), "backend $(get(at, "backend", nothing))")
        chk(fn, get(at, "coordinate_mode", "") in ("lifted", "wrapped"), "coordinate_mode")
        chk(fn, get(at, "chunk", -1) == parse(Int, fn[7:11]), "attribute chunk differs from the file name ordinal")
        chk(fn, status != "failed" || !isempty(get(at, "failure_reasons", String[])), "failed without failure_reasons")
        chk(fn, status != "ok" || isempty(get(at, "failure_reasons", String[])), "status ok with failure_reasons")
        chk(fn, !haskey(f, "grid"), "/grid present (unused by the runner)")
        wrapped = get(at, "coordinate_mode", "") == "wrapped"
        isind = haskey(at, "kind")
        kind = get(at, "kind", "")
        horizon = get(at, "horizon", 0)
        spec, nbins, runlen = "", 0, 0
        id = method
        if isind
            # ---- R15 root attributes
            for nm in NEW_ROOT
                chk(fn, haskey(at, nm), "root attribute $nm missing (indicator file)")
            end
            chk(fn, kind in KINDS, "kind $kind")
            spec = String(get(at, "spec", ""))
            chk(fn, bytes2hex(sha256(spec)) == get(at, "spec_sha256", ""), "SHA-256 of the spec text differs from spec_sha256")
            chk(fn, occursin("indicator.id = \"" * method * "\"", spec), "spec text does not name the id $method")
            chk(fn, occursin(Regex("(?m)^indicator\\.kind = \"" * kind * "\"\$"), spec), "spec text does not name the kind $kind")
            chk(fn, !occursin("budget", spec), "the spec text holds the budget table")
            ck = String(get(at, "chunk_key", ""))
            dirpart = haskey(at, "direction") ? String(at["direction"]) : "-"
            pre = "method=$method kind=$kind spec=$(String(get(at, "spec_sha256", ""))[1:min(16, end)]) precision=$(get(at, "precision", "")) horizon=$horizon tau=$(get(at, "tau", 0)) eps=- direction=$dirpart mode=$(get(at, "coordinate_mode", "")) cell="
            chk(fn, startswith(ck, pre), "chunk_key $ck does not start with the key rebuilt from the attributes")
            members = collect(get(at, "members", String[]))
            chain = collect(get(at, "chain", String[]))
            chk(fn, length(chain) == length(members) && all(i -> startswith(chain[i], members[i]), eachindex(members)), "chain $chain does not follow members $members")
            if kind in ("tangent", "pair")
                chk(fn, all(m -> m in keys(MEMBER_COLS), members) && length(unique(members)) == length(members), "members $members")
                forms = [occursin("replace=true", c) for c in chain]
                chk(fn, count(forms) == 1 && forms[end], "chain $chain: exactly one replacing member and it is last (I19)")
                kind == "tangent" && chk(fn, !("SALI" in members) || members[end] == "SALI", "SALI is not last in $members")
                kind == "pair" && chk(fn, Set(members) <= Set(["RLI", "MEGNO"]), "pair members $members")
            else
                chk(fn, members == [kind] && chain == [kind] || (kind == "stretching" && chain == ["stretching(replace=true)"]), "members/chain of a $kind file: $members / $chain")
            end
            vec_expected = kind in ("frequency", "phase_advance", "birkhoff") ? 1 : 2
            chk(fn, get(at, "vectors", 0) == vec_expected, "vectors $(get(at, "vectors", nothing)) for kind $kind (expected $vec_expected)")
            hasdir = kind in ("tangent", "stretching") || (kind == "pair" && spec_val(spec, "frame_source") == "\"direction\"")
            chk(fn, haskey(at, "direction") == hasdir, "direction attribute presence $(haskey(at, "direction")) for kind $kind (expected $hasdir)")
            if kind == "pair"
                chk(fn, haskey(at, "frame_source") && haskey(at, "direction_rule") && haskey(at, "separation"), "pair root attributes frame_source, direction_rule, separation")
                chk(fn, "\"" * String(get(at, "frame_source", "")) * "\"" == spec_val(spec, "frame_source"), "frame_source differs from the spec text")
            else
                chk(fn, !haskey(at, "frame_source") && !haskey(at, "direction_rule"), "frame_source or direction_rule in a $kind file")
            end
            if kind in ("frequency", "phase_advance", "birkhoff")
                w = collect(get(at, "windows", Int[]))
                runlen = Int(get(at, "run_length", 0))
                chk(fn, !isempty(w) && issorted(w; lt = <=), "windows $w")
                chk(fn, runlen == (kind == "phase_advance" ? 2maximum(w) + 1 : 2maximum(w)), "run_length $runlen for windows $w and kind $kind")
                chk(fn, spec_val(spec, "windows") == "[" * join(w, ", ") * "]", "windows differ from the spec text")
                chk(fn, horizon == runlen, "horizon $horizon differs from run_length $runlen at tau 1")
            else
                chk(fn, !haskey(at, "windows") && !haskey(at, "run_length"), "windows or run_length in a $kind file")
            end
            kind == "birkhoff" ? chk(fn, haskey(at, "observable_names"), "observable_names missing") : chk(fn, !haskey(at, "observable_names"), "observable_names in a $kind file")
            kind == "stretching" && (nbins = parse(Int, spec_val(spec, "bins")))
            # ---- provenance groups
            chk(fn, haskey(f, "provenance"), "provenance group missing")
            if haskey(f, "provenance")
                pg = f["provenance"]
                chk(fn, sort(collect(keys(pg))) == sort(kind in ("tangent", "pair") ? members : [kind == "stretching" ? "stretching" : kind]),
                    "provenance groups $(collect(keys(pg))) differ from the members")
                for m in keys(pg)
                    pa = Dict(nm => HDF5.read_attribute(pg[m], nm) for nm in keys(HDF5.attributes(pg[m])))
                    chk(fn, !isempty(pa), "provenance/$m holds no attribute")
                    for nm in sort(collect(keys(pa)))
                        println(out, "P|", fn, "|", m, "|", nm, "|", attr_val(pa[nm]))
                        note("provenance@" * nm)
                    end
                end
                # a chain whose members all publish the key `indicator`: every member's value is there (E23, D86e)
                if kind in ("tangent", "pair")
                    vals = [Dict(nm => HDF5.read_attribute(pg[m], nm) for nm in keys(HDF5.attributes(pg[m])))["indicator"] for m in members if haskey(HDF5.attributes(pg[m]), "indicator")]
                    chk(fn, length(vals) == length(members) && length(unique(vals)) == length(members), "provenance key `indicator` not kept per member: $vals")
                end
            end
            chk(fn, !haskey(f, "lyap"), "lyap group at the root")
            # ---- pooled counts (stretching with a counts subset)
            if haskey(f, "pooled")
                chk(fn, kind == "stretching" && collect(keys(f["pooled"])) == [id], "pooled groups $(collect(keys(f["pooled"]))) in a $kind file")
                for pid in keys(f["pooled"])
                    g = f["pooled"][pid]
                    chk(fn, haskey(g, "counts") && haskey(HDF5.attributes(g), "orbits"), "pooled/$pid lacks counts or the attribute orbits")
                    if haskey(g, "counts")
                        C = read(g["counts"])
                        println(out, "Q|", fn, "|", pid, "|", eltype(C), "|", join(cdims(g["counts"]), "x"), "|", join(size(C), "x"), "|", rawsha(C), "|", HDF5.read_attribute(g, "orbits"))
                        chk(fn, eltype(C) === Int64 && size(C) == (nbins + 2, 2), "pooled counts type/shape $(eltype(C)) $(size(C)), expected Int64 $(nbins + 2) x 2")
                        no = HDF5.read_attribute(g, "orbits")
                        handoffs = div(horizon, Int(get(at, "tau", 1)))
                        status == "ok" && chk(fn, all(c -> sum(C[:, c]) == no * handoffs, 1:2), "pooled count columns do not sum to orbits x handoffs ($no x $handoffs)")
                    end
                end
            else
                chk(fn, !(kind == "stretching" && !occursin("indicator.counts = \"none\"", spec)), "pooled group missing for a stretching spec with a counts subset")
            end
        else
            chk(fn, !haskey(f, "provenance"), "provenance group in an old-method file")
        end
        haskey(f, "orbits") || (chk(fn, false, "no /orbits group"); return)
        og = f["orbits"]
        for oname in sort(collect(keys(og)); by = s -> parse(Int, s))
            g = og[oname]
            oa = Dict(nm => HDF5.read_attribute(g, nm) for nm in keys(HDF5.attributes(g)))
            for nm in sort(collect(keys(oa)))
                println(out, "O|", fn, "|", oname, "|", nm, "|", typeof(oa[nm]), "|", attr_val(oa[nm]))
                note("orbit@" * nm)
            end
            for nm in ORBIT_ATTRS
                chk(fn, haskey(oa, nm), "orbit $oname attribute $nm missing")
            end
            flags = get(oa, "flags", String[])
            chk(fn, all(x -> x in ("diagnostics exceeded", "drift", "excluded", "breakdown_tau"), flags), "orbit $oname flags $flags")
            t, _ = dump_dataset(fn, oname, "t", g["t"])
            n = length(t)
            chk(fn, eltype(t) === Int64 && ndims(t) == 1, "orbit $oname t not Int64 vector")
            chk(fn, get(oa, "steps_recorded", -1) == n, "orbit $oname steps_recorded differs from the length of t")
            kk = count(>(0), t)
            chk(fn, issorted(t[1:kk]; lt = <=) && all(==(0), t[(kk + 1):end]) && (status != "ok" || (kk == n && (n == 0 || t[end] == horizon))), "orbit $oname t: strictly increasing to the horizon in an ok chunk")
            S, _ = dump_dataset(fn, oname, "state", g["state"])
            chk(fn, eltype(S) === Float64 && size(S) == (n, 2), "orbit $oname state shape/type $(size(S))")
            chk(fn, haskey(g, "winding") == wrapped, "orbit $oname winding presence ($wrapped)")
            if haskey(g, "winding")
                Wd, _ = dump_dataset(fn, oname, "winding", g["winding"])
                chk(fn, eltype(Wd) === Int64 && size(Wd) == (n, 2), "orbit $oname winding shape/type")
            end
            if !isind
                chk(fn, collect(keys(g["lyap"])) == [id], "orbit $oname lyap datasets $(collect(keys(g["lyap"])))")
                L, _ = dump_dataset(fn, oname, "lyap/" * id, g["lyap"][id])
                chk(fn, eltype(L) === Float64 && size(L) == (n, get(at, "vectors", 0) == 1 ? 1 : 2), "orbit $oname lyap shape/type $(size(L))")
            end
            if haskey(g, "diagnostics")
                for dn in sort(collect(keys(g["diagnostics"])))
                    chk(fn, dn in DIAG_NAMES, "orbit $oname diagnostics name $dn")
                    D, da = dump_dataset(fn, oname, "diagnostics/" * dn, g["diagnostics"][dn])
                    chk(fn, eltype(D) === Float64 && size(D) == (n,), "orbit $oname diagnostics/$dn shape/type")
                    status == "ok" && chk(fn, all(isfinite, D) && issorted(D), "orbit $oname diagnostics/$dn not finite and non-decreasing")
                    chk(fn, haskey(da, "threshold") && haskey(da, "first_exceedance_step"), "orbit $oname diagnostics/$dn attributes")
                    note("diagnostics@threshold")
                end
            end
            if isind
                chk(fn, sort(collect(keys(g["diagnostics"]))) == ["det", "symp"], "orbit $oname diagnostics of an indicator file: $(collect(keys(g["diagnostics"])))")
            else
                chk(fn, haskey(g, "lyap"), "orbit $oname lyap missing in an old-method file")
            end
            haskey(g, "cost") || chk(fn, false, "orbit $oname cost missing")
            C, _ = dump_dataset(fn, oname, "cost/" * id, g["cost"][id])
            println(out, "K|", fn, "|", oname, "|cost/", id, "|", rawsha(C[[1, 2, 5]]))
            chk(fn, eltype(C) === Float64 && length(C) == (kind == "frequency" ? 6 : 5), "orbit $oname cost shape $(size(C)) for kind $kind")
            status == "ok" && chk(fn, all(isfinite, C), "orbit $oname non-finite cost")
            isind || continue
            # ---- the datasets of the spec
            chk(fn, collect(keys(g["indicators"])) == [id], "orbit $oname indicators/ holds $(collect(keys(g["indicators"])))")
            ig = g["indicators"][id]
            ds = sort(collect(keys(ig)))
            if kind in ("tangent", "pair", "stretching")
                haskey(ig, "series") || (chk(fn, false, "orbit $oname indicators/$id/series missing"); continue)
                Sr, sa = dump_dataset(fn, oname, "indicators/$id/series", ig["series"])
                cols = kind == "stretching" ? ["stretching.D"] : reduce(vcat, [MEMBER_COLS[m] for m in members])
                chk(fn, eltype(Sr) === Float64 && size(Sr) == (n, length(cols)), "orbit $oname series type/shape $(eltype(Sr)) $(size(Sr)), expected Float64 ($n, $(length(cols)))")
                chk(fn, collect(get(sa, "columns", String[])) == cols, "orbit $oname series columns $(get(sa, "columns", nothing)), expected $cols")
                status == "ok" && chk(fn, all(isfinite, Sr), "orbit $oname non-finite value in series of an ok chunk (E31)")
                chk(fn, !haskey(ig, "final") && !haskey(ig, "flags"), "orbit $oname final/flags in a series kind")
                if kind == "stretching"
                    chk(fn, haskey(ig, "scalars"), "orbit $oname scalars missing")
                    if haskey(ig, "scalars")
                        Sc, ca = dump_dataset(fn, oname, "indicators/$id/scalars", ig["scalars"])
                        chk(fn, eltype(Sc) === Float64 && size(Sc) == (18,) && collect(get(ca, "columns", String[])) == STRETCH_SCALARS, "orbit $oname stretching scalars")
                        chk(fn, haskey(ca, "range") && length(ca["range"]) == 2 && ca["range"][1] < ca["range"][2], "orbit $oname scalars@range")
                        status == "ok" && chk(fn, Sc[1] == div(horizon, Int(get(at, "tau", 1))), "orbit $oname stretching handoffs $(Sc[1])")
                        if haskey(ig, "counts")
                            Ct, _ = dump_dataset(fn, oname, "indicators/$id/counts", ig["counts"])
                            chk(fn, eltype(Ct) === UInt32 && size(Ct) == (nbins + 2, 2), "orbit $oname counts type/shape $(eltype(Ct)) $(size(Ct))")
                            status == "ok" && chk(fn, all(c -> sum(Int.(Ct[:, c])) == Int(Sc[1]), 1:2), "orbit $oname count columns do not sum to the handoff count $(Sc[1])")
                        end
                    end
                else
                    chk(fn, !haskey(ig, "scalars") && !haskey(ig, "counts"), "orbit $oname scalars/counts in a $kind file")
                end
                if kind == "pair"
                    chk(fn, all(haskey(oa, k) for k in ("Dx", "Dx_r", "mode", "xi0")), "orbit $oname pair attributes Dx, Dx_r, mode, xi0")
                    if haskey(oa, "xi0")
                        chk(fn, length(oa["xi0"]) == 2 && eltype(oa["xi0"]) === Float64, "orbit $oname xi0 shape")
                        status == "ok" && chk(fn, abs(sum(abs2, oa["xi0"]) - 1) < 1e-14, "orbit $oname xi0 is not a unit vector")
                    end
                    status == "ok" && haskey(oa, "Dx_r") && chk(fn, 0.5 * oa["Dx"] <= oa["Dx_r"] <= 2 * oa["Dx"], "orbit $oname Dx_r outside [Dx/2, 2 Dx]")
                    chk(fn, get(oa, "mode", "") == get(at, "coordinate_mode", "?"), "orbit $oname pair mode differs from the root coordinate_mode")
                end
            else
                haskey(ig, "final") || (chk(fn, false, "orbit $oname indicators/$id/final missing"); continue)
                chk(fn, !haskey(ig, "series"), "orbit $oname series stored for the window kind $kind (I23: the series is never stored)")
                Fi, fa = dump_dataset(fn, oname, "indicators/$id/final", ig["final"])
                Fl, ga = dump_dataset(fn, oname, "indicators/$id/flags", ig["flags"])
                cols, bitn, scal = kind == "frequency" ? (FREQ_COLS, FREQ_BITS, nothing) : kind == "phase_advance" ? (PA_COLS, PA_BITS, PA_SCAL) : (BK_COLS, BK_BITS, BK_SCAL)
                w = collect(at["windows"])
                nrow = kind == "birkhoff" ? length(w) * length(at["observable_names"]) : length(w)
                chk(fn, eltype(Fi) === Float64 && size(Fi) == (nrow, length(cols)), "orbit $oname final type/shape $(eltype(Fi)) $(size(Fi)), expected Float64 ($nrow, $(length(cols)))")
                chk(fn, collect(get(fa, "columns", String[])) == cols, "orbit $oname final columns")
                rows = collect(get(fa, "rows", String[]))
                chk(fn, length(rows) == nrow && all(startswith.(rows, "W=")), "orbit $oname final rows $rows")
                chk(fn, kind == "birkhoff" || rows == ["W=$x" for x in w], "orbit $oname final rows $rows differ from the windows $w")
                chk(fn, eltype(Fl) === Int64 && size(Fl) == (nrow,) && collect(get(ga, "bits", String[])) == bitn, "orbit $oname flags type/shape/bits")
                if kind == "frequency"
                    chk(fn, !haskey(ig, "scalars"), "orbit $oname scalars stored for frequency (E25: none)")
                    chk(fn, length(Fi) * 8 + nrow * 8 == nrow * 96, "96 bytes per window")
                else
                    chk(fn, haskey(ig, "scalars"), "orbit $oname scalars missing")
                    if haskey(ig, "scalars")
                        Sc, ca = dump_dataset(fn, oname, "indicators/$id/scalars", ig["scalars"])
                        chk(fn, eltype(Sc) === Float64 && size(Sc) == (length(scal),) && collect(get(ca, "columns", String[])) == scal, "orbit $oname scalars type/shape/columns")
                        status == "ok" && chk(fn, Sc[1] == runlen && Sc[end] == 1.0, "orbit $oname scalars: handoffs $(Sc[1]) vs run_length $runlen, complete $(Sc[end])")
                    end
                end
                chk(fn, !haskey(ig, "counts"), "orbit $oname counts in a window kind")
                if status == "ok"
                    # E31: a NaN only where its flag says undefined, never an infinity
                    bit(name) = Int64(1) << (findfirst(==(name), bitn) - 1)
                    for r in 1:nrow, c in 1:length(cols)
                        x = Fi[r, c]
                        isinf(x) && chk(fn, false, "orbit $oname final[$r,$c] ($(cols[c])) is infinite")
                        isnan(x) || continue
                        # E31 as the module applies it: `reldig` is NaN-tolerated under `reldig_undefined` or `zero` (S = 0 sets `zero` only); the other columns under `zero`
                        okflag = kind == "frequency" ? (Fl[r] & bit("degenerate")) != 0 :
                                 kind == "phase_advance" ? (Fl[r] & bit("undefined window")) != 0 :
                                 cols[c] == "reldig" ? (Fl[r] & (bit("reldig_undefined") | bit("zero"))) != 0 : (Fl[r] & bit("zero")) != 0
                        cols[c] == "reldig" && (Fl[r] & bit("reldig_undefined")) == 0 && (Fl[r] & bit("zero")) != 0 &&
                            println(out, "C|", fn, "|NOTE|orbit $oname row $r: reldig is NaN with `zero` set and `reldig_undefined` clear (E31 wording says reldig_undefined)")
                        chk(fn, okflag, "orbit $oname final[$r,$c] ($(cols[c])) is NaN without its flag (E31)")
                    end
                    for r in 1:nrow
                        chk(fn, Fl[r] >= 0 && Fl[r] < (Int64(1) << length(bitn)), "orbit $oname flags[$r] = $(Fl[r]) outside the named bits")
                        kind == "frequency" && chk(fn, (Fl[r] & bit("failed")) == 0, "orbit $oname frequency status failed in an ok chunk")
                    end
                end
            end
            for d in ds
                note("indicators/" * d)
            end
        end
    end
    println(out, "C|", fn, "|DONE|read")
end
for (nm, c) in sort(collect(names_seen); by = first)
    println(out, "N|", nm, "|", c)
end
println(out, "C|store|", nfail[] == 0 ? "PASS" : "FAIL", "|", length(files), " entries, ", nfail[], " failed checks")
close(out)
println(nfail[] == 0 ? "READER PASS" : "READER FAIL $(nfail[])")
exit(nfail[] == 0 ? 0 : 1)
