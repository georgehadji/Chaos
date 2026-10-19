# Independent reader of the run store (D53m, oracle O13). Uses HDF5.jl and SHA only: no Chaoslab code is loaded.
# Usage: julia --project=verify verify/experiments-runner-reader.jl <store dir> <out file>
# Reads every file of the store, checks the version 2 layout of schema/run-store.md and writes one line per attribute and dataset:
#   A | file | -         | name | type | value
#   O | file | <orbit>   | name | type | value            (orbit attributes)
#   D | file | <orbit>   | path | eltype | file dims (C order) | julia dims | sha256 of the raw bytes
#   X | file | <orbit>   | path | name | value             (dataset attributes: threshold, first_exceedance_step)
#   K | file | <orbit>   | path | sha256 of cost entries 1, 2, 5 (steps, vectors, resident bytes; the two time entries left out)
#   C | file | PASS/FAIL | message                        (one line per failed check, and per file a PASS line)
using HDF5, SHA

const ROOT_ATTRS = ["schema_version", "study", "chunk", "chunk_key", "chunk_hash", "plan_sha256", "panel_sha256", "git_sha", "code_version",
                    "prereg_sha256", "config_sha256", "julia_version", "packages", "cpu", "threads", "backend_request", "backend", "device",
                    "backend_fallback_reason", "precision", "method", "horizon", "tau", "vectors", "coordinate_mode", "direction", "grid_sha256",
                    "status", "attempt", "setup_check", "first_attempt_reasons", "failure_reasons", "pilot", "arms", "renorm_interval_tau",
                    "deviation_vectors_k", "seed", "chunk_size", "gpu_chunks", "steps_per_launch", "launches_per_chunk", "orthonormalizer",
                    "implementation", "created_utc", "complete"]
const ORBIT_ATTRS = ["system", "params", "ic", "class_label", "seed", "item_hash", "design_class", "flags", "steps_recorded"]
const DIAG_NAMES = Set(["det", "symp", "orth", "sum", "pair"])

store, outfile = ARGS[1], ARGS[2]
out = open(outfile, "w")
nfail = Ref(0)
chk(file, ok, msg) = ok || (nfail[] += 1; println(out, "C|", file, "|FAIL|", msg))

attr_val(v) = replace(v isa AbstractArray ? join(string.(v), ";") : string(v), '\n' => "\\n")
rawsha(A) = bytes2hex(sha256(collect(reinterpret(UInt8, vec(A)))))
cdims(d) = Tuple(Int.(HDF5.API.h5s_get_simple_extent_dims(HDF5.dataspace(d))[1]))   # file (row-major) extents, not reversed

function dump_dataset(file, orbit, path, d)
    A = read(d)
    println(out, "D|", file, "|", orbit, "|", path, "|", eltype(A), "|", join(cdims(d), "x"), "|", join(size(A), "x"), "|", rawsha(A))
    for nm in keys(HDF5.attributes(d))
        println(out, "X|", file, "|", orbit, "|", path, "|", nm, "|", attr_val(HDF5.read_attribute(d, nm)))
    end
    return A
end

files = sort(readdir(store))
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
            v = at[nm]
            println(out, "A|", fn, "|-|", nm, "|", typeof(v), "|", attr_val(v))
        end
        method = get(at, "method", "")
        for nm in ROOT_ATTRS
            chk(fn, haskey(at, nm), "root attribute $nm missing")
        end
        chk(fn, !haskey(at, "perturbation") || method == "L2", "perturbation present for method $method")
        chk(fn, method != "L2" || haskey(at, "perturbation"), "perturbation missing for L2")
        get(at, "schema_version", 0) == 2 || chk(fn, false, "schema_version is $(get(at, "schema_version", nothing))")
        chk(fn, Bool(get(at, "complete", false)), "complete is not true in a file with a final name")
        chk(fn, !Bool(get(at, "pilot", true)), "pilot is true in a study file")
        chk(fn, get(at, "status", "") in ("ok", "failed"), "status $(get(at, "status", nothing))")
        chk(fn, get(at, "attempt", 0) in (1, 2), "attempt $(get(at, "attempt", nothing))")
        chk(fn, get(at, "backend", "") in ("cpu", "gpu"), "backend $(get(at, "backend", nothing))")
        chk(fn, get(at, "coordinate_mode", "") in ("lifted", "wrapped"), "coordinate_mode")
        chk(fn, get(at, "chunk", -1) == parse(Int, fn[7:11]), "attribute chunk differs from the file name ordinal")
        chk(fn, get(at, "status", "") != "failed" || !isempty(get(at, "failure_reasons", String[])), "failed without failure_reasons")
        chk(fn, get(at, "attempt", 1) != 2 || get(at, "status", "") != "ok" || !isempty(get(at, "first_attempt_reasons", String[])),
            "attempt 2 with status ok but no first_attempt_reasons")
        chk(fn, get(at, "attempt", 1) != 1 || isempty(get(at, "first_attempt_reasons", String[])), "attempt 1 with first_attempt_reasons")
        chk(fn, get(at, "status", "") != "ok" || isempty(get(at, "failure_reasons", String[])), "status ok with failure_reasons")
        chk(fn, !haskey(f, "grid"), "/grid present (unused by the runner)")
        k = get(at, "vectors", 0)
        k == (method in ("L1", "L2") ? 1 : 2) || chk(fn, false, "vectors $k for method $method")
        horizon, wrapped, status = get(at, "horizon", 0), get(at, "coordinate_mode", "") == "wrapped", get(at, "status", "")
        haskey(f, "orbits") || (chk(fn, false, "no /orbits group"); return)
        og = f["orbits"]
        for oname in sort(collect(keys(og)); by = s -> parse(Int, s))
            g = og[oname]
            oa = Dict(nm => HDF5.read_attribute(g, nm) for nm in keys(HDF5.attributes(g)))
            for nm in sort(collect(keys(oa)))
                println(out, "O|", fn, "|", oname, "|", nm, "|", typeof(oa[nm]), "|", attr_val(oa[nm]))
            end
            for nm in ORBIT_ATTRS
                chk(fn, haskey(oa, nm), "orbit $oname attribute $nm missing")
            end
            flags = get(oa, "flags", String[])
            chk(fn, all(x -> x in ("diagnostics exceeded", "drift", "excluded"), flags), "orbit $oname flags $flags")
            t = dump_dataset(fn, oname, "t", g["t"])
            n = length(t)
            chk(fn, eltype(t) === Int64 && ndims(t) == 1, "orbit $oname t not Int64 vector")
            chk(fn, get(oa, "steps_recorded", -1) == n, "orbit $oname steps_recorded differs from the length of t")
            kk = count(>(0), t)       # a failed chunk blanks the rows after the failure point (t = 0, values NaN)
            chk(fn, issorted(t[1:kk]; lt = <=) && all(==(0), t[(kk + 1):end]) && (status != "ok" || (kk == n && (n == 0 || t[end] == horizon))), "orbit $oname t: strictly increasing to the horizon in an ok chunk, blank tail only in a failed one")
            S = dump_dataset(fn, oname, "state", g["state"])
            chk(fn, eltype(S) === Float64 && size(S) == (n, 2), "orbit $oname state shape/type")
            chk(fn, haskey(g, "winding") == wrapped, "orbit $oname winding presence ($wrapped)")
            if haskey(g, "winding")
                W = dump_dataset(fn, oname, "winding", g["winding"])
                chk(fn, eltype(W) === Int64 && size(W) == (n, 2), "orbit $oname winding shape/type")
            end
            chk(fn, collect(keys(g["lyap"])) == [method], "orbit $oname lyap datasets $(collect(keys(g["lyap"])))")
            L = dump_dataset(fn, oname, "lyap/" * method, g["lyap"][method])
            chk(fn, eltype(L) === Float64 && size(L) == (n, k), "orbit $oname lyap shape/type $(size(L))")
            if haskey(g, "diagnostics")
                for dn in sort(collect(keys(g["diagnostics"])))
                    chk(fn, dn in DIAG_NAMES, "orbit $oname diagnostics name $dn")
                    D = dump_dataset(fn, oname, "diagnostics/" * dn, g["diagnostics"][dn])
                    chk(fn, eltype(D) === Float64 && size(D) == (n,), "orbit $oname diagnostics/$dn shape/type")
                    status == "ok" && chk(fn, all(isfinite, D) && issorted(D), "orbit $oname diagnostics/$dn not finite and non-decreasing")
                    da = g["diagnostics"][dn]
                    chk(fn, haskey(HDF5.attributes(da), "threshold") && haskey(HDF5.attributes(da), "first_exceedance_step"), "orbit $oname diagnostics/$dn attributes")
                end
            end
            chk(fn, (method == "L2") == !haskey(g, "diagnostics"), "orbit $oname diagnostics presence for $method")
            C = dump_dataset(fn, oname, "cost/" * method, g["cost"][method])
            println(out, "K|", fn, "|", oname, "|cost/", method, "|", rawsha(C[[1, 2, 5]]))      # cost entries without the two time values
            chk(fn, eltype(C) === Float64 && length(C) == 5, "orbit $oname cost shape/type")
            if status == "ok"
                chk(fn, all(isfinite, S) && all(isfinite, L) && all(isfinite, C), "orbit $oname non-finite value in an ok chunk")
            end
        end
    end
    println(out, "C|", fn, "|DONE|read")
end
println(out, "C|store|", nfail[] == 0 ? "PASS" : "FAIL", "|", length(files), " entries, ", nfail[], " failed checks")
close(out)
println(nfail[] == 0 ? "READER PASS" : "READER FAIL $(nfail[])")
exit(nfail[] == 0 ? 0 : 1)
