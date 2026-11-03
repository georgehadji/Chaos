# Independent reader of a test-run store of study p3-abstention (oracle O20 reader, Amendment D95 of the analysis spec).
# Uses HDF5, SHA, TOML only: no Chaoslab code is loaded; no file of src/, test/, bench/, ext/ was read. Written from schema/run-store.md,
# analysis-spec.md (D95 E59 to E65, D89 E23 to E26) and experiments-spec.md (D94 R23, E39 to E43), and from verify/analysis-p3-export.jl (same author line).
#
# Usage (from the repository root):
#   julia --project=verify verify/analysis-p3t-export.jl <store folder> <output folder> [<N_1> <N_c> <n_s>]
# Output: <output folder>/p3t-export.jsonl, one JSON object per line (NaN and +-Infinity written as such):
#   {"k":"chunk", "fn", "sha256" (file), "chunk_key", "chunk_hash", "arms", "status", "complete", "study", "panel_sha256", "prereg_sha256", "method",
#    "precision", "horizon", "direction", "mode", "kind", "orbits":[ids], "attempt"}
#   {"k":"it", "fn", "o", "m" (method or spec id), "kind", "prec", "hor", "dir", "mode", "status", "flags":[orbit flags], "ic", "params", "design_class",
#    "t":[recorded steps <= TRUNC], and by kind:
#      lyap (L1, L5): "lam":[first column up to TRUNC], "lamfull": {"<step>": value} for every recorded step (all horizons), "xl": {"<step>": lifted first coordinate} (L5)
#      tangent: "c":{column: values up to TRUNC}
#      window kinds: "rows", "c":{D | dig}, "fw":[flag word per row], "bits":[names]}
using HDF5, SHA, TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
length(ARGS) >= 2 || error("usage: julia --project=verify verify/analysis-p3t-export.jl <store folder> <output folder>")
const STORE = ARGS[1]
const OUTD = ARGS[2]
occursin(r"/data(/|$)", replace(abspath(OUTD), "\\" => "/")) && error("output folder inside data/")
mkpath(OUTD)
const TRUNC = 200_000

jnum(x::Real) = x isa Integer ? string(x) : (isnan(x) ? "NaN" : (isinf(x) ? (x > 0 ? "Infinity" : "-Infinity") : repr(Float64(x))))
jstr(s::AbstractString) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n") * "\""
function jval(v)
    v === nothing && return "null"
    v isa Bool && return string(v)
    v isa AbstractString && return jstr(v)
    v isa Real && return jnum(v)
    v isa AbstractVector && return "[" * join((jval(x) for x in v), ",") * "]"
    v isa AbstractDict && return "{" * join((jstr(string(k)) * ":" * jval(x) for (k, x) in sort(collect(v); by = first)), ",") * "}"
    error("jval: $(typeof(v))")
end

chunk_attrs(f) = Dict(nm => HDF5.read_attribute(f, nm) for nm in keys(HDF5.attributes(f)))
strs(x) = x isa AbstractString ? [String(x)] : String.(collect(x))
sattr(a, k) = haskey(a, k) ? String(a[k]) : ""

files = sort(filter(fn -> occursin(r"^chunk-\d{5}\.h5$", fn), readdir(STORE)))
println("reading ", length(files), " chunk files of ", STORE)
nit = 0
open(joinpath(OUTD, "p3t-export.jsonl"), "w") do io
    for fn in files
        path = joinpath(STORE, fn)
        fsha = bytes2hex(sha256(read(path)))
        h5open(path, "r") do f
            a = chunk_attrs(f)
            prec = startswith(sattr(a, "precision"), "BigFloat") ? "BigFloat" : sattr(a, "precision")
            kind = sattr(a, "kind")
            m = sattr(a, "method")
            hor = haskey(a, "horizon") ? Int(a["horizon"]) : 0
            dir = sattr(a, "direction")
            mode = sattr(a, "coordinate_mode")
            ids = haskey(f, "orbits") ? sort(parse.(Int, collect(keys(f["orbits"])))) : Int[]
            println(io, jval(Dict("k" => "chunk", "fn" => fn, "sha256" => fsha, "chunk_key" => sattr(a, "chunk_key"), "chunk_hash" => sattr(a, "chunk_hash"),
                                  "arms" => haskey(a, "arms") ? strs(a["arms"]) : String[], "status" => sattr(a, "status"), "complete" => haskey(a, "complete") ? Bool(a["complete"]) : false,
                                  "study" => sattr(a, "study"), "panel_sha256" => sattr(a, "panel_sha256"), "prereg_sha256" => sattr(a, "prereg_sha256"), "method" => m, "precision" => prec,
                                  "horizon" => hor, "direction" => dir, "mode" => mode, "kind" => kind, "orbits" => ids, "attempt" => haskey(a, "attempt") ? Int(a["attempt"]) : 0)))
            for id in ids
                g = f["orbits"][string(id)]
                d = Dict{String,Any}("k" => "it", "fn" => fn, "o" => id, "m" => m, "kind" => kind, "prec" => prec, "hor" => hor, "dir" => dir, "mode" => mode, "status" => sattr(a, "status"),
                                     "flags" => (haskey(HDF5.attributes(g), "flags") ? strs(HDF5.read_attribute(g, "flags")) : String[]),
                                     "ic" => (haskey(HDF5.attributes(g), "ic") ? String(HDF5.read_attribute(g, "ic")) : ""),
                                     "params" => (haskey(HDF5.attributes(g), "params") ? String(HDF5.read_attribute(g, "params")) : ""),
                                     "design_class" => (haskey(HDF5.attributes(g), "design_class") ? String(HDF5.read_attribute(g, "design_class")) : ""))
                t = collect(Int64, read(g["t"]))
                tr = findlast(<=(TRUNC), t); tr === nothing && (tr = 0)
                d["t"] = t[1:tr]
                try
                    if kind == "" && haskey(g, "lyap") && (m == "L1" || m == "L5")
                        L = read(g["lyap"][m]); lam = ndims(L) == 1 ? collect(L) : collect(L[:, 1])   # file extent (k, n): julia (n, k) after the schema note; HDF5.jl returns (k,n)->see below
                        if ndims(L) == 2
                            # HDF5.jl returns the Julia shape (n, k) as written by Julia; guard for the (k, n) order
                            lam = size(L, 1) == length(t) ? collect(L[:, 1]) : collect(L[1, :])
                        end
                        d["lam"] = lam[1:tr]
                        d["lamfull"] = Dict(string(t[i]) => lam[i] for i in eachindex(t))
                        if m == "L5"
                            S = read(g["state"]); S1 = size(S, 1) == length(t) ? collect(S[:, 1]) : collect(S[1, :])
                            if mode == "wrapped"
                                W = read(g["winding"]); W1 = size(W, 1) == length(t) ? collect(W[:, 1]) : collect(W[1, :])
                                S1 = S1 .+ 2pi .* Float64.(W1)
                            end
                            d["xl"] = Dict(string(t[i]) => S1[i] for i in eachindex(t))
                        end
                    elseif kind == "tangent" && haskey(g, "indicators")
                        ds = g["indicators"][m]["series"]
                        cols = String.(HDF5.read_attribute(ds, "columns")); M = read(ds)
                        M = size(M, 1) == length(t) ? M : permutedims(M)
                        c = Dict{String,Vector{Float64}}()
                        for (j, cn) in enumerate(cols)
                            c[cn] = collect(M[1:tr, j])
                        end
                        d["c"] = c
                    elseif kind in ("frequency", "phase_advance", "birkhoff") && haskey(g, "indicators")
                        grp = g["indicators"][m]
                        ds = grp["final"]; cols = String.(HDF5.read_attribute(ds, "columns")); rows = String.(HDF5.read_attribute(ds, "rows"))
                        M = read(ds)
                        M = size(M, 1) == length(rows) ? M : permutedims(M)
                        c = Dict{String,Vector{Float64}}()
                        for cn in ("D", "dig")
                            j = findfirst(==(cn), cols); j === nothing || (c[cn] = collect(M[:, j]))
                        end
                        fds = grp["flags"]
                        d["rows"] = rows; d["c"] = c; d["fw"] = collect(Int64, read(fds)); d["bits"] = String.(HDF5.read_attribute(fds, "bits"))
                    end
                catch e
                    d["read_error"] = sprint(showerror, e)
                end
                println(io, jval(d))
                global nit += 1
            end
        end
    end
end
println("items exported: ", nit)
