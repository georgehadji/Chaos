# Shared set-up of the D86 verify scripts (included, not run): paths, the reduced preregistration copy, the reduced panel copy, kill and fault hooks,
# the independent reader driver. Nothing here is written outside the scratch folder EXI_SCRATCH. The study data folder `data/` is never written (D86y).
using Printf, Random, LinearAlgebra, Statistics, SHA, Dates, TOML
using Chaoslab, StaticArrays

using DoubleFloats
const EX = Chaoslab.Experiments
# FINDING F1 (stage 04, D86): the provenance attribute writer `_h5val` of src/Experiments/run_indicator_store.jl (line 86) lets a Double64 pass as an
# AbstractFloat and HDF5.jl has no method to write one, so every Double64 chunk of a `stretching` spec raises a MethodError (spectrum_lo, spectrum_hi,
# spectrum_cell_width of the provenance). The campaign scripts therefore define, in this verification script only and never in src/, the method below
# (the value as its decimal text, which is what `_h5val` already does for every other type). EXI_PATCH=0 runs without it (section F1 of
# verify/experiments-indicators.jl shows the unpatched failure).
if get(ENV, "EXI_PATCH", "1") == "1"
    Chaoslab.Experiments._h5val(v::DoubleFloats.Double64) = string(v)
end
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const PREREG = joinpath(ROOT, "workspaces", "study", "stages", "01-design", "output", "p1-le-benchmark-preregistration.md")
const PANEL = joinpath(ROOT, "workspaces", "study", "stages", "02-panel", "output", "p1-le-benchmark-panel.toml")
const BLOCK = joinpath(@__DIR__, "experiments-indicators-block.toml")
const SCR = get(ENV, "EXI_SCRATCH", raw"C:\Users\tesse\AppData\Local\Temp\claude\C--Users-tesse-Documents-Vibe-Coding-Chaos\5df15ee1-47df-4875-85bd-1bd4d7c665a6\scratchpad\d86-verify")
const VERIFY_ENV = joinpath(ROOT, "verify")
const NT = Threads.nthreads()
const PER_K = 60
mkpath(SCR)

lf(s) = replace(s, "\r\n" => "\n")
bits(A) = collect(reinterpret(UInt8, vec(A)))
rawsha(A) = bytes2hex(sha256(bits(A)))
sha_file(p) = bytes2hex(open(sha256, p))
function filecount(dir)
    isdir(dir) || return 0
    return sum((length(fs) for (_, _, fs) in walkdir(dir)); init = 0)
end
function tree_listing(dir)
    isdir(dir) || return String[]
    return sort!([relpath(joinpath(r, f), dir) * "|" * string(filesize(joinpath(r, f))) * "|" * string(mtime(joinpath(r, f))) for (r, _, fs) in walkdir(dir) for f in fs])
end

# ---------------------------------------------------------------- reduced panel: PER_K drawn orbits per K (round robin over the classes) plus all controls;
# in_SL is true for the first selected orbit of each cell only (the grammar offers no per-cell selection; S_L is the one-orbit-per-cell subset of the check)
const CLASS_ORDER = ["T", "R", "C", "S", "X"]
function reduced_panel(path; per_k = PER_K)
    t = TOML.parsefile(PANEL)
    orb = t["orbit"]
    keep = Dict{String, Any}[]
    for k in sort(unique(o["k_index"] for o in orb if !o["control"]))
        cells = Dict(c => sort!([o for o in orb if !o["control"] && o["k_index"] == k && o["class"] == c]; by = o -> o["draw_index"]) for c in CLASS_ORDER)
        cells = Dict(c => v for (c, v) in cells if !isempty(v))
        pos = Dict(c => 0 for c in keys(cells))
        taken = 0
        while taken < per_k && any(c -> pos[c] < length(cells[c]), keys(cells))
            for c in CLASS_ORDER
                haskey(cells, c) && pos[c] < length(cells[c]) && taken < per_k || continue
                pos[c] += 1
                o = copy(cells[c][pos[c]])
                o["in_SL"] = pos[c] == 1
                push!(keep, o)
                taken += 1
            end
        end
    end
    append!(keep, [copy(o) for o in orb if o["control"]])
    sort!(keep; by = o -> o["id"])
    t2 = copy(t)
    t2["orbit"] = keep
    open(io -> TOML.print(io, t2; sorted = true), path, "w")
    return path, (selected = length(keep) - count(o -> o["control"], keep), controls = count(o -> o["control"], keep), sl = count(o -> o["in_SL"], keep))
end

# ---------------------------------------------------------------- reduced preregistration: the real file with its first [run] block replaced by BLOCK
function run_toml_blocks(text)
    out = []
    for m in eachmatch(r"```toml\n(.*?)\n```"s, text)
        d = try TOML.parse(String(m.captures[1])) catch; continue end
        haskey(d, "run") && push!(out, (m, d))
    end
    isempty(out) && error("no [run] block")
    return out
end
function make_prereg(path, run; real = PREREG)
    text = lf(read(real, String))
    new = text
    for (i, (m, d)) in enumerate(run_toml_blocks(text))
        if i == 1
            io = IOBuffer(); TOML.print(io, Dict("run" => run); sorted = true)
            blk = "```toml\n" * rstrip(String(take!(io))) * "\n```"
        else
            blk = "```toml\n# emptied in the verification copy\n```"
        end
        new = replace(new, String(m.match) => blk; count = 1)
    end
    h = EX.prereg_hash(new)
    new = replace(new, r"(?m)^(\| SHA-256 of this file[^|]*\|[ \t]*)[0-9a-f]{64}" => SubstitutionString("\\g<1>" * h); count = 1)
    write(path, new)
    EX.check_prereg(path)
    return path
end
block_run() = TOML.parsefile(BLOCK)["run"]

function prepare(dir = SCR)
    mkpath(dir)
    pan, info = reduced_panel(joinpath(dir, "panel-reduced.toml"))
    pre = make_prereg(joinpath(dir, "prereg-reduced.md"), block_run())
    return (prereg = pre, panel = pan, info = info)
end

# ---------------------------------------------------------------- kill and fault hooks (as verify/experiments-runner.jl)
struct VKill <: Exception
    point::Symbol
    n::Int
end
has_kill(::VKill) = true
has_kill(e::CompositeException) = any(has_kill, e.exceptions)
has_kill(e::TaskFailedException) = has_kill(e.task.result)
has_kill(e::LoadError) = has_kill(e.error)
has_kill(::Any) = false
function killed(f)
    try
        f()
    catch e
        has_kill(e) && return true
        rethrow()
    end
    return false
end
kill_after(k) = (p, info) -> (p === :chunk_done && info.chunks_done == k && throw(VKill(p, k)); nothing)
kill_during(ord) = (p, info) -> (p === :during_data && info.ordinal == ord && throw(VKill(p, ord)); nothing)
# a NaN in the first column of the series of the second orbit of chunk `ord`, on the listed attempts (a failure of the chunk, R18 item 1)
function nan_fault(ord, id; attempts)
    return (info, orbits) -> begin
        if info.ordinal == ord && info.attempt in attempts
            t = orbits[min(2, end)].indicators[id]
            t.series[min(3, end), 1] = NaN
        end
        nothing
    end
end
chunkname(ord) = @sprintf("chunk-%05d.h5", ord)

# ---------------------------------------------------------------- independent reader (separate process, HDF5.jl and SHA only, no Chaoslab)
const READER = joinpath(@__DIR__, "experiments-indicators-reader.jl")
function run_reader(store, tag)
    out = joinpath(SCR, "reader-$tag.txt")
    t = @elapsed ok = success(pipeline(`$(Base.julia_cmd()) --project=$VERIFY_ENV -t 1 $READER $store $out`; stdout = joinpath(SCR, "reader-$tag.out"), stderr = joinpath(SCR, "reader-$tag.err")))
    return (path = out, ok = ok, seconds = t)
end
function parse_reader(path)
    R = (attrs = Dict{Tuple{String, String}, String}(), oattrs = Dict{Tuple{String, String, String}, String}(),
         data = Dict{Tuple{String, String, String}, NamedTuple}(), dattrs = Dict{Tuple{String, String, String, String}, String}(),
         pattrs = Dict{Tuple{String, String, String}, String}(), pooled = Dict{Tuple{String, String}, NamedTuple}(),
         ksha = Dict{Tuple{String, String, String}, String}(), checks = String[], files = Set{String}())
    for line in eachline(path)
        isempty(line) && continue
        k = line[1]
        if k == 'A'
            f = split(line, '|'; limit = 6); R.attrs[(f[2], f[4])] = f[6]; push!(R.files, f[2])
        elseif k == 'O'
            f = split(line, '|'; limit = 6); R.oattrs[(f[2], f[3], f[4])] = f[6]
        elseif k == 'D'
            f = split(line, '|'); R.data[(f[2], f[3], f[4])] = (eltype = f[5], cdims = f[6], jdims = f[7], sha = f[8])
        elseif k == 'X'
            f = split(line, '|'; limit = 6); R.dattrs[(f[2], f[3], f[4], f[5])] = f[6]
        elseif k == 'P'
            f = split(line, '|'; limit = 5); R.pattrs[(f[2], f[3], f[4])] = f[5]
        elseif k == 'Q'
            f = split(line, '|'); R.pooled[(f[2], f[3])] = (eltype = f[4], cdims = f[5], jdims = f[6], sha = f[7], orbits = f[8])
        elseif k == 'K'
            f = split(line, '|'); R.ksha[(f[2], f[3], f[4])] = f[5]
        elseif k == 'C'
            push!(R.checks, line)
        end
    end
    return R
end
attr_str(v) = replace(v isa AbstractArray ? join(string.(v), ";") : string(v), '\n' => "\\n")
