# Golden capture of D125a (oracle O28): definitions only, shared by the generating script and the test, so both run the same code.
# Reads the tracked inputs of the three existing studies under workspaces/study/ (preregistrations, panels, round log); never data/.
using SHA, Random

const D125_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const D125_STUDY = Dict(   # study => (preregistration, panel, round log or nothing)
    "p1-le-benchmark" => ("p1-le-benchmark-preregistration.md", "p1-le-benchmark-panel.toml", nothing),
    "p2-indicators" => ("p2-indicators-preregistration.md", "p1-le-benchmark-panel.toml", nothing),     # p2 runs on the p1 panel (same panel SHA-256 in both manifests)
    "p3-abstention" => ("p3-abstention-preregistration.md", "p3-abstention-panel.toml", "p3-abstention-round-log.txt"))   # the tracked log names the code version of its day: see d125_synth_log
d125_study_paths(n) = (joinpath(D125_ROOT, "workspaces", "study", "stages", "01-design", "output", D125_STUDY[n][1]),
                       joinpath(D125_ROOT, "workspaces", "study", "stages", "02-panel", "output", D125_STUDY[n][2]),
                       D125_STUDY[n][3] === nothing ? nothing : joinpath(D125_ROOT, "workspaces", "study", "stages", "03-run", "output", D125_STUDY[n][3]))
const D125_CODE = "d125-golden-code"          # fixed code version string: the plan must not depend on the source
d125_sha(s::AbstractString) = bytes2hex(sha256(s))
d125_sha(v::AbstractVector{<:AbstractString}) = bytes2hex(sha256(join(v, "\n")))

# canonical text of a nested table (D96): keys sorted, one `path = repr(value)` line per leaf
function d125_canon_table(x, path = "")
    io = IOBuffer()
    if x isa AbstractDict
        for k in sort!(collect(string.(keys(x))))
            write(io, d125_canon_table(x[k], path * "/" * k))
        end
    elseif x isa AbstractVector
        for (i, v) in enumerate(x)
            write(io, d125_canon_table(v, path * "[$i]"))
        end
    else
        println(io, path, " = ", repr(x))
    end
    return String(take!(io))
end

# the values of one study: file hashes, plan hash, counts, SHA-256 of the lists of item hashes, chunk keys and chunk hashes, of the plan text, of the
# item texts of a seeded sample of min(10^4, n) items, the block hash, and (specs) the indicator manifest table
# `log`: path of a round log valid for the fixed code version (p3-abstention only; see d125_synth_log); the tracked log of the run records
# the digest of the code version of its day, which no later source can reproduce, so the golden holds a log regenerated for D125_CODE.
function d125_golden_capture(study::AbstractString; plan = nothing, log = nothing)
    pre, pan, _ = d125_study_paths(study)
    plan === nothing && (plan = log === nothing ? EX.plan_run(pre, pan; code_version = D125_CODE) :
                                                 EX.plan_run(pre, pan; code_version = D125_CODE, round_log = log))
    n = length(plan.items)
    idx = n <= 10_000 ? collect(1:n) : sort!(randperm(StableRNG(EX_SEED + 63050), n)[1:10_000])
    texts = [EX.item_text(plan.items[i]) for i in idx]
    ind = isempty(plan.not_applicable) && !any(it -> it.kind !== nothing, plan.items) ? "none" : d125_sha(d125_canon_table(EX.indicator_manifest(plan)))
    return (; prereg_sha = bytes2hex(sha256(read(pre))), panel_sha = bytes2hex(sha256(read(pan))), plan_hash = plan.hash, items = n,
            chunks = length(plan.chunks), item_hashes = d125_sha([it.hash for it in plan.items]), chunk_keys = d125_sha([c.key for c in plan.chunks]),
            chunk_hashes = d125_sha([c.hash for c in plan.chunks]), plan_text = d125_sha(plan.text), sample = length(idx), sample_texts = d125_sha(texts),
            block_hash = EX.block_hash(plan.block), indicators = ind)
end

# The log of p3-abstention regenerated for D125_CODE: the rounds and decisions of the tracked log, the panel and preregistration hashes of the tracked files, and the
# ref digests computed round by round from the plan of the log so far (the ref chunks of rounds <= r do not depend on later rounds, D94l).
function d125_synth_log(path::AbstractString)
    pre, pan, tracked = d125_study_paths("p3-abstention")
    tracked = joinpath(D125_ROOT, "workspaces", "study", "stages", "03-run", "output", "p3-abstention-round-log.txt")
    specs = [(round = e.round, decision = e.decision, counts = e.counts) for e in EX.read_round_log(tracked).entries]
    psha = bytes2hex(sha256(read(pan))); hsha = ex_header_hash(read(pre, String))
    texts = String[]; digests = String[]
    for r in 1:length(specs)
        write(path, join(texts))
        plan = EX.plan_run(pre, pan; code_version = D125_CODE, round_log = path)
        push!(digests, EX._ref_digest(plan.chunks, plan.block, r))
        prev = isempty(texts) ? D94_EMPTY_SHA : bytes2hex(sha256(texts[end]))
        push!(texts, d94_entry(; round = specs[r].round, decision = specs[r].decision, panel = psha, prereg = hsha, digest = digests[r], prev = prev, counts = specs[r].counts))
    end
    write(path, join(texts))
    return join(texts)
end
