# O24 (V), runner side: the round log written by count-round for the toy campaign is accepted by the log reader and checker of the runner (Experiments.read_round_log, check_round_log).
# The toy store was written by the chunk writer, not by the runner, so the full run_campaign / `run --round-log` call (which re-plans the chunks and compares them with the store) is not
# possible on it; this script applies the runner's own parser and checker, with the digest taken from the chunk hashes of the store, to the log and to the 39 alterations of the
# independent writer (verify/experiments-d94-log.py). Run: julia --project=. verify/analysis-p3t-runner.jl [<scratch root>]
using Chaoslab, SHA
const XP = Chaoslab.Experiments
const ROOT = length(ARGS) >= 1 ? ARGS[1] : "C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/d95v"
camp = joinpath(ROOT, "campaign")
files = sort([joinpath(camp, "store", f) for f in readdir(joinpath(camp, "store")) if endswith(f, ".h5")])
chunks = [(a = XP.read_chunk_attrs(f)) for f in files]
key(a) = String(a["chunk_key"]); hsh(a) = String(a["chunk_hash"])
rnd(k) = (m = match(r":r(\d+)$", k); m === nothing ? nothing : parse(Int, m.captures[1]))
ref(r) = sort([hsh(a) for a in chunks if (q = rnd(key(a))) !== nothing && q <= r && occursin(r"method=L[15] ", key(a))])
digest(r) = bytes2hex(sha256(join(h * "\n" for h in ref(r))))
panel = bytes2hex(sha256(read(joinpath(camp, "inputs", "p3-abstention-test-panel.toml"))))
prereg = match(r"^prereg_sha256=(.*)$"m, read(joinpath(camp, "round-log.txt"), String)).captures[1]
function runner_check(path)
    try
        lg = XP.read_round_log(path)
        st = XP.check_round_log(lg.entries; path, r_max = 3, panel, prereg, digest)
        return (:accepted, string(st))
    catch e
        return (e isa ArgumentError ? :refused : :other, sprint(showerror, e))
    end
end
x = runner_check(joinpath(camp, "round-log.txt"))
println("O24 | runner reader and checker on the log of count-round (3 entries, r_max 3, panel and prereg hashes, ref digest from the store chunk hashes) | ", x[1], " ", x[2], " | accepted | ", x[1] === :accepted ? "PASS" : "FAIL")
dir = joinpath(ROOT, "o24")
nacc = 0; nref = 0; accepted_names = String[]; other = 0
for f in sort(readdir(dir; join = true))
    occursin("round-log-bad-", f) || continue
    r = runner_check(f)
    r[1] === :accepted ? (global nacc += 1; push!(accepted_names, basename(f))) : r[1] === :refused ? (global nref += 1) : (global other += 1)
end
println("O24 | runner on the 39 one-byte alterations of the independent writer | refused ", nref, ", accepted ", nacc, ", other exceptions ", other, " | accepted are the 7 count lines of the last entry (the runner does not recount; a count line of entries 1 and 2 breaks the chain), as the 2 of D94 O26 | ", nacc == 7 && other == 0 ? "PASS" : "CHECK")
println("      accepted: ", join(accepted_names, ", "))
