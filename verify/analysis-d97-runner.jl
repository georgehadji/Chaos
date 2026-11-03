# O28 (V), runner side, on the store of the runner (Amendment D97; the D95 script verify/analysis-p3t-runner.jl is for the toy store and is not changed): the round log written by count_round / count-round
# is accepted by the log reader and checker of the runner (Experiments.read_round_log, check_round_log), with the digest taken from the chunk files of the arm ref of the store, and the 39 one-byte
# alterations of the independent writer (verify/experiments-d94-log.py, output folder of this script's second argument) are refused except the count lines of the last entry (the runner does not recount).
# Run: julia --project=. verify/analysis-d97-runner.jl <fixture dir> <out dir> <tag> <o24 folder with round-log-bad-*.txt>
using Chaoslab, SHA
const XP = Chaoslab.Experiments
FIX, OUTD, TAG, O24 = ARGS[1:4]
store = joinpath(OUTD, "store-$TAG")
files = sort([joinpath(store, f) for f in readdir(store) if endswith(f, ".h5")])
chunks = [XP.read_chunk_attrs(f) for f in files]
key(a) = String(a["chunk_key"]); hsh(a) = String(a["chunk_hash"]); arms(a) = String.(collect(a["arms"]))
rnd(k) = (m = match(r":r(\d+)$", k); m === nothing ? nothing : parse(Int, m.captures[1]))
ref(r) = sort([hsh(a) for a in chunks if arms(a) == ["ref"] && (q = rnd(key(a))) !== nothing && q <= r])
digest(r) = bytes2hex(sha256(join(h * "\n" for h in ref(r))))
panel = bytes2hex(sha256(read(joinpath(FIX, "p3-abstention-test-panel.toml"))))
logf = joinpath(OUTD, "round-log-$TAG.txt")
prereg = match(r"^prereg_sha256=(.*)$"m, read(logf, String)).captures[1]
function runner_check(path)
    try
        lg = XP.read_round_log(path)
        st = XP.check_round_log(lg.entries; path, r_max = 3, panel, prereg, digest)
        return (:accepted, string(st))
    catch e
        return (e isa ArgumentError ? :refused : :other, sprint(showerror, e))
    end
end
x = runner_check(logf)
println("O28 | runner reader and checker on the log of count_round (3 entries, r_max 3, panel and prereg hashes, ref digest from the chunk files of the arm ref) | ", x[1], " ", x[2], " | accepted | ", x[1] === :accepted ? "PASS" : "FAIL")
nacc = 0; nref = 0; other = 0; accepted_names = String[]
for f in sort(readdir(O24; join = true))
    occursin("round-log-bad-", f) || continue
    r = runner_check(f)
    r[1] === :accepted ? (global nacc += 1; push!(accepted_names, basename(f))) : r[1] === :refused ? (global nref += 1) : (global other += 1)
end
println("O28 | runner on the 39 one-byte alterations of the independent writer | refused ", nref, ", accepted ", nacc, ", other exceptions ", other, " | accepted are the 7 count lines of the last entry, as D94 O26 and D95 O24 | ", nacc == 7 && other == 0 ? "PASS" : "CHECK")
println("      accepted: ", join(accepted_names, ", "))
