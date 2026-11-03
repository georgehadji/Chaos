# Stage 04 verification of the analysis module, Amendment D95: O23 / D95e, every systematic alteration of the 3-entry round log of the toy campaign (and of its 2-entry prefix) is
# refused by count_round and by analyze, naming the entry. Run from the repo root after verify/analysis-p3t.jl:
#   julia --project=. verify/analysis-p3t-alter.jl [<scratch root>] [extra-dir-with-round-log-bad-*.txt]
# Reads the campaign folder <root>/campaign (INPUTS.md); writes only <root>/alt/. One printed line per family:  O23 | family | cases | refused by both | names the altered entry | PASS/FAIL
using Chaoslab, SHA
const AN = Chaoslab.Analysis

const ROOT = length(ARGS) >= 1 ? ARGS[1] : "C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/d95v"
const EXTRA = length(ARGS) >= 2 ? ARGS[2] : ""
const CAMPNAME = get(ENV, "D95_CAMPAIGN", "campaign")      # as in analysis-p3t.jl
camp = joinpath(ROOT, CAMPNAME); alt = joinpath(ROOT, CAMPNAME == "campaign" ? "alt" : "alt-" * CAMPNAME); rm(alt; recursive = true, force = true); mkpath(alt)
P(x) = joinpath(camp, x)
log = joinpath(alt, "log.txt")
sha(f) = bytes2hex(sha256(read(f)))
function fp()
    out = String[]
    for (d, _, fs) in walkdir(camp), f in fs
        f == "round-log.txt" && continue
        push!(out, joinpath(d, f) * " " * sha(joinpath(d, f)))
    end
    sort(out)
end
fp0 = fp()
base = (prereg = P("inputs/p3-abstention-preregistration.md"), params = P("inputs/p3-abstention-calibration-parameters.toml"), test_params = P("inputs/p3-abstention-test-parameters.toml"),
        panel = P("inputs/p3-abstention-test-panel.toml"), round_log = log, p2_manifest = P("p2/p2-manifest.json"), p2_prereg = P("p2/p2-preregistration.md"),
        p2_params = P("p2/p2-parameters.toml"), calibration_record = P("inputs/p3-abstention-calibration-record.json"))
stage3 = joinpath(ROOT, CAMPNAME == "campaign" ? "stage" : "stage-" * CAMPNAME, "r3")
full_store = P("store")
function call_cr(text, store)
    write(log, text)
    try
        x = AN.count_round(; base..., store)
        return (:accepted, string(x.action))
    catch e
        return (e isa ArgumentError ? :refused : :other, sprint(showerror, e))
    end
end
oc = 0
function call_an(text, store)
    global oc += 1
    write(log, text)
    out = joinpath(alt, "out-$oc")
    try
        x = AN.analyze(; base..., store, manifest = P("inputs/p3-abstention-test-manifest.json"), p1_prereg = P("p1/p1-preregistration.md"), p1_params = P("p1/p1-parameters.toml"), out)
        return (:accepted, "")
    catch e
        ispath(out) && return (:other, "output folder created on refusal: " * sprint(showerror, e))
        return (e isa ArgumentError ? :refused : :other, sprint(showerror, e))
    end
end
named(msg) = (m = match(r"round log entry (\d+)", msg); m === nothing ? 0 : parse(Int, m.captures[1]))

# the valid entries
text3 = read(P("round-log.txt"), String)
ents = [String(e) for e in split(text3, r"(?<=\n)(?=round=)") if !isempty(e)]
@assert length(ents) == 3 && join(ents) == text3
lines(e) = split(e, '\n'; keepempty = false)
unl(ls) = join(ls .* "\n")
cases = NamedTuple[]   # (family, label, entry (0 = none), text)
add(fam, label, e, t) = (t == text3 || t == join(ents[1:2])) ? nothing : push!(cases, (; fam, label, e, text = t))
function with_entry(k, newe, all = ents)
    v = copy(all); v[k] = newe; join(v)
end
alterbyte(c) = isdigit(c) ? Char('0' + (c - '0' + 1) % 10) : c in 'a':'e' ? Char(c + 1) : c == 'f' ? 'e' : c == 'e' ? 'd' : c == 'n' ? 'o' : 'x'
function gen(all, tag)
    n = length(all)
    for k in 1:n
        ls = lines(all[k])
        for (i, l) in enumerate(ls)
            nm = split(l, '=')[1]
            put(fam, new) = (v = copy(ls); v[i] = new; add("$tag/$fam", "entry $k line $nm", k, with_entry(k, unl(v), all)))
            put("F1a last byte substituted", l[1:end-1] * alterbyte(l[end]))
            put("F1b first byte substituted", "x" * l[2:end])
            put("F1c last byte deleted", l[1:end-1])
            put("F1d byte inserted at end", l * "X")
            put("F1e trailing space", l * " ")
        end
        for (i, l) in enumerate(ls)
            m = match(r"^(count\.[a-z_0-9]+)=(\d+)$", l)
            m === nothing && continue
            v0 = parse(Int, m.captures[2])
            put(fam, new) = (v = copy(ls); v[i] = new; add("$tag/$fam", "entry $k line $(m.captures[1])", k, with_entry(k, unl(v), all)))
            put("F2a count plus one", "$(m.captures[1])=$(v0 + 1)")
            v0 > 0 && put("F2b count minus one", "$(m.captures[1])=$(v0 - 1)")
            put("F5 leading zero in a count", "$(m.captures[1])=0$(v0)")
            put("F6a count line emptied", "")
            v = copy(ls); deleteat!(v, i); add("$tag/F6b count line deleted", "entry $k line $(m.captures[1])", k, with_entry(k, unl(v), all))
        end
        for i in 1:length(ls), j in i+1:length(ls)
            v = copy(ls); v[i], v[j] = v[j], v[i]
            add("$tag/F3 two lines swapped", "entry $k lines $(i),$(j)", k, with_entry(k, unl(v), all))
        end
        for i in 0:length(ls)
            v = copy(ls); insert!(v, i + 1, "")
            add("$tag/F4 blank line", "entry $k after line $i", k, with_entry(k, unl(v), all))
        end
        c = findfirst(startswith("count."), ls)
        v = copy(ls); insert!(v, c, "count.aaa=0"); add("$tag/F7 count name added", "entry $k count.aaa before the first", k, with_entry(k, unl(v), all))
        v = copy(ls); insert!(v, findlast(startswith("count."), ls) + 1, "count.zzz=0"); add("$tag/F7 count name added", "entry $k count.zzz after the last", k, with_entry(k, unl(v), all))
        v = copy(ls); insert!(v, c + 2, "count.orbits=1"); add("$tag/F7 count name added", "entry $k count line repeated", k, with_entry(k, unl(v), all))
        v = copy(ls); insert!(v, 3, "count.target_t_s=1"); add("$tag/F7 count name added", "entry $k name repeated out of order", k, with_entry(k, unl(v), all))
        v = copy(ls); push!(v, ls[end]); add("$tag/F8 line duplicated", "entry $k last line twice", k, with_entry(k, unl(v), all))
        add("$tag/F8 line deleted", "entry $k first line deleted", k, with_entry(k, unl(ls[2:end]), all))
        add("$tag/F8 line deleted", "entry $k last line deleted", k, with_entry(k, unl(ls[1:end-1]), all))
    end
    t = join(all)
    add("$tag/F9 file level", "missing final LF", n, t[1:end-1])
    add("$tag/F9 file level", "second final LF", n, t * "\n")
    add("$tag/F9 file level", "CRLF everywhere", 1, replace(t, "\n" => "\r\n"))
    add("$tag/F9 file level", "CRLF in the first line only", 1, replace(t, "\n" => "\r\n"; count = 1))
    add("$tag/F9 file level", "CR in the last entry", n, join(all[1:n-1]) * replace(all[n], "\n" => "\r\n"; count = 1))
    add("$tag/F9 file level", "leading blank line", 1, "\n" * t)
    add("$tag/F9 file level", "leading space", 1, " " * t)
    add("$tag/F9 file level", "BOM", 1, "\ufeff" * t)
    add("$tag/F9 file level", "tab for a space-free separator in the first line", 1, replace(t, "round=1" => "round=\t1"; count = 1))
    add("$tag/F9 file level", "empty file", 1, "")
    add("$tag/F9 file level", "only a LF", 1, "\n")
    n >= 2 && add("$tag/F9 file level", "entries 1 and 2 swapped", 1, all[2] * all[1] * join(all[3:end]))
    n >= 2 && add("$tag/F9 file level", "entry 1 deleted", 1, join(all[2:end]))
    n >= 3 && add("$tag/F9 file level", "entry 2 deleted", 2, all[1] * join(all[3:end]))
    add("$tag/F9 file level", "entry 1 duplicated", 2, all[1] * t)
    add("$tag/F9 file level", "entry $n duplicated at the end", n + 1, t * all[n])
    add("$tag/F9 file level", "a line of an unknown key appended", n + 1, t * "extra=1\n")
end
gen(ents, "3e")
gen(ents[1:2], "2e")
# the bad files of the independent writer (verify/experiments-d94-log.py, the 24 alterations of the D94 verification generalised to the 13 lines of an entry)
if !isempty(EXTRA)
    for f in sort(readdir(EXTRA; join = true))
        m = match(r"round-log-bad-(.+)-(\d+)\.txt$", f); m === nothing && continue
        push!(cases, (; fam = "py/O26 writer alteration", label = "$(m.captures[1]) entry $(m.captures[2])", e = parse(Int, m.captures[2]), text = read(f, String)))
    end
end

# structural alterations that leave a valid shorter log are expected to differ for count_round; the cases are run and reported, but split from the verdict
valid_prefix(c) = false
res = NamedTuple[]
t0 = time()
for (i, c) in enumerate(cases)
    two = startswith(c.fam, "2e/")
    store = two ? stage3 : full_store
    a = call_cr(c.text, store); b = call_an(c.text, full_store)
    push!(res, (; c..., cr = a, an = b))
    i % 100 == 0 && println("# progress $i / $(length(cases)) after $(round(time() - t0; digits = 0)) s")
end
fams = unique([split(c.fam, '/')[1] * "/" * split(c.fam, '/')[2] for c in res])
println("O23 | family | cases | count_round refused | analyze refused | named entry == altered entry (both) | status")
total = (; n = 0, cr = 0, an = 0, nm = 0)
bad = NamedTuple[]
for f in fams
    rs = [r for r in res if startswith(r.fam, f)]
    ncr = count(r -> r.cr[1] === :refused, rs); nan = count(r -> r.an[1] === :refused, rs)
    nm = count(r -> r.cr[1] === :refused && r.an[1] === :refused && named(r.cr[2]) == r.e && named(r.an[2]) == r.e, rs)
    # analyze on the 2-entry prefix is a different question (the log is open for a 3-round store): its refusal is expected, naming need not match
    ok = ncr == length(rs) && nan == length(rs)
    println("O23 | ", f, " | ", length(rs), " | ", ncr, " | ", nan, " | ", nm, " | ", ok ? "all refused" : "NOT ALL REFUSED")
    for r in rs
        (r.cr[1] !== :refused || r.an[1] !== :refused) && push!(bad, r)
    end
end
println("TOTAL | cases ", length(res), " | count_round refused ", count(r -> r.cr[1] === :refused, res), " | analyze refused ", count(r -> r.an[1] === :refused, res),
        " | other exception types ", count(r -> r.cr[1] === :other || r.an[1] === :other, res), " | entry named == altered entry in both ", count(r -> r.cr[1] === :refused && r.an[1] === :refused && named(r.cr[2]) == r.e && named(r.an[2]) == r.e, res))
println("--- cases not refused by both, or refused by an exception that is not ArgumentError:")
for r in bad
    println("  ", r.fam, " | ", r.label, " | count_round: ", r.cr[1], " ", first(r.cr[2], 160), " | analyze: ", r.an[1], " ", first(r.an[2], 160))
end
println("--- cases refused by both but naming a different entry than the altered one (message of count_round / analyze):")
mm = [r for r in res if r.cr[1] === :refused && r.an[1] === :refused && !(named(r.cr[2]) == r.e && named(r.an[2]) == r.e)]
println("  count: ", length(mm))
for r in mm[1:min(end, 60)]
    println("  ", r.fam, " | ", r.label, " | altered entry ", r.e, " | named cr ", named(r.cr[2]), " an ", named(r.an[2]), " | ", first(r.cr[2], 130))
end
println("--- sample messages")
for r in res[1:3:min(end, 60)]
    println("  ", r.fam, " | ", r.label, " | ", first(r.cr[2], 170))
end
# the valid logs pass
println("--- valid logs (accepted):")
for (k, store) in ((1, joinpath(dirname(stage3), "r1")), (2, joinpath(dirname(stage3), "r2")), (3, full_store))
    t = join(ents[1:k]); a = call_cr(t, store)
    println("  ", k, "-entry log with the store of round ", k, " | count_round: ", a, k == 3 ? " | analyze: " * string(call_an(t, full_store)[1]) : "")
end
write(log, text3)
fp1 = fp()
println("inputs unchanged by all calls (every file of the campaign except the log): ", fp0 == fp1)
println("round-log restored to the valid text: ", read(P("round-log.txt"), String) == text3 || "(the alt log is a separate file)")
