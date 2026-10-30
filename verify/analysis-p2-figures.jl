# Stage 04 verification of the analysis module, Amendment D91 (figure drawers F1 to F6 of study p2-indicators). Run from the repo root:
#   julia --project=verify -t 23 verify/analysis-p2-figures.jl
# Written from the spec only (workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D91": E38 to E45, I22 to I29, O11 to O14,
# criteria D91a to D91o) and from the accepted tables file of study stage 04. Nothing of src/, ext/ or test/ was read; the only look at the
# implementation is run time introspection of the documented hook (the method signature of the function stored in Chaoslab.Analysis.FIGURE_DRAWER).
# PDF measurements are made by verify/analysis-p2-figures.py (PyMuPDF), called as `python -I`.
# Sections (ANV_ONLY="HASH0,BLOCKS" runs only those): HASH0, BLOCKS, CLIFIG, CLINOFIG, CLIFIG2, API, HOOK, FIXT, P1, P1PRE, MISC, CONTACT, HASH1.
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL/INFO
# Nothing in data/, src/, ext/, the studies or the specs is written; scratch outputs go to the scratch folder SCR; the log and the contact sheets go to
# the stage 04 output folder of the module-build workspace.
using Printf, SHA, TOML, Dates
using Chaoslab, CairoMakie

const AN = Chaoslab.Analysis
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const STUDY = joinpath(ROOT, "workspaces", "study", "stages")
const PREREG = joinpath(STUDY, "01-design", "output", "p2-indicators-preregistration.md")
const MAN2 = joinpath(STUDY, "03-run", "output", "p2-indicators-manifest.json")
const PANEL = joinpath(STUDY, "02-panel", "output", "p1-le-benchmark-panel.toml")
const PARAMS = joinpath(STUDY, "04-analyze", "references", "p2-indicators-analysis-parameters.toml")
const PREREG1 = joinpath(STUDY, "01-design", "output", "p1-le-benchmark-preregistration.md")
const MAN1 = joinpath(STUDY, "03-run", "output", "p1-le-benchmark-manifest.json")
const PARAMS1 = joinpath(STUDY, "04-analyze", "references", "p1-le-benchmark-analysis-parameters.toml")
const ACC1 = joinpath(STUDY, "04-analyze", "output", "p1-le-benchmark-tables.md")
const ACC2 = joinpath(STUDY, "04-analyze", "output", "p2-indicators-tables.md")
const TRACKED1 = joinpath(STUDY, "04-analyze", "output")
const PAPERFIG = joinpath(ROOT, "workspaces", "paper", "stages", "04-finalize", "output", "p1-le-benchmark", "figures")
const DEFAULTS = joinpath(ROOT, "_config", "numerics-defaults.toml")
const STORE2 = joinpath(ROOT, "data", "p2-indicators")
const STORE1 = joinpath(ROOT, "data", "p1-le-benchmark")
const BIN = joinpath(ROOT, "bin", "chaoslab.jl")
const SCR = raw"C:\Users\tesse\AppData\Local\Temp\claude\C--Users-tesse-Documents-Vibe-Coding-Chaos\5df15ee1-47df-4875-85bd-1bd4d7c665a6\scratchpad\d91v"
const PY = "python"
const HELPER = joinpath(@__DIR__, "analysis-p2-figures.py")
const JL = first(Base.julia_cmd().exec)
const NT = Threads.nthreads()
const ONLY = filter(!isempty, split(get(ENV, "ANV_ONLY", ""), ","))
const LOG = open(get(ENV, "ANV_LOG", joinpath(OUT, "analysis-d91-verify-log.txt")), get(ENV, "ANV_APPEND", "") == "1" ? "a" : "w")
const IND = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA", "RLI"]
const FIGS = ["F$k" for k in 1:6]
const EXPECTED_ROWS = Dict("F1" => 112, "F2" => 528, "F3" => 320, "F4" => 2700, "F5" => 2304, "F6" => 489)
const D91C_DIGESTS = Dict("F1" => "36beb051854affd330ee00e0c5b6edfdf651f6e582fe7f8bb037db3ff6ad5ecb", "F2" => "eafb235298146784c7dfc14100f6bfe5bf88ddae022a1bb403a2cda261f807ad",
    "F3" => "3157b5b5592456c7847b93ad5bc00552486113090d6525c076837bae0f8298fa", "F4" => "7b96868bde4c93a1dcb62289eafc8b867fe82606f710a3e55c76531ced707c21",
    "F5" => "63babbf4bd5a960af9314a79d8886752dceeac9f8399bf3165c3c401cb71f05e", "F6" => "f2d0dabc09bba4137f90c9aa4ff7492946445625f196f1a341f2c5a7c7e23967")
const PDFNAMES2 = ["p2-indicators-fig-F$k.pdf" for k in 1:6]
const PDFNAMES1 = vcat(["p1-le-benchmark-fig-F$k.pdf" for k in 1:7], ["p1-le-benchmark-fig-$v.pdf" for v in ("F1a", "F1b", "F2a", "F2b")])
mkpath(SCR)
mkpath(OUT)

fmt(x::AbstractFloat) = @sprintf("%.6g", x)
fmt(x) = string(x)
const COUNTS = Dict("PASS" => 0, "FAIL" => 0, "INFO" => 0)
const VERD = Dict{String,Vector{String}}()
function say(args...)
    s = string(args...)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(id, quantity, value, bound, verdict)
    v = startswith(verdict, "FAIL") ? "FAIL" : startswith(verdict, "PASS") ? "PASS" : "INFO"
    COUNTS[v] += 1
    push!(get!(VERD, string(id), String[]), v)
    say(string(id, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
want(sec) = isempty(ONLY) || sec in ONLY
function section(f, name)
    want(name) || return
    say("== ", name, " ==")
    t0 = time()
    try
        f()
    catch e
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 1800)], "-", "FAIL")
    end
    say(@sprintf("== %s done in %.1f s", name, time() - t0))
end

# ---------------------------------------------------------------------------------------------------------------- helpers: python bridge, tables parse, digests
function py(args...)
    out = IOBuffer()
    run(pipeline(ignorestatus(Cmd(String[PY, "-I", HELPER, string.(args)...])); stdout = out, stderr = devnull))
    String(take!(out))
end
function pyemit(text)
    n = 0
    for l in split(text, '\n')
        startswith(l, "@@ ") || continue
        parts = split(l[4:end], " | ")
        if length(parts) < 5
            emit("PYERR", "unparsable line of the PDF reader", l, "-", "FAIL"); continue
        end
        emit(parts[1], join(parts[2:end-3], " | "), parts[end-2], parts[end-1], parts[end])
        n += 1
    end
    n
end

function parse_blocks(path)
    L = split(replace(read(path, String), "\r\n" => "\n"), '\n')
    blocks = Dict{String,Any}()
    cells(l) = String[strip(c) for c in split(l, '|')[2:end-1]]
    i = 1
    while i <= length(L)
        m = match(r"^## (F[1-6]) data", L[i])
        if m !== nothing
            j = i + 1
            while isempty(L[j]); j += 1; end
            k = j
            while startswith(L[k], "|"); k += 1; end
            d = k
            while !startswith(L[d], "digest:"); d += 1; end
            blocks[m[1]] = (cols = cells(L[j]), raw = [cells(l) for l in L[j+2:k-1]], digest = String(match(r"digest: ([0-9a-f]{64})", L[d])[1]))
            i = d
        end
        i += 1
    end
    blocks
end
typedv(s) = occursin(r"^-?\d+$", s) ? parse(Int, s) : s == "nothing" ? nothing : (x = tryparse(Float64, s); x === nothing ? String(s) : x)
fmtcell(v) = v === nothing ? "nothing" : v isa Integer ? string(v) : v isa AbstractFloat ? string(round(v; sigdigits = 12)) : string(v)
# E37 rule: SHA-256 of the block written as a markdown table, numbers with 12 significant digits, LF line ends, no final line end
function digest_of(cols, rows)
    lines = String["| " * join(cols, " | ") * " |", "|" * join(fill("---", length(cols)), "|") * "|"]
    for r in rows
        push!(lines, "| " * join((fmtcell(v) for v in r), " | ") * " |")
    end
    bytes2hex(sha256(join(lines, "\n")))
end
const BLK = parse_blocks(ACC2)
typedrows(k) = Vector{Any}[Any[typedv(c) for c in r] for r in BLK[k].raw]
const TOML2 = TOML.parsefile(PARAMS)
const PRM = (budgets = Vector{Int}(TOML2["budgets"]["value"]), b_primary = Int(TOML2["b_primary"]["value"]), N_c = Int(TOML2["N_c"]["value"]))
basefigs() = Dict{String,Any}(k => (columns = copy(BLK[k].cols), rows = typedrows(k), extra = String[]) for k in FIGS)
resultlike(figs) = (study = "p2-indicators", figures = figs, params = PRM)

function tdiff(newf, oldf)
    A = split(replace(read(newf, String), "\r\n" => "\n"), '\n')
    B = split(replace(read(oldf, String), "\r\n" => "\n"), '\n')
    diffs = [i for i in 1:max(length(A), length(B)) if get(A, i, nothing) != get(B, i, nothing)]
    hdr(l) = startswith(l, "- module version:") || startswith(l, "- git commit:")
    ok = length(A) == length(B) && all(i -> hdr(get(B, i, "")) && hdr(get(A, i, "")), diffs)
    (na = length(A), nb = length(B), diffs = diffs, ok = ok, old = B, new = A)
end
digestlines(f) = [m.captures[1] for l in split(replace(read(f, String), "\r\n" => "\n"), '\n') for m in (match(r"^digest: ([0-9a-f]{64})", l),) if m !== nothing]

function runcli(name, proj, out, args; bin = BIN)
    rm(out; recursive = true, force = true)
    logf = joinpath(SCR, name * "-stdout.txt"); resf = joinpath(SCR, name * "-result.txt")
    rm(resf; force = true)
    cmd = String[JL, "--project=$proj", "-t", string(NT), bin, "analyze", args..., "--out", out]
    run(pipeline(ignorestatus(Cmd(String[PY, "-I", HELPER, "run", logf, resf, ROOT, "--", cmd...])); stderr = devnull))
    kv = Dict(String(split(l, '=')[1]) => String(split(l, '=')[2]) for l in readlines(resf) if occursin('=', l))
    txt = read(logf, String)
    occursin("Bad Image", txt) || occursin("0xc0e90002", txt) ? say("NOTE: Windows application control block seen in the output of run ", name, ": delete the blocked cache file, wait one minute, rerun") : nothing
    (exit = parse(Int, kv["exit"]), wall = parse(Float64, kv["wall"]), peak = parse(Int, kv["peak_wset"]), text = txt, out = out)
end
args2() = String[PREREG, MAN2, PANEL, PARAMS, "--store", STORE2]
args1() = String[PREREG1, MAN1, PANEL, PARAMS1, "--store", STORE1]
gb(x) = @sprintf("%.2f GB", x / 2^30)

# the state of the runs, shared by the sections (a section run alone reads what an earlier run left in SCR)
const RUN = Dict{String,Any}()

# ---------------------------------------------------------------------------------------------------------------- fingerprints (D91o, I26)
files_under(d) = sort([joinpath(r, f) for (r, _, fs) in walkdir(d) for f in fs])
function snapshot()
    data = vcat(files_under(STORE2), files_under(STORE1))
    other = vcat([MAN2, PANEL, PREREG, PARAMS, PREREG1, MAN1, PARAMS1, ACC1, ACC2], isfile(DEFAULTS) ? [DEFAULTS] : String[],
        [joinpath(TRACKED1, n) for n in PDFNAMES1[1:7]], [joinpath(PAPERFIG, n) for n in PDFNAMES1])
    h(files) = Dict(f => open(io -> bytes2hex(sha256(io)), f) for f in files if isfile(f))
    tree = String[]
    for (r, _, fs) in walkdir(joinpath(ROOT, "data")), f in fs
        p = joinpath(r, f)
        push!(tree, string(relpath(p, ROOT), " ", filesize(p), " ", floor(Int, mtime(p))))
    end
    (data = h(data), other = h(other), tree = bytes2hex(sha256(join(sort(tree), "\n"))), ntree = length(tree))
end
const SNAP0 = Ref{Any}(nothing)

section("HASH0") do
    t = @timed snapshot()
    SNAP0[] = t.value
    s = t.value
    emit("D91o", "fingerprints before the runs: files under data/p2-indicators and data/p1-le-benchmark (SHA-256 each), the other inputs, entries in the listing of the data tree and its SHA-256 (path, size, mtime)",
         "$(length(s.data)) data files, $(length(s.other)) other files, $(s.ntree) tree entries, tree $(s.tree[1:16]); $(round(t.time; digits = 1)) s", "recorded", "INFO")
end

# ---------------------------------------------------------------------------------------------------------------- BLOCKS (D91b, O11: digests recomputed by my own code)
section("BLOCKS") do
    for k in FIGS
        b = BLK[k]
        rows = typedrows(k)
        mine = digest_of(b.cols, rows)
        dif = sum(count(j -> fmtcell(typedv(b.raw[i][j])) != b.raw[i][j], eachindex(b.raw[i])) for i in eachindex(b.raw))
        emit("D91b", "$k: digest recomputed from the numbers parsed out of the accepted tables file (own code, E37 rule) against the digest line of the file and the value of D91c; rows; cells whose 12-digit rewriting differs from the file text",
             "$(mine[1:16]), file $(b.digest[1:16]), D91c $(D91C_DIGESTS[k][1:16]); $(length(rows)) rows; $dif cells differ", "equal; $(EXPECTED_ROWS[k]) rows; 0",
             pf(mine == b.digest == D91C_DIGESTS[k] && length(rows) == EXPECTED_ROWS[k] && dif == 0))
    end
end

# ---------------------------------------------------------------------------------------------------------------- CLIFIG (D91a, D91c, D91d, D91f to D91j: the real run, verify project)
function check_run_files(out, names, label)
    have = sort(readdir(out))
    expect = sort(vcat([split(names[1], "-fig-")[1] * "-tables.md"], names))
    emit("D91o", "$label: files in the output folder", join(have, " "), "tables file and the figure files only", pf(have == expect))
    have == expect
end
function check_tables(newf, label)
    d = tdiff(newf, ACC2)
    emit("D91c", "$label: tables file against the accepted stage 04 file, line by line (lines in new and old; lines that differ)", "$(d.na) and $(d.nb); $(length(d.diffs)): " * join(["line $i: \"$(d.old[i])\" -> \"$(get(d.new, i, ""))\"" for i in d.diffs], "; "), "only the lines module version and git commit", pf(d.ok))
    dg = digestlines(newf)
    emit("D91c", "$label: the six digest lines equal the six values of D91c in order", "$(length(dg)) lines; equal $(dg == [D91C_DIGESTS[k] for k in FIGS])", "equal", pf(dg == [D91C_DIGESTS[k] for k in FIGS]))
end

section("CLIFIG") do
    out = joinpath(SCR, "run1")
    r = runcli("run1", joinpath(ROOT, "verify"), out, args2())
    RUN["run1"] = r
    emit("D91a", "bin/chaoslab.jl analyze on the real P2 inputs under the verify project (23 threads): exit code", r.exit, "0", pf(r.exit == 0))
    emit("D91d", "real inputs: the output has neither the line 'figures skipped' nor 'figures failed'", "skipped $(occursin("figures skipped", r.text)), failed $(occursin("figures failed", r.text)); output: " * replace(strip(r.text), '\n' => " / ")[1:min(end, 400)], "false, false", pf(!occursin("figures skipped", r.text) && !occursin("figures failed", r.text)))
    check_run_files(out, PDFNAMES2, "run 1")
    check_tables(joinpath(out, "p2-indicators-tables.md"), "run 1 (with figures)")
    emit("D91j", "real run: wall time of the whole command (package loading included), peak working set of the process", "$(round(r.wall; digits = 1)) s, $(gb(r.peak))", "whole run <= 8 GB", pf(r.peak <= 8 * 2^30))
    tbl = joinpath(out, "p2-indicators-tables.md")
    pyemit(py("figs", out, tbl, PARAMS))
end

section("CLINOFIG") do
    out = joinpath(SCR, "nofig")
    r = runcli("nofig", ROOT, out, args2())
    RUN["nofig"] = r
    files = isdir(out) ? sort(readdir(out)) : String[]
    emit("D91d", "default project (no plotting package): exit code; the output says 'figures skipped'; files written", "$(r.exit); $(occursin("figures skipped", r.text)); $(join(files, " "))", "2; true; the tables file only", pf(r.exit == 2 && occursin("figures skipped", r.text) && files == ["p2-indicators-tables.md"]))
    isfile(joinpath(out, "p2-indicators-tables.md")) && check_tables(joinpath(out, "p2-indicators-tables.md"), "run without figures")
    r1 = get(RUN, "run1", nothing)
    if isfile(joinpath(SCR, "run1", "p2-indicators-tables.md")) && isfile(joinpath(out, "p2-indicators-tables.md"))
        d = tdiff(joinpath(SCR, "run1", "p2-indicators-tables.md"), joinpath(out, "p2-indicators-tables.md"))
        emit("I23", "tables file of the run with figures against the one without figures: lines that differ", "$(length(d.diffs)) lines: " * join(["line $i" for i in d.diffs], ", "), "only module version and git commit", pf(d.ok))
    end
    r1 !== nothing && emit("D91j", "CLI wall time with figures (verify project) and without (default project, no plotting package: the difference includes the loading of the plotting package); peak working set without figures",
         "$(round(r1.wall; digits = 1)) s and $(round(r.wall; digits = 1)) s, difference $(round(r1.wall - r.wall; digits = 1)) s; $(gb(r.peak))", "recorded", "INFO")
end

section("CLIFIG2") do
    out = joinpath(SCR, "run2")
    r = runcli("run2", joinpath(ROOT, "verify"), out, args2())
    RUN["run2"] = r
    emit("D91k", "second run into a fresh folder: exit code, wall time, peak working set", "$(r.exit), $(round(r.wall; digits = 1)) s, $(gb(r.peak))", "0", pf(r.exit == 0))
    check_run_files(out, PDFNAMES2, "run 2")
    t1, t2 = joinpath(SCR, "run1", "p2-indicators-tables.md"), joinpath(out, "p2-indicators-tables.md")
    same = isfile(t1) && isfile(t2) && read(t1) == read(t2)
    emit("D91k", "tables file of run 1 and run 2: byte-identical (including module version and git commit); six digests equal", "$same; $(digestlines(t1) == digestlines(t2))", "true; true", pf(same && digestlines(t1) == digestlines(t2)))
    res = py("cmp", joinpath(SCR, "run1"), out, PDFNAMES2...)
    for l in split(strip(res), '\n')
        p = split(l, '\t')
        length(p) >= 3 || continue
        nd = parse(Int, p[2])
        emit("D91k", "$(p[1]): render at 72 dpi of run 1 against run 2, differing pixels; $(p[3])", nd, "0 (a difference only in embedded metadata is listed, not a failure)", pf(nd == 0))
    end
end

# ---------------------------------------------------------------------------------------------------------------- API: the real command in-process with a recording hook (D91b, I22, D91d, D91j)
function cli_inproc(out, args)
    rm(out; recursive = true, force = true)
    logf = joinpath(SCR, basename(out) * "-inproc.txt")
    local code
    t = @timed open(logf, "w") do io
        redirect_stdout(io) do
            redirect_stderr(io) do
                code = AN.cli(vcat(["analyze"], args, ["--out", out, "--threads", string(NT)]))
            end
        end
    end
    (code = code, time = t.time, text = read(logf, String), out = out)
end

section("API") do
    orig = AN.FIGURE_DRAWER[]
    emit("D91d", "documented hook Chaoslab.Analysis.FIGURE_DRAWER: set when the plotting package is loaded; the function and its methods (introspection only)", string(orig, " ", [string(m.sig) for m in methods(orig)]), "set", pf(orig !== nothing))
    cap = Ref{Any}(nothing); before = Ref{Any}(nothing); dtime = Ref(NaN)
    AN.FIGURE_DRAWER[] = (res, o) -> begin
        cap[] = res; before[] = deepcopy(res.figures)
        t = @timed orig(res, o)
        dtime[] = t.time
        t.value
    end
    local r
    try
        r = cli_inproc(joinpath(SCR, "api1"), args2())
    finally
        AN.FIGURE_DRAWER[] = orig
    end
    res = cap[]
    emit("D91a", "in-process command (the same entry point as the CLI) with a recording hook: exit code; the hook was called", "$(r.code); $(res !== nothing)", "0; true", pf(r.code == 0 && res !== nothing))
    emit("D91d", "in-process run output", replace(strip(r.text), '\n' => " / ")[1:min(end, 300)], "no 'figures skipped' or 'figures failed'", pf(!occursin("figures skipped", r.text) && !occursin("figures failed", r.text)))
    filetab = parse_blocks(joinpath(r.out, "p2-indicators-tables.md"))
    for k in FIGS
        fg = res.figures[k]
        handed = res.digests[k]
        recv = digest_of(fg.columns, fg.rows)
        emit("D91b", "$k: digest handed to the drawer (res.digests), digest recomputed from the rows the drawer received (own code), digest line of the tables file, digest recomputed from the numbers parsed out of the accepted file; rows received",
             "$(handed[1:12]), $(recv[1:12]), $(filetab[k].digest[1:12]), $(digest_of(BLK[k].cols, typedrows(k))[1:12]); $(length(fg.rows))", "all equal; $(EXPECTED_ROWS[k]) rows",
             pf(handed == recv == filetab[k].digest == digest_of(BLK[k].cols, typedrows(k)) && length(fg.rows) == EXPECTED_ROWS[k]))
        nd = sum(count(j -> fmtcell(fg.rows[i][j]) != BLK[k].raw[i][j], eachindex(fg.rows[i])) for i in eachindex(fg.rows))
        unchanged = isequal(before[][k].rows, fg.rows) && before[][k].columns == fg.columns
        emit("I22", "$k: the figure is the block: cells the drawer received that differ from the text of the accepted block (12 digits); rows and columns unchanged by the drawing", "$nd cells differ; unchanged $unchanged; columns equal $(fg.columns == BLK[k].cols)", "0; true; true", pf(nd == 0 && unchanged && fg.columns == BLK[k].cols))
    end
    emit("D91j", "figure step in the real run (hook call, plotting package already loaded, first call so compilation included), 23 threads; the whole in-process command", "$(round(dtime[]; digits = 1)) s; $(round(r.time; digits = 1)) s", "<= 120 s", pf(dtime[] <= 120))
    out2 = joinpath(SCR, "api-warm"); rm(out2; recursive = true, force = true); mkpath(out2)
    tw = @timed orig(res, out2)
    emit("D91j", "figure step, second call (warm)", "$(round(tw.time; digits = 1)) s", "recorded", "INFO")
    # failing drawer through the whole command: one non-positive D_str in F6 (the hook is wrapped so that the block it receives has the value)
    AN.FIGURE_DRAWER[] = (res, o) -> begin
        figs = deepcopy(res.figures)
        figs["F6"].rows[100][6] = -1.0   # column p50 of the 100th row of the block of F6
        orig(merge(res, (figures = figs,)), o)
    end
    local rf
    try
        rf = cli_inproc(joinpath(SCR, "apifail"), args2())
    finally
        AN.FIGURE_DRAWER[] = orig
    end
    files = isdir(rf.out) ? sort(readdir(rf.out)) : String[]
    line = [l for l in split(rf.text, '\n') if startswith(l, "figures failed")]
    emit("D91d", "command with a non-positive D_str (-1.0) in row 100 of F6 (value p50): exit code; the line 'figures failed'; the F6 file; the files written", "$(rf.code); $(isempty(line) ? "no line" : strip(line[1])); F6 file $(("p2-indicators-fig-F6.pdf" in files)); $(join(files, " "))",
         "2; names F6 and the row; no F6 file", pf(rf.code == 2 && !isempty(line) && occursin("F6", line[1]) && occursin(r"\d", line[1]) && !("p2-indicators-fig-F6.pdf" in files)))
    if isfile(joinpath(rf.out, "p2-indicators-tables.md"))
        d = tdiff(joinpath(rf.out, "p2-indicators-tables.md"), ACC2)
        emit("D91d", "failed-figure run: the tables file is intact (against the accepted file)", "$(length(d.diffs)) differing lines", "only module version and git commit", pf(d.ok))
    else
        emit("D91d", "failed-figure run: the tables file is intact", "missing", "present", "FAIL")
    end
    emit("I26", "failed-figure run: files written (the figure files drawn before the failure stay or not; the spec says tables stay and the line names the figure)", join(files, " "), "reported", "INFO")
end

# ---------------------------------------------------------------------------------------------------------------- HOOK and fixture cases (D91d direct hook call, D91m)
function drawcase(name, res)
    out = joinpath(SCR, "fix-" * name)
    rm(out; recursive = true, force = true); mkpath(out)
    hook = AN.FIGURE_DRAWER[]
    err = nothing
    t = @timed try
        hook(res, out)
    catch e
        err = sprint(showerror, e)
    end
    files = sort(readdir(out))
    sig = [f for f in files if endswith(f, ".pdf") && read(joinpath(out, f), 5) == b"%PDF-" && filesize(joinpath(out, f)) > 1000]
    (err = err, files = files, sig = sig, out = out, time = t.time)
end

section("HOOK") do
    figs = basefigs()
    c = drawcase("own-parse", resultlike(figs))
    emit("D91d", "hook called with a result built from my own parse of the accepted blocks (study, figures, params = budgets, b_primary, N_c): error; the six files with signature", "$(c.err === nothing ? "none" : c.err); $(length(c.sig)) files; $(round(c.time; digits = 1)) s", "no error; 6 files", pf(c.err === nothing && length(c.sig) == 6))
    figs2 = basefigs()
    figs2["F6"].rows[100][6] = -1.0
    c2 = drawcase("f6-nonpositive", resultlike(figs2))
    msg = c2.err === nothing ? "no error" : c2.err
    emit("D91d", "hook called with a non-positive D_str (-1.0, column p50) in row 100 of F6: the exception message; the F6 file; the files present", "$(msg[1:min(end, 300)]); F6 file $("p2-indicators-fig-F6.pdf" in c2.files); files $(join(c2.files, " "))",
         "an exception that names F6 and a row; no F6 file", pf(c2.err !== nothing && occursin("F6", c2.err) && occursin(r"\d", c2.err) && !("p2-indicators-fig-F6.pdf" in c2.files)))
    figs3 = basefigs()
    figs3["F6"].rows[200][4] = 0.0
    c3 = drawcase("f6-zero", resultlike(figs3))
    emit("D91d", "hook with a zero D_str (0.0, column p5) in row 200 of F6: exception names F6 and a row; no F6 file", "$((c3.err === nothing ? "no error" : c3.err)[1:min(end, 250)]); F6 file $("p2-indicators-fig-F6.pdf" in c3.files)", "exception; no F6 file", pf(c3.err !== nothing && occursin("F6", c3.err) && !("p2-indicators-fig-F6.pdf" in c3.files)))
end

colidx(figs, k, name) = findfirst(==(name), figs[k].columns)
function setrows!(figs, k, rows)
    figs[k] = (columns = figs[k].columns, rows = rows, extra = figs[k].extra)
end
function fixture_cases()
    cases = Pair{String,Function}[]
    push!(cases, "no-rli-in-S2" => f -> begin
        for (k, cs) in (("F1", ("X", "Y")), ("F2", ("indicator",)), ("F3", ("indicator",)))
            si = colidx(f, k, "set"); ci = [colidx(f, k, c) for c in cs]
            setrows!(f, k, filter(r -> !(r[si] == "S2" && any(r[i] == "RLI" for i in ci)), f[k].rows))
        end
    end)
    push!(cases, "no-rli-anywhere" => f -> begin
        for (k, cs) in (("F1", ("X", "Y")), ("F2", ("indicator",)), ("F3", ("indicator",)), ("F4", ("indicator",)))
            ci = [colidx(f, k, c) for c in cs]
            setrows!(f, k, filter(r -> !any(r[i] == "RLI" for i in ci), f[k].rows))
        end
    end)
    push!(cases, "k3-all-auc-nothing" => f -> begin
        ki, ai = colidx(f, "F4", "K"), colidx(f, "F4", "AUC")
        for r in f["F4"].rows
            r[ki] == 3.0 && (r[ai] = nothing)
        end
    end)
    push!(cases, "min-size-set-30-orbits" => f -> begin
        si = colidx(f, "F2", "set"); ii = colidx(f, "F2", "indicator"); ki = colidx(f, "F2", "kind"); vi = colidx(f, "F2", "value"); ci = colidx(f, "F2", "count")
        ts = Dict((r[ii], r[ki], r[vi]) => r[ci] for r in f["F2"].rows if r[si] == "T_s")
        for r in f["F2"].rows
            r[si] == "W2" && (r[ci] = ts[(r[ii], r[ki], r[vi])])
        end
        si, ii, ki2, di, fi = colidx(f, "F3", "set"), colidx(f, "F3", "indicator"), colidx(f, "F3", "k"), colidx(f, "F3", "decided"), colidx(f, "F3", "fraction")
        t3 = Dict((r[ii], r[ki2]) => (r[di], r[fi]) for r in f["F3"].rows if r[si] == "T_s")
        for r in f["F3"].rows
            r[si] == "W2" && ((r[di], r[fi]) = t3[(r[ii], r[ki2])])
        end
        si, ni = colidx(f, "F1", "set"), colidx(f, "F1", "N")
        for r in f["F1"].rows
            if r[si] == "W2" && r[ni] > 30
                r[ni] = 30; r[5], r[6], r[7], r[8] = 20, 4, 2, 4
            end
        end
    end)
    push!(cases, "f5-single-k-index" => f -> begin
        ki = colidx(f, "F5", "K index")
        setrows!(f, "F5", filter(r -> r[ki] == 1, f["F5"].rows))
    end)
    push!(cases, "f6-class-count-changes-with-step" => f -> begin
        ci, si, ni = colidx(f, "F6", "class"), colidx(f, "F6", "step"), colidx(f, "F6", "count")
        for r in f["F6"].rows
            r[ci] == "R" && r[si] > 1000 && (r[ni] = 1100)
        end
    end)
    push!(cases, "f1-single-budget" => f -> begin
        bi = colidx(f, "F1", "budget")
        setrows!(f, "F1", filter(r -> r[bi] == 12800, f["F1"].rows))
    end)
    cases
end

section("FIXT") do
    for (name, mut) in fixture_cases()
        figs = basefigs()
        mut(figs)
        c = drawcase(name, resultlike(figs))
        extra = ""
        if name == "f6-class-count-changes-with-step" && "p2-indicators-fig-F6.pdf" in c.files
            txt = py("text", joinpath(c.out, "p2-indicators-fig-F6.pdf"))
            extra = "; F6 text with the changed count: " * join(unique([l for l in split(txt, '\n') if occursin("1100", l) || occursin("count", l) || occursin("changes", l) || occursin("varies", l)]), " / ")
        elseif name == "f1-single-budget" && "p2-indicators-fig-F1.pdf" in c.files
            txt = py("text", joinpath(c.out, "p2-indicators-fig-F1.pdf"))
            extra = "; F1 titles: " * join([l for l in split(txt, '\n') if occursin(" = ", l)], " / ")
        end
        probe = ""
        if c.err === nothing && length(c.sig) == 6
            tx(k) = strip.(split(py("text", joinpath(c.out, "p2-indicators-fig-F$k.pdf")), '\n'))
            if name == "k3-all-auc-nothing"
                probe = "; F4 spans 'AUC not computed:' = $(count(==("AUC not computed:"), tx(4))) (K = 0.0, 3.0, 10.0 expected 3)"
            elseif name == "f5-single-k-index"
                probe = "; F5 row titles 'K index' = $(count(startswith("K index"), tx(5))) (one per group expected 2)"
            elseif name == "no-rli-anywhere"
                probe = "; F1 spans 'RLI*' = $(count(==("RLI*"), tx(1))), F2 spans starting 'RLI' = $(count(startswith("RLI"), tx(2)))"
            elseif name == "no-rli-in-S2"
                probe = "; F2 RLI labels: $(join([l for l in tx(2) if startswith(l, "RLI (")], ", "))"
            elseif name == "min-size-set-30-orbits"
                probe = "; F2 titles: $(join([l for l in tx(2) if occursin(r"^(S2|W2|T_s) \(N =", l)], ", "))"
            end
        end
        emit("D91m", "fixture '$name': all six figures drawn from the fixture blocks: error; files with the %PDF- signature", "$(c.err === nothing ? "none" : c.err[1:min(end, 300)]); $(length(c.sig)) of 6 ($(join(setdiff(["p2-indicators-fig-F$k.pdf" for k in 1:6], c.sig), " ")) missing)$extra$probe", "no error; 6 files with signature", pf(c.err === nothing && length(c.sig) == 6))
    end
end

# ---------------------------------------------------------------------------------------------------------------- P1 (D91e, I24)
section("P1") do
    out = joinpath(SCR, "p1")
    r = runcli("p1", joinpath(ROOT, "verify"), out, args1())
    emit("D91e", "real P1 inputs under the verify project: exit code; wall time", "$(r.exit); $(round(r.wall; digits = 1)) s; output: " * replace(strip(r.text), '\n' => " / ")[1:min(end, 250)], "0", pf(r.exit == 0))
    have = isdir(out) ? sort(readdir(out)) : String[]
    expect = sort(vcat(["p1-le-benchmark-tables.md"], PDFNAMES1))
    emit("D91e", "P1 output folder: the tables file and the eleven figure files (F1 to F7, F1a, F1b, F2a, F2b)", join(have, " "), join(expect, " "), pf(have == expect))
    if isfile(joinpath(out, "p1-le-benchmark-tables.md"))
        d = tdiff(joinpath(out, "p1-le-benchmark-tables.md"), ACC1)
        emit("D91e", "P1 tables file against the tracked accepted file, line by line (lines in new and old; lines that differ)", "$(d.na) and $(d.nb); $(length(d.diffs)): " * join(["line $i: \"$(d.old[i])\" -> \"$(get(d.new, i, ""))\"" for i in d.diffs], "; "), "only module version and git commit", pf(d.ok))
        dn, do_ = digestlines(joinpath(out, "p1-le-benchmark-tables.md")), digestlines(ACC1)
        emit("D91e", "P1 digests F1 to F7 equal the accepted ones", "$(length(dn)) digests, equal $(dn == do_)", "7, equal", pf(length(dn) == 7 && dn == do_))
    end
    # signature of all eleven
    for l in split(strip(py("sig", [joinpath(out, n) for n in PDFNAMES1]...)), '\n')
        p = split(l, '\t')
        ok = length(p) >= 4 && lowercase(p[2]) == "true" && parse(Int, p[3]) > 1000 && strip(p[4]) == "1"
        emit("D91e", "P1 file $(p[1]): signature, size (bytes), pages, page size", join(p[2:end], ", "), "%PDF-, > 1000, 1 page", pf(ok))
    end
    # pixel equality with the tracked P1 files of the study folder (F1 to F7); renderer PyMuPDF (MuPDF), 72 dpi
    for l in split(strip(py("cmp", out, TRACKED1, PDFNAMES1[1:7]...)), '\n')
        p = split(l, '\t')
        length(p) >= 3 || continue
        nd = parse(Int, p[2])
        emit("D91e", "$(p[1]): render at 72 dpi (PyMuPDF 2.0, MuPDF 1.28.5) of the new file against the tracked file of workspaces/study/stages/04-analyze/output, differing pixels; $(p[3])", nd, "0", pf(nd == 0))
    end
    # the paper folder holds tracked copies of all eleven (F1a, F1b, F2a, F2b included): reported
    for l in split(strip(py("cmp", out, PAPERFIG, PDFNAMES1...)), '\n')
        p = split(l, '\t')
        length(p) >= 3 || continue
        emit("D91e", "$(p[1]): render at 72 dpi of the new file against the copy tracked in workspaces/paper/stages/04-finalize (the half-width variants have no file in the study folder), differing pixels; $(p[3])", parse(Int, p[2]), "0 (reported)", parse(Int, p[2]) == 0 ? "PASS" : "INFO")
    end
end

# the P1 pictures of the module before the amendment (the commit before the spec of D91), drawn from a copy of that tree: the cause of any difference with the tracked files
section("P1PRE") do
    pre = joinpath(SCR, "pre-d91")
    rm(pre; recursive = true, force = true); mkpath(pre)
    zipf = joinpath(SCR, "pre-d91.zip")
    run(`git -C $ROOT archive --format=zip -o $zipf fc4cace`)
    py("unzip", zipf, pre)
    rm(zipf; force = true)
    out = joinpath(SCR, "p1pre")
    r = runcli("p1pre", joinpath(pre, "verify"), out, args1(); bin = joinpath(pre, "bin", "chaoslab.jl"))
    emit("D91e", "P1 command of the tree of commit fc4cace (the module before D91; copy made with git archive): exit code; wall time (this run precompiles that copy)", "$(r.exit); $(round(r.wall; digits = 1)) s", "0", pf(r.exit == 0))
    cur = joinpath(SCR, "p1")
    if isdir(out) && isdir(cur)
        d = tdiff(joinpath(cur, "p1-le-benchmark-tables.md"), joinpath(out, "p1-le-benchmark-tables.md"))
        emit("D91e", "P1 tables file of the current module against that of the module before D91: lines that differ", "$(length(d.diffs)): " * join(["line $i" for i in d.diffs], ", "), "only module version and git commit", pf(d.ok))
        for l in strip.(split(strip(py("cmp", cur, out, PDFNAMES1...)), '\n'))
            p = split(l, '\t')
            length(p) >= 3 || continue
            nd = parse(Int, p[2])
            emit("D91e", "$(p[1]): render at 72 dpi of the current module's file against the file of the module before D91, differing pixels; $(p[3])", nd, "0", pf(nd == 0))
        end
        for n in ("p1-le-benchmark-fig-F1.pdf", "p1-le-benchmark-fig-F3.pdf")
            tr = strip(read(`git -C $ROOT log -1 --format=%h,%ad --date=short -- workspaces/study/stages/04-analyze/output/$n`, String))
            emit("D91e", "cause of the difference with the tracked study file $n: last commit that changed the tracked file (hash, date); the figure code was redrawn after it (caaf95a, 2026-10-21, owner visual check; 812d484)", tr, "reported", "INFO")
        end
    end
end

# ---------------------------------------------------------------------------------------------------------------- MISC (I27, I28, I29, D91l, D91n, D91o folders)
section("MISC") do
    run1 = joinpath(SCR, "run1")
    if isdir(run1)
        pyemit(py("enc", run1))
        hv = get(VERD, "D91h", String[])
        emit("I27", "no silent alteration: every row of a block is drawn (D91h: all rows of F2 to F6 and all cells of F1 and F5 re-read) and the gaps are empty and stated (F4 gaps, text 'AUC not computed'); the verdicts of the D91h re-read checks of the real run", "$(count(==("PASS"), hv)) PASS of $(length(hv))", "all PASS", pf(!isempty(hv) && all(==("PASS"), hv)))
    end
    gd = read(`git -C $ROOT diff --stat fc4cace HEAD -- Project.toml Manifest.toml verify/Project.toml verify/Manifest.toml`, String)
    emit("I29", "dependency files unchanged by the D91 commits (git diff --stat of Project.toml, Manifest.toml, verify/Project.toml, verify/Manifest.toml between the commit before the amendment, fc4cace, and HEAD)", isempty(strip(gd)) ? "no difference" : strip(gd), "no difference", pf(isempty(strip(gd))))
    pk = read(pipeline(ignorestatus(Cmd(String[JL, "--project=$ROOT", "-e", "using Chaoslab; print(sort([string(k.name) for k in keys(Base.loaded_modules) if occursin(r\"Makie|Cairo\", string(k.name))]))"])); stderr = devnull), String)
    emit("I29", "the analysis loads no plotting package: modules loaded by 'using Chaoslab' under the default project that match Makie or Cairo", pk, "[]", pf(strip(pk) == "String[]" || strip(pk) == "[]"))
    emit("D91n", "layering test and dependency list (package suite)", "not run here: criterion D91n is (S), covered by the package suite; I29 above checks the dependency files and the plotting package only from outside", "(S)", "INFO")
    emit("D91l", "owner visual check of the six figures of the real run", "human judgement: not measured; contact sheets written by the CONTACT section (analysis-d91-F1.png to F6.png)", "owner: accepted or the change asked for", "INFO")
    emit("D91e", "(S) part of D91e: the P1 fixture run of the base suite", "not run here: criterion part (S), package suite", "(S)", "INFO")
end

section("CONTACT") do
    isdir(joinpath(SCR, "run1")) && for l in strip.(split(strip(py("contact", joinpath(SCR, "run1"), OUT)), '\n'))
        emit("D91l", "contact sheet for the owner visual check (72 dpi render of the figure of the real run)", joinpath(OUT, l), "written", isfile(joinpath(OUT, l)) ? "INFO" : "FAIL")
    end
end

section("HASH1") do
    s0 = SNAP0[]
    if s0 === nothing
        emit("D91o", "fingerprints after the runs", "no snapshot from the HASH0 section in this invocation", "-", "INFO")
    else
        t = @timed snapshot()
        s1 = t.value
        dd = [f for f in union(keys(s0.data), keys(s1.data)) if get(s0.data, f, nothing) != get(s1.data, f, nothing)]
        od = [f for f in union(keys(s0.other), keys(s1.other)) if get(s0.other, f, nothing) != get(s1.other, f, nothing)]
        emit("D91o", "files under data/p2-indicators and data/p1-le-benchmark after the runs (SHA-256 each) against before", "$(length(s1.data)) files, $(length(dd)) differ", "0 differ, same files", pf(isempty(dd) && length(s0.data) == length(s1.data)))
        emit("D91o", "the other inputs (preregistrations, manifests, panel, parameter tables, accepted tables files, tracked P1 figure files, defaults) after the runs", "$(length(s1.other)) files, $(length(od)) differ", "0 differ", pf(isempty(od) && length(s0.other) == length(s1.other)))
        emit("D91o", "listing of the whole data tree (path, size, mtime) after the runs: entries and SHA-256 against before", "$(s1.ntree) vs $(s0.ntree); $(s1.tree[1:16]) vs $(s0.tree[1:16])", "same", pf(s1.tree == s0.tree))
    end
    for (d, names) in (("run1", PDFNAMES2), ("run2", PDFNAMES2), ("api1", PDFNAMES2), ("p1", PDFNAMES1))
        p = joinpath(SCR, d)
        isdir(p) || continue
        tabs = d == "p1" ? "p1-le-benchmark-tables.md" : "p2-indicators-tables.md"
        have = sort(readdir(p))
        emit("D91o", "output folder $d holds only the tables file and the figure files of the study", join(have, " "), join(sort(vcat([tabs], names)), " "), pf(have == sort(vcat([tabs], names))))
    end
end

say(@sprintf("TOTAL | PASS %d | FAIL %d | INFO %d | %s", COUNTS["PASS"], COUNTS["FAIL"], COUNTS["INFO"], Dates.format(now(), "yyyy-mm-dd HH:MM")))
say("done")
close(LOG)
