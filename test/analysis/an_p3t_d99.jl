# Amendment D99, stage 02: corrections of the p3 test side and of the calibrate default. Spec: analysis-spec.md, Amendment D99, E83 to E86, I52, I53, O31, O32, criteria D99a to D99i.
# Helpers of an_p3t_helpers.jl and of an_p3t_d97.jl (the d97_* readers of the tables file) are reused; this file adds d99_* only. The D99 testsets are RED because the text, one status cell, one
# target cell, one added line or the default folder is the old one, never because of a fixture. Guards that pass before the module changes are named "(guard; passes now)".
# D99g (the chain with defaults on the campaign of the runner) and D99h (O31, before and after on the committed inputs) need the real store: stage 04, see the test plan.
# No random draw is made here (the toy campaign draws with p3t_seed('a')). Hand values: the fixtures of D99a are written in the spec (delta_a = 0.1).

# a number as the tables print it: a decimal, or a rational a/b
d99_num(s) = (m = match(r"^(-?\d+)/(\d+)$", s); m !== nothing ? parse(Int, m.captures[1]) / parse(Int, m.captures[2]) : tryparse(Float64, s))

# the rows of table T5 (cells stripped) and the decision line (the line that begins `HA:`)
d99_t5(tab) = (t = d97_tables(d97_sec(tab, 5)); isempty(t) ? Vector{Vector{String}}() : first(t)[2:end])
d99_ha_line(tab) = (l = [l for l in split(d97_sec(tab, 5), '\n') if startswith(l, "HA:")]; isempty(l) ? "" : first(l))
d99_i3_fails(line) = occursin(r"item 3\W{0,3}(?:\w+\W+){0,2}fails", line)         # "item 3 fails", "item 3: fails", "item 3 as a whole fails"

# the own condition of a row of item 3 from the cells it prints (E83): point at least delta_a and the lower bound above 0
function d99_own(row)
    pt = d99_num(row[3]); m = match(r"lower bound (\S+) ", row[4]); da = d99_num(replace(row[5], "delta_a = " => ""))
    lb = m === nothing ? nothing : d99_num(m.captures[1])
    return (pt !== nothing && lb !== nothing && da !== nothing && pt >= da && lb > 0) ? "ok" : "fails"
end

# a hand-built analysis state for the T5 writer (O32): only the fields that T5 prints, with delta_a = 0.1 and no overlap rows. Items 1, 2, 4 are given by `failed`.
function d99_ha(; decision, failed, dw, dw_lo, dt, dt_lo, e = 0.05, e_up = 0.08, g = 0.0)
    b(v) = (; value = v, dropped = 0)
    point = (; e, aE = 0.2, dw, dt, g, a = Dict(:E => 0.2, :W2 => 0.1, :Ts => 0.05, :S2 => 0.3))
    bounds = (; e_up = b(e_up), aE_up = b(0.25), dw_lo = b(dw_lo), dt_lo = b(dt_lo))
    return (; ha = (; decision, failed, point, bounds, strata = [5, 5, 7, 1, 7]), qt = (; e_max = 0.1, a_max = 0.3, delta_a = 0.1, g_max = 0.05), O = (; rows = NamedTuple[]))
end
d99_t5_fixture(A) = d97_call(() -> join(AN._p3t_t5(A), '\n'))

# the four fixtures of D99a: name => (state, statuses of the two rows of item 3 in the order of the table, statuses of the rows of items 1, 2, 4, decision, item 3 as a whole)
const D99_FIXTURES = [
    ("1 first row fails, second ok", d99_ha(decision = "not supported", failed = [3], dw = 0.0527, dw_lo = 0.043, dt = 0.2069, dt_lo = 0.158), ["fails", "ok"], ["ok", "ok", "ok"], "not supported", false),
    ("2 both ok", d99_ha(decision = "supported", failed = Int[], dw = 0.12, dw_lo = 0.05, dt = 0.2069, dt_lo = 0.158), ["ok", "ok"], ["ok", "ok", "ok"], "supported", true),
    ("3 both points below delta_a (item 1 fails too)", d99_ha(decision = "not supported", failed = [1, 3], dw = 0.05, dw_lo = 0.01, dt = 0.08, dt_lo = 0.02, e = 0.15, e_up = 0.2), ["fails", "fails"], ["fails", "ok", "ok"], "not supported", false),
    ("4 point 0.15 with lower bound -0.01 (item 4 fails too)", d99_ha(decision = "not supported", failed = [3, 4], dw = 0.15, dw_lo = -0.01, dt = 0.2069, dt_lo = 0.158, g = 0.5), ["fails", "ok"], ["ok", "ok", "fails"], "not supported", false)]

# ---------------------------------------------------------------- D99a and D99b: T5 (E83, O32)
@testset "D99a T5: each row of item 3 shows the status of its own condition; the decision line states item 3 as a whole; the decision and the numbers are unchanged (hand fixtures, O32)" begin
    for (name, A, rows3, rows124, dec, whole) in D99_FIXTURES
        txt = d99_t5_fixture(A)
        @testset "fixture $name" begin
            @test d97_ok(txt)
            d97_ok(txt) || continue
            rows = [r for t in d97_tables(txt) for r in t[2:end] if length(r) == 6 && r[1] in ("1", "2", "3", "4")]
            @test [r[1] for r in rows] == ["1", "2", "3", "3", "4"] && [r[2] for r in rows][3:4] == ["a(W2) - a(S2)", "a(T_s) - a(S2)"]
            @test [r[6] for r in rows[3:4]] == rows3                                                      # D99a: the status of the row, not of the item
            @test [r[6] for r in rows[[1, 2, 5]]] == rows124                                              # D99b: items 1, 2, 4 as before
            ha = first(l for l in split(txt, '\n') if startswith(l, "HA:"))
            @test startswith(ha, "HA: $dec.") && d99_i3_fails(ha) == !whole                               # the decision is the hand decision; item 3 as a whole in the decision line
            @test [r[3] for r in rows] == [string(A.ha.point.e), string(A.ha.point.aE), string(A.ha.point.dw), string(A.ha.point.dt), string(A.ha.point.g)] &&
                  [r[5] for r in rows] == ["e_max = 0.1", "a_max = 0.3", "delta_a = 0.1", "delta_a = 0.1", "g_max = 0.05"]      # points and thresholds unchanged
            @test occursin("status of a row of item 3", lowercase(txt)) && occursin("own condition", lowercase(txt))                      # E83: the added sentence
        end
    end
end

@testset "D99a/D99b the T5 of the run: every row of item 3 shows its own condition from its printed point, bound and delta_a; items 1, 2, 4 follow the failed list (guard on the toy run)" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; m, tab) = R
    rows = d99_t5(tab)
    failed = Int.(m["decisions"]["HA"]["failed"])
    @test length(rows) == 5 && [r[1] for r in rows] == ["1", "2", "3", "3", "4"]
    for r in rows
        r[1] in ("1", "2", "4") && @test r[6] == (parse(Int, r[1]) in failed ? "fails" : "ok")
        r[1] == "3" && @test r[6] == d99_own(r)
    end
    @test d99_i3_fails(d99_ha_line(tab)) == (3 in failed) && startswith(d99_ha_line(tab), "HA: " * m["decisions"]["HA"]["verdict"])
end

# ---------------------------------------------------------------- D99c: the omission text names no study amendment (E84)
@testset "D99c P120 and the T3 introduction: `deviation of the analysis`, the cause, no `A2`, no `study amendment`, in the tables file and in the machine file" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; text, tab) = R
    p120_tab = [l for l in split(tab, '\n') if occursin("P120", l)]
    p120_json = [l for l in split(text, '\n') if occursin("P120", l)]
    sec3 = d97_sec(tab, 3)
    intro = join([l for l in split(sec3, '\n') if !startswith(l, "|") && !startswith(l, "## ") && !startswith(l, "Table T") && !isempty(strip(l)) && !startswith(l, "Orbits on which")], ' ')       # the introduction of T3, not the count sentence
    @test length(p120_tab) == 1 && length(p120_json) == 1 && !isempty(intro) && occursin("labels agree", intro)
    for (what, s) in (("P120 in the tables file", join(p120_tab)), ("P120 in the machine file", join(p120_json)), ("the T3 introduction", intro))
        @testset "$what" begin
            @test occursin("deviation of the analysis", s)
            @test all(w -> occursin(w, s), ["chk", "N_B", "three precisions"])                                 # the cause of E84
            @test !occursin("A2", s) && !occursin("study amendment", lowercase(s))
        end
    end
    @test occursin("BigFloat(256)", join(p120_tab)) && occursin("BigFloat(256)", join(p120_json))              # the precision label of the cause, in the deviation
    @test !occursin(r"\bA2\b", text) && !occursin("study amendment", lowercase(text)) && !occursin(r"\bA2\b", tab) && !occursin("study amendment", lowercase(tab))   # 0 matches in both files
end

# ---------------------------------------------------------------- D99d: T2 target of D (E85)
@testset "D99d T2: the row D has `n.a.` as target; the six other rows keep their targets" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; tab, I) = R
    t2 = d97_tables(d97_sec(tab, 2))
    set = [t for t in t2 if first(t) == ["set", "size", "target"]]
    @test length(set) == 1
    isempty(set) && return
    rows = Dict(r[1] => r for r in first(set)[2:end])
    @test [r[1] for r in first(set)[2:end]] == ["E_test", "W2", "T_s", "S2", "D", "D_B", "gate"]
    @test rows["D"][3] == "n.a." && rows["D"][3] != string(I.t.target_f)                                         # 45 in the study is n_f, not a size target of D
    @test (rows["W2"][3], rows["T_s"][3]) == (string(I.t.target_w2), string(I.t.target_t_s)) == ("11", "5")        # hand: the toy table of the fixture (p3t_inputs)
    @test all(k -> rows[k][3] == "n.a.", ["E_test", "S2", "D_B", "gate"])                                       # the sets without an orbit target, as before
end

# ---------------------------------------------------------------- D99e: T8 added line (E85)
function d99_nf()
    get!(D97_CACHE, "d99nf") do
        S, st = p3t_cached(), p3_cached("d95")
        I = p3t_inputs(S, st, "d99nf"; target_f = 37)                  # a value that no count of the toy takes: the line cannot match by chance
        inst = p3t_instance(S, I, "d99nf")
        write(inst.log, join(first(p3t_log_entries(S, I))))
        res = d97_run(S, I, st, inst)
        (; S, I, inst, res, X = p3t_expected(S, I.qt), d97_read(res)...)
    end
end

function d99_t8_check(m, tab, nf, X)
    sec = d97_sec(tab, 8)
    lines = split(sec, '\n')
    k = findfirst(l -> startswith(l, "HC:"), lines)
    k === nothing && return false
    opening = lines[k]
    tbl = findfirst(l -> startswith(l, "|"), lines)
    added = [l for l in lines[(k + 1):(tbl === nothing ? end : tbl - 1)] if startswith(l, "Target failures n_f:")]        # after the opening sentence, before the first table
    fl = parse(Int, match(r"failures (\d+),", opening).captures[1])
    @test fl == X.hc.failures == Int(m["hc"]["failures"])
    @test length(added) == 1
    isempty(added) && return false
    mm = match(r"^Target failures n_f: (\S+); achieved failures: (\S+)$", only(added))
    @test mm !== nothing
    mm === nothing && return false
    @test d99_num(mm.captures[1]) == nf && d99_num(mm.captures[2]) == fl && occursin(r"underpowered: " * string(fl < nf), opening)       # the flag compares failures with n_f, unchanged
    @test m["hc"]["underpowered"] == (fl < nf) == X.hc.underpowered
    # the machine file carries the same two numbers in the data of T8: the achieved failures under `failures`, n_f under a key that names a target (the spec fixes no key name, see the plan)
    tk = [k for (k, v) in m["hc"] if occursin("target", lowercase(k)) && v isa Real && v == nf]
    @test length(tk) == 1
    return true
end

@testset "D99e T8 has the line `Target failures n_f: <n_f>; achieved failures: <failures>`: values equal the parameter table, the opening sentence and the machine file (main toy run, n_f = 6)" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; m, tab, I, X) = R
    @test I.t.target_f == 6
    d99_t8_check(m, tab, 6, X)
end

@testset "D99e the same with n_f = 37 (the failure target is read from the parameter table, not from the count; the underpowered flag is true)" begin
    R = d99_nf()
    @test d97_ok(R.res)
    d97_ok(R.res) || return
    @test R.I.t.target_f == 37 && R.m["hc"]["underpowered"] == true
    d99_t8_check(R.m, R.tab, 37, R.X)
end

# ---------------------------------------------------------------- D99f: the default folder of calibrate (E86, O32)
# A hand-built layout of the study: <root>/stages/01-design/output/<prereg> (a byte copy of the preregistration of the toy study) and the stage 04 folder next to it. The path of the stage 04
# folder is written by hand here, not computed by the module. The inputs of the study stay where they are and are given by option (p3_run).
function d99_layout(tag, st)
    root = an_path("d99f-$tag"); rm(root; recursive = true, force = true)
    s01 = joinpath(root, "stages", "01-design", "output"); mkpath(s01)
    prereg = joinpath(s01, basename(st.prereg)); cp(st.prereg, prereg)
    return (; root, s01, s04 = joinpath(root, "stages", "04-analyze", "output"), prereg)
end
d99_same(a, b) = normpath(abspath(a)) == normpath(abspath(b))

@testset "D99f calibrate without --out writes the record and the tables to the folder from which count-round and analyze read the record" begin
    st = p3_cached("d95")
    L = d99_layout("default", st)
    r = d97_call(() -> p3_run(st; prereg = L.prereg, out = nothing))
    @test d97_ok(r) && r.status === :ok
    d97_ok(r) || return
    kw = d97_call(() -> AN._p3t_keywords(L.prereg, Dict{String, String}()))           # the readers of count-round and analyze (one keyword function, cli.jl)
    @test d97_ok(kw) && d99_same(dirname(kw.calibration_record), L.s04)                  # hand path: stages/04-analyze/output
    @test d99_same(dirname(r.record), L.s04) && d99_same(dirname(r.tables), L.s04)
    d97_ok(kw) && @test d99_same(dirname(r.record), dirname(kw.calibration_record)) && basename(r.record) == basename(kw.calibration_record)      # string equality of the folder, the file name is the one read
    @test isfile(r.record) && isfile(r.tables) && sort(readdir(L.s01)) == [basename(L.prereg)]                           # nothing is written to the stage 01 folder
    # the command line without --out: the same folder
    L2 = d99_layout("cli", st)
    io = IOBuffer(); args = vcat(["calibrate", L2.prereg], p3_cli_args(st)[3:end]); k = findfirst(==("--out"), args)
    deleteat!(args, k:(k + 1))
    @test d97_call(() -> AN.cli(args; io = io)) == 0 && isfile(joinpath(L2.s04, basename(r.record))) && isfile(joinpath(L2.s04, basename(r.tables))) && sort(readdir(L2.s01)) == [basename(L2.prereg)]
end

@testset "D99f an explicit output folder is used as given, nowhere else is written, and the record has the same bytes wherever it is written (guard; passes now for the folder, I53)" begin
    st = p3_cached("d95")
    L = d99_layout("explicit", st)
    out = joinpath(L.root, "elsewhere")
    r = d97_call(() -> p3_run(st; prereg = L.prereg, out = out))
    @test d97_ok(r) && r.status === :ok
    d97_ok(r) || return
    @test d99_same(dirname(r.record), out) && d99_same(dirname(r.tables), out) && sort(readdir(out)) == sort([basename(r.record), basename(r.tables)])
    @test sort(readdir(L.s01)) == [basename(L.prereg)] && !ispath(joinpath(L.root, "stages", "04-analyze"))
    L2 = d99_layout("bytes", st)
    r2 = d97_call(() -> p3_run(st; prereg = L2.prereg, out = nothing))
    @test d97_ok(r2) && isfile(r2.record) && read(r.record) == read(r2.record)             # same bytes in the default folder and in an explicit one
end
