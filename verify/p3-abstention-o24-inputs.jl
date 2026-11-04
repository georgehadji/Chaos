# Study p3-abstention, gate O24 (round-loop.md): the fixture inputs of the dry campaign. The panel input block and run block are read from section
# "Amendment A2" of the committed preregistration and shrunk (panel seed 777, 4 K, 3 rounds; horizons and windows to the toy p2 of the test helpers:
# N_1 80000, N_c 20480, chk 1000; k19 threads 4; pilot horizon 1000); arms, specs, release keys, health and budgets are those of A2.
# Writes <out dir>/inputs (test-side inputs, build copy X and sealed copy Y, two-file seal), p1/, p2/ and p3-abstention-test-panel.toml.
# Run from the repo root:  julia --project -t 4 verify/p3-abstention-o24-inputs.jl <out dir>
# Then:  julia --project -t 4 verify/experiments-d96.jl lib <out dir> <run dir>;  bash verify/experiments-d96-cli.sh <out dir> <run dir>;
#        julia --project -t 4 verify/experiments-d96.jl analyze <out dir> <run dir> lib
using Test, Chaoslab, StableRNGs, TOML
const REPO = pwd()
const WORK = abspath(ARGS[1])
occursin("/data/", replace(WORK, "\\" => "/")) && error("refused: $WORK is inside data/")
const TA = joinpath(REPO, "test", "analysis")
for f in ("an_helpers.jl", "an_p2_ref.jl", "an_p2_helpers.jl", "an_p2_study.jl", "an_p2_build.jl", "an_p2_refstore.jl", "an_p3_ref.jl", "an_p3_helpers.jl", "an_p3t_ref.jl", "an_p3t_helpers.jl")
    include(joinpath(TA, f))
end
const E = Chaoslab.Experiments
mkpath(WORK)
# A2 blocks read from the committed preregistration (e5d97f8), section "Amendment A2"
pre = replace(read(joinpath(REPO, "workspaces/study/stages/01-design/output/p3-abstention-preregistration.md"), String), (Char(13) * Char(10)) => Char(10))
a2 = pre[findfirst("## Amendment A2 (2026-11-04)", pre)[1]:end]
blocks = [m.captures[1] for m in eachmatch(r"```toml\n(.*?)```"s, a2)]
@assert length(blocks) == 2
pan, run = blocks
rep(t, a, b) = (n = count(a, t); n >= 1 || error("no match $a"); replace(t, a => b))
# ---- tiny panel block (same keys, same structure)
pan = rep(pan, "K = [\"0\", \"0.5\", \"0.971635406\", \"pi/2\", \"2pi\", \"10\"]", "K = [\"0\", \"0.5\", \"0.971635406\", \"2pi\"]")
pan = rep(pan, "N_scr = 100000", "N_scr = 30000")
pan = rep(pan, "seed = 20261030", "seed = 777")
pan = rep(pan, "screen_points = 262144", "screen_points = 400")
pan = rep(pan, "r_max = 5", "r_max = 3")
pan = rep(pan, "db_size = 100", "db_size = 5")
pan = replace(pan, r"\[panel\.quota\].*"s => """
[panel.quota]
"1" = { T = 1, R = 1 }
"2" = { T = 2, R = 1, C = 2, S = 1, X = 1 }
"3" = { T = 2, R = 1, C = 2, S = 2, X = 1 }
"4" = { T = 1, R = 1, C = 1, S = 1, X = 1 }

[panel.d_fraction]
"1" = { T = 0.5, R = 0.5 }
"2" = { T = 0.5, R = 0.5, C = 0.5, S = 0.5, X = 0.5 }
"3" = { T = 0.5, R = 0.5, C = 0.5, S = 0.5, X = 0.5 }
"4" = { T = 0.5, R = 0.5, C = 0.5, S = 0.5, X = 0.5 }
""")
# ---- tiny run block (arms, specs, release keys identical; horizons and windows shrunk to the toy p2 values N_1 = 80000, N_c = 20480)
win = "[20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240]"
run = rep(run, "[100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]", win)
run = rep(run, "horizons = [10000000]", "horizons = [80000]")
run = rep(run, "horizons = [102400]", "horizons = [20480]")
run = rep(run, "horizons = [102401]", "horizons = [20481]")
run = rep(run, "horizons = [1000000]", "horizons = [1000]")   # chk horizon = N_B of the P1 table (real: 1e6 = N_B, with N_1 = 1e7)
run = rep(run, "threads = 23", "threads = 4")
run = rep(run, "horizon = 100000", "horizon = 1000")
a2x = "## Amendment A2 (fixture, o24 dry campaign)\n\nFixture copy of the drafted A2 blocks with tiny sizes. Not a study.\n\n```toml\n" * pan * "```\n\n```toml\n" * run * "```\n\n"
write(joinpath(WORK, "a2-fixture-section.md"), a2x)

# ---- redefine the preregistration builder: the same text as p3t_prereg_text with the A2 section before the A1 section
function p3t_prereg_text(q, t; slug = P3_SLUG, hashrow = true, a1 = true, texts = Dict{Symbol, String}())
    qq = p3_quotes(q); tq = p3t_quotes(t; texts)
    body = IOBuffer()
    println(body, "# Preregistration: ", slug, "\n\n## Header\n")
    println(body, "| Field | Value |\n|-------|-------|\n| Study slug | ", slug, " |")
    hashrow && println(body, "| SHA-256 of this file (current version) | ", "0"^64, " (rule as p2: fixture) |")
    println(body, "| Amendment log | A1 |\n\nFixture. Not a study.\n")
    for sk in P3_SEC_ORDER
        println(body, "## ", getfield(P3_SEC, sk), "\n")
        for k in P3_PARAM_KEYS
            qq[k][1] == getfield(P3_SEC, sk) && println(body, "- ", qq[k][2])
        end
        if getfield(P3_SEC, sk) == P3_SEC.COST
            for k in P3T_PILOT; println(body, "- ", tq[k][2]); end
        end
        println(body)
    end
    print(body, a2x)
    if a1
        println(body, "## ", P3T_A1, "\n")
        for k in P3T_KEYS
            k in P3T_PILOT || println(body, "- ", tq[k][2])
        end
        println(body)
    end
    text = String(take!(body))
    prep(s) = replace(replace(s, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
    h = p3ref_sha256(prep(text))
    hashrow && (text = replace(text, "0"^64 => h; count = 1))
    return text, h
end

st = p3_cached("o24dry"; p = p2_params(; n_s = 1000))   # n_s = 1000 lies on the recording grid (20 per decade); the real study has n_s = 1e5
I = p3t_inputs((q = p3t_q(; n_s = 1000),), st, "o24dry")
inp = joinpath(WORK, "inputs"); rm(inp; recursive = true, force = true); cp(I.dir, inp)
for (nm, d) in (("p2", st.s2.dir), ("p1", st.p1.dir)); end
# the p2 and P1 inputs the analysis reads
p2 = joinpath(WORK, "p2"); mkpath(p2)
for f in (st.s2.manifest, st.s2.prereg, st.s2.params); cp(f, joinpath(p2, basename(f)); force = true); end
p1 = joinpath(WORK, "p1"); mkpath(p1)
for f in (st.p1.prereg, st.p1.params); cp(f, joinpath(p1, basename(f)); force = true); end
X = joinpath(inp, "p3-abstention-preregistration.md")
chk = E.check_prereg(X)
println("X header ", chk.header_sha256, " computed ", chk.computed_sha256, " amendments ", chk.amendments)
@assert chk.header_sha256 == chk.computed_sha256
@assert chk.header_sha256 == I.prereg_hash
# panel from X (build copy)
PANEL = joinpath(WORK, "p3-abstention-test-panel.toml")
t = isfile(PANEL) ? 0.0 : @elapsed E.build_panel(X; out = PANEL)   # reuse the panel of a previous run only if the inputs are byte-identical (checked below)
println("panel built in ", t, " s")
# seal: Y = X with the accepted list replaced by H_X and the header recomputed
tx = read(X, String)
ty = replace(tx, "0"^64 * "\"" => chk.header_sha256 * "\""; count = 1)   # list entry (first quoted run of zeros after the header cell)
@assert ty != tx
hx = chk.header_sha256
ty = replace(ty, hx => E.prereg_hash(ty); count = 1)
Y = joinpath(inp, "p3-abstention-preregistration.sealed.md")
write(Y, ty)
cy = E.check_prereg(Y)
@assert cy.header_sha256 == cy.computed_sha256
println("Y header ", cy.header_sha256, " accepted ", E.read_run_block(Y).raw["panel_prereg_sha256_accepted"])
open(joinpath(WORK, "inputs.txt"), "w") do io
    println(io, "X_hash=", hx); println(io, "Y_hash=", cy.header_sha256)
end
pn = TOML.parsefile(joinpath(WORK, "p3-abstention-test-panel.toml"))
println(pn["manifest"])
for o in pn["orbit"]; println(o["id"], " k", o["k_index"], " ", o["class"], " r", get(o, "round", "-"), " D=", get(o, "in_D", "-"), " K=", o["K"]); end
