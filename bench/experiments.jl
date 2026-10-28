using BenchmarkTools, Chaoslab, TOML

const SUITE = BenchmarkGroup()
const EXP = Chaoslab.Experiments

# One screening point at the preregistered size order (N_scr = 1e5, w = 1000): ns per orbit-step is time / 1e5.
const BUF = screen_buffer(100_000, 1000)
SUITE["screen_point! 1e5 steps K=5"] = @benchmarkable screen_point!($BUF, 5.0, 1.25, 2.5, 100_000, 1000, $(pi / 4), 1)
SUITE["screen_point! 2000 steps K=1.5"] = @benchmarkable screen_point($(1.5), 1.25, 2.5; N_scr = 2000, w = 100, alpha = $(pi / 4))

# One Newton solve at 256 bits (period 2, K = 5, seed at a grid centre).
SUITE["newton solve q=2 256 bit"] = @benchmarkable EXP._newton_seed(K, x, p, 2, 1e-12, 50, 1e-8) setup = (
    K = BigFloat(5.0); x = BigFloat(1.25); p = BigFloat(2.5)
)

# One small panel build (3 K values, 8 x 8 grid, N_scr = 2000, q_max = 2).
const BLOCK = Dict{String, Any}("K" => ["0", "1.5", "5"], "n_c" => 3, "grid_side" => 8, "N_scr" => 2000, "alpha" => "pi/4",
    "theta_T" => 0.01, "g_T" => 5.0, "r_T" => 0.01, "theta_C" => 0.05, "r_C" => 0.05, "theta_S" => 0.02, "w" => 100,
    "f_S" => 0.2, "L_S" => 500, "log10_r" => [-4.0, -2.0], "q_max" => 2, "rho" => 1e-10, "delta_R" => 1e-3, "newton_max" => 50,
    "delta_dup" => 1e-8, "n_ctrl_max" => 4, "mu" => "(3 + sqrt(5))/2", "s_D" => 2, "cap_D" => 100, "s_L" => 1, "cap_L" => 50)

function bench_prereg(dir)
    io = IOBuffer()
    print(io, "# Preregistration: bench fixture\n\n| Field | Value |\n|-------|-------|\n| SHA-256 of this file (current version) | ", "0"^64,
          " (rule: fixture) |\n| Amendment log | A1 (2026-01-01, fixture) |\n\n```toml\n")
    TOML.print(io, Dict("panel" => BLOCK))
    print(io, "```\n")
    text = String(take!(io))
    path = joinpath(dir, "bench-prereg.md")
    write(path, replace(text, "0"^64 => prereg_hash(text); count = 1))
    return path
end

const DIR = mktempdir()
const PREREG = bench_prereg(DIR)
SUITE["build_panel small (8x8, 3 K, q_max 2)"] = @benchmarkable build_panel($PREREG; out = $(joinpath(DIR, "panel.toml"))) samples = 3 evals = 1

# Study runner (D53): one work item per method family at N = 1e5 steps on the standard map; ns per orbit-step is time / 1e5.
# The A4 budget lines are in preregistration Amendment A4 ([run.budget]); the numbers measured here go in the implementation notes.
using StaticArrays

const RUN_N = 100_000
const RUN_GRID = record_grid(RUN_N, 1; points_per_decade = 20, first_steps = 100, extra_fractions = [0.1, 0.5])

function run_item_bench(method)
    sys, x0 = StandardMap(1.5), SVector(1.0, 2.0)
    k = method in ("L1", "L2") ? 1 : 2
    W0 = deviation_frame(Float64, "pi/4", k)
    b = EXP._new_buf(RUN_GRID, k, false)
    if method == "L2"
        u0 = SVector(W0[1, 1], W0[2, 1])
        return @benchmarkable EXP._record_shadow!($b, $sys, $x0, $u0, 1.0e-8, 1, RUN_N, 1) samples = 5 evals = 1
    end
    reducers = (EXP._estimator(method, 8.0), EXP.RunMonitor(sys))
    st0 = EXP.Tangent._init_states(reducers, Float64, Val(2), Val(k))
    return @benchmarkable EXP._record_loop!($b, $sys, $x0, $W0, $reducers, $st0, 1, RUN_N, 1) samples = 5 evals = 1
end

for (method, family) in (("L1", "Benettin"), ("L2", "two-trajectory"), ("L5", "QR family (Householder)"), ("L6", "QR library (LAPACK)"), ("L9", "SVD"))
    SUITE["run item $method ($family), 1e5 steps"] = run_item_bench(method)
end

# A small campaign chunk: the panel of the bench fixture, one S_L orbit per cell, L5 at 1e4 steps, written to a store.
const PANEL_HASH = prereg_hash(read(PREREG, String))
const RUN_BLOCK = """
[run]
panel_prereg_sha256_accepted = ["$PANEL_HASH"]

[run.primary]
tau = 1
eps = 1.0e-8
direction = "pi/4"
coordinate_mode = "lifted"

[[run.arm]]
name = "bench"
orbits = "S_L"
methods = ["L5"]
precisions = ["Float64"]
horizons = [10000]
h1_check = true

[run.grid]
points_per_decade = 20
first_steps = 100
extra_fractions = [0.1, 0.5]

[run.backend]
Float64 = "cpu"

[run.health]
fail_fraction = 0.05
control_bound = 0.01
control_methods = ["L5"]
setup_check_steps = { T = 5, R = 5, control_elliptic = 5, control_hyperbolic = 5, C = 5, S = 5, X = 5, h4 = 5 }
setup_check_tol = 1.0e-9
drift_bound = 1000000.0
rerun_attempts = 1
stop_failed_fraction = 0.05
first_chunk_cost_factor = 1000.0
wall_cap_hours = 1.0

[run.budget.Float64]
ns_cpu_thread = 100.0
safety = 2.0

[run.k19]
threads = $(Threads.nthreads())

[run.verify]
backend_agreement = { per_orbit_sigma = 6.0, mean_standard_errors = 4.0 }
rerun_sample = 2

[run.pilot]
orbits_per_k = 1
horizon = 1000
seed_label = "bench"
"""

function bench_run_prereg(dir)
    text = replace(replace(read(PREREG, String), "\r\n" => "\n"), "|-------|-------|\n" => "|-------|-------|\n| Study slug | bench |\n"; count = 1)
    text *= "\n## Amendment A4\n\n```toml\n" * RUN_BLOCK * "```\n"
    path = joinpath(dir, "bench-run-prereg.md")
    row = r"(?m)^(\| SHA-256 of this file[^|]*\|)[^|\n]*\|"
    write(path, replace(text, row => m -> match(r"^\| SHA-256 of this file[^|]*\|", m).match * " " * prereg_hash(text) * " |"))
    return path
end

const RUN_PREREG = bench_run_prereg(DIR)
const RUN_PANEL = joinpath(DIR, "run-panel.toml")
build_panel(PREREG; out = RUN_PANEL)
const RUN_STORE = joinpath(DIR, "run-store")
SUITE["run_campaign chunk (S_L, L5, 1e4 steps, written)"] = @benchmarkable run_campaign($RUN_PREREG, $RUN_PANEL; store = $RUN_STORE,
    manifest = $(joinpath(DIR, "run-manifest.json")), threads = Threads.nthreads(), allow_dirty = true, code_version = "bench") setup = (rm($RUN_STORE; recursive = true, force = true)) samples = 3 evals = 1

# Indicator recorders (D86): one orbit of the standard map per kind at N = 1e5 steps (the window kinds at windows 1000 and 4000, run length 8000 or 8001
# steps); ns per orbit-step is time / steps. The budget lines of the study are the study's; the numbers measured here go in the implementation notes.
const IND_SPECS = [
    ("tangent FLI MEGNO SALI", Dict{String, Any}("kind" => "tangent", "id" => "tg", "members" => ["FLI", "MEGNO", "SALI"])),
    ("tangent SALI", Dict{String, Any}("kind" => "tangent", "id" => "ts", "members" => ["SALI"])),
    ("stretching 1000 bins", Dict{String, Any}("kind" => "stretching", "id" => "st", "range" => "derived", "bins" => 1000, "interval" => 1, "min_per_cell" => 1, "counts" => "none")),
    ("frequency 1000 4000 (orbit pass)", Dict{String, Any}("kind" => "frequency", "id" => "fq", "observable" => "q_plus_i_p", "windows" => [1000, 4000], "order" => 1,
        "positions" => "half", "pad" => 4, "mean_removal" => true, "mirror_ratio" => 0.5, "floor_multiplier" => 8.0)),
    ("phase advance 1000 4000", Dict{String, Any}("kind" => "phase_advance", "id" => "pa", "observable" => "planar", "windows" => [1000, 4000], "birkhoff_order" => 2)),
    ("birkhoff K = 4, 1000 4000", Dict{String, Any}("kind" => "birkhoff", "id" => "bk", "windows" => [1000, 4000], "birkhoff_order" => 2,
        "observables" => [Dict{String, Any}("kind" => k, "coordinate" => 1, "frequency" => f) for k in ("cos", "sin") for f in (1, 2)])),
]

function ind_item_bench(spec)
    item = (system = "standard", param = "1.5", x0 = "1.0", p0 = "2.0", coordinate_mode = "wrapped")
    sys, x0 = WrappedStandardMap(1.5), torus_state(SVector(1.0, 2.0))
    red = EXP.ind_chain(Float64, item, spec, sys)
    L = EXP.ind_run_length(spec)
    N = L === nothing ? RUN_N : L
    b = EXP._new_ind_buf(record_grid(N, 1; points_per_decade = 20, first_steps = 100, extra_fractions = [0.1, 0.5]), true, length(EXP.ind_series_columns(spec)))
    if EXP.ind_window_kind(spec["kind"])
        W1 = SMatrix{2, 1, Float64}(1.0, 0.0)
        buf = Vector{ComplexF64}(undef, 2 * maximum(spec["windows"]))
        feed = spec["kind"] == "frequency" ? (j, x) -> (buf[j] = ComplexF64(x[1], x[2]); nothing) : (j, x) -> nothing
        st0 = EXP.Tangent._init_states(red, Float64, Val(2), Val(1))
        return @benchmarkable EXP._ind_window_loop!($b, $sys, $x0, $W1, $red, $st0, 1, $N, 1, true, $feed) samples = 5 evals = 1
    end
    W0 = deviation_frame(Float64, "pi/4", 2)
    st0 = EXP.Tangent._init_states(red, Float64, Val(2), Val(2))
    nws = EXP.Tangent._chain_workspace(red, Float64, Val(2), Val(2))
    ws = nws == 0 ? EXP.Tangent._NO_WS : zeros(UInt32, nws)
    return @benchmarkable EXP._ind_series_loop!($b, $sys, $x0, $W0, $red, $st0, $ws, 1, $N, 1, true) samples = 5 evals = 1
end

for (name, spec) in IND_SPECS
    SUITE["run indicator $name"] = ind_item_bench(spec)
end
