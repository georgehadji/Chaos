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
