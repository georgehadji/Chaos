# Helpers for the experiments suite (prefix ex_). Definitions only: nothing here calls Chaoslab.Experiments at load time.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, TOML

# EX.name resolves Chaoslab.Experiments.name when called, so a missing module is an UndefVarError inside the testset.
struct ExLazy end
Base.getproperty(::ExLazy, s::Symbol) = getfield(Chaoslab.Experiments, s)
const EX = ExLazy()

const EX_SEED = 20261012                    # _config/numerics-defaults.toml, seed
const EX_BASE = EX_SEED + 27000             # this slug uses 20261012 + 27000 .. 27999; testset AC<n> uses 50 n + sub (sub < 50)
ex_seed(ac::Integer, sub::Integer = 0) = EX_BASE + 50 * ac + sub

const EX_RHO = 1e-10                        # fixture residual tolerance (the study value is not used in any test)
const EX_DUP = 1e-8
const EX_EPS = ldexp(1.0, -52)
const EX_EPS_BF = ldexp(1.0, -255)          # spacing at 1 of a 256-bit BigFloat
const EX_BITS = 1024
const EX_TWO_PI = 2 * pi

# ---- independent SHA-256 (platform command-line tool, spec O3); the module's value is compared against it
function ex_sha256(bytes::AbstractVector{UInt8})
    mktemp() do path, io
        write(io, bytes)
        close(io)
        for cmd in (`sha256sum $path`, `certutil -hashfile $path SHA256`, `shasum -a 256 $path`)
            out = try
                read(cmd, String)
            catch
                continue
            end
            m = match(r"\b[0-9a-fA-F]{64}\b", out)
            m === nothing || return lowercase(m.match)
        end
        error("no command-line SHA-256 tool found")
    end
end
ex_sha256(s::AbstractString) = ex_sha256(Vector{UInt8}(codeunits(s)))

# E18: CRLF to LF, then the header value cell replaced by SELF (the spaces around the cell text stay)
ex_prepare(text::AbstractString) =
    replace(replace(text, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
ex_header_hash(text::AbstractString) = match(r"(?m)^\| SHA-256 of this file \(current version\) \| *([0-9a-f]{64})", text).captures[1]

# ---- evaluation of closed-form input strings ("pi/2", "(3 + sqrt(5))/2", "0.5") at the current BigFloat precision
function ex_eval(x)
    x isa Expr && x.head === :call && return getfield(Base, x.args[1])(map(ex_eval, x.args[2:end])...)
    x === :pi && return BigFloat(pi)
    x isa Number && return BigFloat(x)
    x isa Expr && x.head === :block && return ex_eval(last(filter(a -> !(a isa LineNumberNode), x.args)))
    error("ex_eval: unsupported $x")
end
function ex_eval(s::AbstractString)
    e = Meta.parse(replace(s, r"(\d)pi" => s"\1*pi"))
    e isa Number && return parse(BigFloat, s)           # decimals: the exact decimal, not a Float64 first
    return ex_eval(e)
end
ex_f64(s::AbstractString) = setprecision(() -> Float64(ex_eval(s)), BigFloat, 256)

# ---- fixture preregistration: header with the SELF hash rule, an amendment log, one ```toml panel block
function ex_block_default()
    return Dict{String, Any}(
        "K" => ["0", "1.5", "5"], "n_c" => 3, "grid_side" => 8, "N_scr" => 2000, "alpha" => "pi/4",
        "theta_T" => 0.01, "g_T" => 5.0, "r_T" => 0.01, "theta_C" => 0.05, "r_C" => 0.05,
        "theta_S" => 0.02, "w" => 100, "f_S" => 0.2, "L_S" => 500, "log10_r" => [-4.0, -2.0],
        "q_max" => 2, "rho" => EX_RHO, "delta_R" => 1e-3, "newton_max" => 50, "delta_dup" => EX_DUP,
        "n_ctrl_max" => 4, "mu" => "(3 + sqrt(5))/2", "s_D" => 2, "cap_D" => 100, "s_L" => 1, "cap_L" => 50)
end

function ex_block(over::AbstractDict = Dict{String, Any}(); drop = String[])
    d = merge(ex_block_default(), Dict{String, Any}(over))
    foreach(k -> delete!(d, k), drop)
    return d
end

function ex_prereg_text(block::AbstractDict; amendments::Int = 2, hash_rows::Int = 1)
    io = IOBuffer()
    println(io, "# Preregistration: ex-fixture\n\nFixture for the experiments suite. Not a study.\n\n## Header\n")
    println(io, "| Field | Value |\n|-------|-------|\n| Study slug | ex-fixture |")
    for _ in 1:hash_rows
        println(io, "| SHA-256 of this file (current version) | ", "0"^64, " (rule: fixture) |")
    end
    println(io, "| Amendment log | ", join(("A$i (2026-01-0$i, fixture)" for i in 1:amendments), "; "), " |\n")
    println(io, "## Panel\n\nProse that is not read by the module.\n\n## Amendment A1\n\n```toml")
    TOML.print(io, Dict("panel" => block))
    println(io, "```\n\nEnd of file.")
    return String(take!(io))
end

# Writes the fixture and puts the true hash into the header cell (the cell is replaced by SELF before hashing, so the
# placeholder digits do not matter).
function ex_write_prereg(path::AbstractString, block::AbstractDict = ex_block(); kw...)
    text = ex_prereg_text(block; kw...)
    h = ex_sha256(ex_prepare(text))
    text = replace(text, "0"^64 => h; count = 1)
    write(path, text)
    return path
end

const EX_DIR = mktempdir()
ex_path(name) = joinpath(EX_DIR, name)

# ---- 1024-bit reference algebra for the standard map (call inside setprecision(BigFloat, EX_BITS))
ex_w(x) = zero(BigFloat) + x                                   # widen to the current precision, exactly
ex_T(K, z) = (p = z[2] + K * sin(z[1]); (z[1] + p, p))
function ex_Tq(K, z, q)
    for _ in 1:q
        z = ex_T(K, z)
    end
    return z
end
ex_DT(K, x) = (c = K * cos(x); [1+c 1; c 1])
function ex_J(K, z, q)
    J = Matrix{typeof(z[1])}(I, 2, 2)
    for _ in 1:q
        J = ex_DT(K, z[1]) * J
        z = ex_T(K, z)
    end
    return J
end
ex_delta(s) = (tp = 2 * oftype(s, pi); abs(s - tp * round(s / tp)))                  # E8
ex_tdist(u, v) = sqrt(ex_delta(u[1] - v[1])^2 + ex_delta(u[2] - v[2])^2)
ex_wrap(s) = (tp = 2 * oftype(s, pi); s - tp * floor(s / tp))
ex_residual(K, z, q) = (zq = ex_Tq(K, z, q); max(ex_delta(zq[1] - z[1]), ex_delta(zq[2] - z[2])))
ex_orbit_trace(K, z, q) = (J = ex_J(K, z, q); J[1, 1] + J[2, 2])
ex_exponent(t, q) = abs(t) <= 2 ? zero(t) : log((abs(t) + sqrt(t^2 - 4)) / 2) / q    # E4 (parabolic as 0 handled by callers)

# ---- cached periodic sets (the search is the slow part of the suite)
const EX_PSETS = Dict{Any, Any}()
ex_pset(K::Real; q_max = 2, N_g = 32) = get!(EX_PSETS, (Float64(K), q_max, N_g)) do
    EX.periodic_set(Float64(K); q_max = q_max, N_g = N_g, rho = EX_RHO, newton_max = 50, delta_dup = EX_DUP)
end

# an orbit of the set with a point within tol of z (BigFloat or Float64 pair), period q; nothing when absent
function ex_orbit_with(ps, q, z; tol = 1e4 * EX_RHO^2)
    for o in ps.orbits
        o.q == q || continue
        any(pt -> ex_tdist(ex_w.(pt), z) <= tol, o.points) && return o
    end
    return nothing
end

# ---- cached panel builds: tag -> (parsed file, path)
const EX_BUILDS = Dict{String, Any}()
# The manifest records the preregistration path, so equal blocks share one file: byte-identical panels need an identical path.
function ex_prereg_for(over::AbstractDict = Dict{String, Any}())
    block = ex_block(over)
    path = ex_path("prereg-" * ex_sha256(ex_prereg_text(block))[1:12] * ".md")
    isfile(path) || ex_write_prereg(path, block)
    return path
end

function ex_build(tag::AbstractString, over::AbstractDict = Dict{String, Any}(); kw...)
    get!(EX_BUILDS, tag) do
        pre = ex_prereg_for(over)
        out = ex_path("$tag-panel.toml")
        EX.build_panel(pre; out = out, kw...)
        (panel = TOML.parsefile(out), path = out, prereg = pre)
    end
end

# ---- E15 and E14 reference implementations written from the spec
ex_cell_seed(run_seed, i, class) = reinterpret(UInt64, reverse(hex2bytes(ex_sha256("chaoslab-panel-1|$(run_seed)|$(i)|$(class)"))[1:8]))[1]
function ex_ref_draw(rng, cands, n_c)
    c = collect(Int, cands)
    M = length(c)
    out = Int[]
    for t in 1:min(n_c, M)
        j = t + floor(Int, rand(rng) * (M - t + 1))
        c[t], c[j] = c[j], c[t]
        push!(out, c[t])
    end
    return out
end

# one-sample Kolmogorov-Smirnov statistic against the uniform law on [0, 1]
function ex_ks_uniform(u)
    s = sort(u)
    n = length(s)
    return maximum(max(abs(s[i] - (i - 1) / n), abs(s[i] - i / n)) for i in 1:n)
end

# every key of a nested TOML value, recursively
function ex_all_keys(x, acc = String[])
    if x isa AbstractDict
        for (k, v) in x
            push!(acc, String(k))
            ex_all_keys(v, acc)
        end
    elseif x isa AbstractVector
        foreach(v -> ex_all_keys(v, acc), x)
    end
    return acc
end

# ---- device plumbing (copy of the tangent-core pattern): the run-wide device is set and restored explicitly
const ex_F = Chaoslab.Foundation
const ex_Tg = Chaoslab.Tangent
struct ExFakeDevice <: ex_F.AbstractDevice
    loadable::Bool
    functional::Bool
end
ex_F.probe(d::ExFakeDevice, ::Type) = (loadable = d.loadable, functional = d.functional, device = "FakeGPU", compiles = true, gpu_faster = true)
ex_Tg.device_backend(::ExFakeDevice) = ex_Tg.KernelAbstractions.CPU()
function ex_with_device(f, dev)
    old = ex_F.DEVICE[]
    ex_F.DEVICE[] = dev
    try
        return f()
    finally
        ex_F.DEVICE[] = old
    end
end

# true when the call throws an ArgumentError whose message contains one of the keys
function ex_rejects(f, keys)
    try
        f()
    catch e
        return e isa ArgumentError && any(k -> occursin(k, e.msg), keys)
    end
    return false
end

# lifted starts of a parsed panel, by orbit
ex_orbits(panel) = panel["orbit"]
ex_start(o) = (Float64(o["x0"]), Float64(o["p0"]))

# spec S3 as amended (D35): B_S(n) = 8 eps (n + G_first) + 16 |K| n^2 eps, eps = eps_53 spacing; G_first of the first handoff window (tau = 1)
ex_BS(n, K, G_first) = 8 * EX_EPS * (n + G_first) + 16 * abs(K) * n^2 * EX_EPS
ex_Gfirst(K, x) = (c = K * cos(x); (1 + c)^2 + 1 + c^2 + 1)
