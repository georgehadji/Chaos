# AC17 to AC23 and the API testsets: element-type agreement, backends (CPU kernel against a plain loop; GPU is stage 04), the Taylor reference R1 (extension, skipped
# cleanly when the weak dependency is not loaded; the verify environment covers it), health, layering. AC20 is reserved (R2 deferred) and has no test.

# step counts and taus of the structure flows (n_H tau = 10 time units)
const IC_NH = Dict(:f4 => (0.05, 200), :f2 => (0.01, 1000))

@testset "AC17 element-type agreement: binary64 and double-double against 256 bits" begin
    _ = (advance_n, Gauss, effective_eps)
    @test ic_eps(Double64) <= 2.0^(-100)                                                  # effective double-double epsilon (precision spec); otherwise report to the owner
    rng = StableRNG(ic_seed(17, 1))
    for M in IC_ALLX, kind in (:f4, :f2)
        tau, n = IC_NH[kind]
        for _ in 1:ic_size(2, 200)                                                        # D38: 2 orbits per method and flow (spec 200); the 256-bit orbit is shared by both types
            z0 = ic_draw(rng, kind); Mn = ic_mmax(z0)
            zb = setprecision(BigFloat, IC_BITS) do
                Vector(advance_n(M.mk(BigFloat), ic_flow(BigFloat, kind), ic_to(BigFloat, z0), BigFloat(tau), n))
            end
            for T in (Float64, Double64)
                zt = advance_n(M.mk(T), ic_flow(T, kind), ic_to(T, z0), T(tau), n)
                @test ic_err(zt, zb) <= ic_B(T, Mn, n, tau)
            end
        end
    end
    for T in (Float64, Double64), (s, mk) in IC_GAUSS[1:2]                                 # implicit: plus the solver tolerance per step
        for _ in 1:ic_size(2, 200)
            z0 = ic_draw(rng, :f4); n = 100
            zt = advance_n(mk(T), ic_flow(T, :f4), ic_to(T, z0), T(0.1), n)
            zb = setprecision(BigFloat, IC_BITS) do
                Vector(advance_n(mk(BigFloat), ic_flow(BigFloat, :f4), ic_to(BigFloat, z0), BigFloat(0.1), n))
            end
            @test ic_err(zt, zb) <= ic_B(T, ic_mmax(z0), n, 0.1) + 8 * ic_eps(T) * ic_mmax(z0) * n
        end
    end
end

@testset "AC18 CPU kernel against a plain loop, orbit results independent of the chunking" begin
    _ = (advance_ensemble, advance_n, advance_checked)
    rng = StableRNG(ic_seed(18, 1))
    cases = ((ICM("S2", T -> Leapfrog(T, :DKD), 2, true, Inf)), IC_S[7], IC_S[14], IC_S[19], IC_P[3])   # S2, S5, S8, S12-6 (a table of S12), SABA3
    norb = ic_size(500, 100000)
    for M in cases, kind in (:f4, :f2)
        tau, n = IC_NH[kind]; n = min(n, 200)
        x0 = [ic_to(Float64, kind === :f2 ? ic_kepler_state(0.5, 2 * Float64(pi) * rand(rng)) : ic_draw(rng, kind)) for _ in 1:norb]
        m = M.mk(Float64); flow = ic_flow(Float64, kind)
        states, ok = advance_ensemble(m, flow, x0, tau, n; backend = :cpu)
        @test all(ok) && length(states) == norb
        loop = [advance_n(m, flow, x, tau, n) for x in x0]
        @test maximum(maximum(abs, collect(states[j]) .- collect(loop[j])) for j in 1:norb) <= ic_B(Float64, 4.0, n, tau)
        # result of orbit j does not depend on how the orbits are distributed (chunk sizes 7 and 1000; the thread count is the run's)
        s7, _ = advance_ensemble(m, flow, x0, tau, n; backend = :cpu, chunk = 7)
        s1000, _ = advance_ensemble(m, flow, x0, tau, n; backend = :cpu, chunk = 1000)
        @test all(s7[j] == s1000[j] == states[j] for j in 1:norb)                          # bitwise
    end
    # implicit: per-orbit failure statuses; a benign orbit and an orbit started at r = 0.02 (diverging iteration) in one ensemble
    benign = ic_to(Float64, ic_kepler_state(0.1, 0.7)); bad = SVector(0.02, 0.0, 0.0, 7.0)
    m = Gauss(Float64, 2); flow = ic_flow(Float64, :f2)
    states, ok = advance_ensemble(m, flow, [benign, bad, benign], 0.05, 10; backend = :cpu)
    @test ok == [true, false, true]
    @test states[1] == states[3] == advance_n(m, flow, benign, 0.05, 10)
    @test !any(isfinite, states[2])                                                         # never a silent value
end

@testset "AC19 R1: Taylor reference declared in Integrators, methods from the extension" begin
    _ = (taylor_reference,)
    proj = read(joinpath(@__DIR__, "..", "..", "Project.toml"), String)
    @test occursin(r"\[weakdeps\][^\[]*TaylorIntegration"s, proj) && occursin(r"\[extensions\][^\[]*ChaoslabTaylorExt\s*=\s*\"TaylorIntegration\""s, proj)
    @test isfile(joinpath(@__DIR__, "..", "..", "ext", "ChaoslabTaylorExt.jl"))
    ext = Base.get_extension(Chaoslab, :ChaoslabTaylorExt)
    if ext === nothing
        @test isempty(methods(taylor_reference))                                            # the library never depends on the Taylor package
        @test_skip ext !== nothing                                                          # criteria 19 and 21 run in the verify environment, where TaylorIntegration is loaded
    else
        setprecision(BigFloat, IC_BITS) do
            z0 = IC_START[:f1_3]
            zT = taylor_reference(ic_flow(BigFloat, :f1_3), ic_to(BigFloat, z0), BigFloat(100))
            @test ic_err(zT, ic_f1_exact(:f1_3, z0, BigFloat(100))) <= 2.0^(-180)
            k0 = ic_to(BigFloat, IC_START[:f2]); fl = ic_flow(BigFloat, :f2)
            kT = taylor_reference(fl, k0, BigFloat(100))
            @test abs(hamiltonian(fl, kT) - hamiltonian(fl, k0)) <= 2.0^(-180)
            @test abs(angular_momentum(fl, kT) - angular_momentum(fl, k0)) <= 2.0^(-180)
            # reference tangent by central finite difference of orbits (step 2^-60) against the exact matrix exponential on F1
            f1 = ic_flow(BigFloat, :f1_3); dl = BigFloat(2)^(-60); worst = 0.0
            for j in 1:6
                e = zeros(BigFloat, 6); e[j] = dl
                zp = taylor_reference(f1, ic_to(BigFloat, z0) .+ SVector{6}(e), BigFloat(100)); zm = taylor_reference(f1, ic_to(BigFloat, z0) .- SVector{6}(e), BigFloat(100))
                ex = ic_f1_exact(:f1_3, [i == j ? 1.0 : 0.0 for i in 1:6], BigFloat(100))
                worst = max(worst, Float64(maximum(abs, Vector(zp .- zm) ./ (2 * dl) .- ex)))
            end
            @test worst <= 2.0^(-100)
        end
    end
end

@testset "AC21 methods against R1 (two oracles agree on the error), extension only" begin
    _ = (taylor_reference, advance_n)
    ext = Base.get_extension(Chaoslab, :ChaoslabTaylorExt)
    if ext === nothing
        @test_skip ext !== nothing
    else
        for M in IC_S
            (M.p >= 2 && M.sym !== false && M.p <= 4) || continue
            z0 = IC_START[:f4]; ref = ic_ref(:f4, z0, 8); tau = 1 / 16; n = 128
            zt = advance_n(M.mk(Float64), ic_flow(Float64, :f4), ic_to(Float64, z0), tau, n)
            tay = setprecision(BigFloat, IC_BITS) do
                Vector(taylor_reference(ic_flow(BigFloat, :f4), ic_to(BigFloat, z0), BigFloat(8)))
            end
            @test abs(ic_err(zt, tay) - ic_err(zt, ref)) <= ic_B(Float64, 1.0, n, tau)
        end
    end
end

@testset "AC22 failures are visible" begin
    _ = (SplittingMethod, advance, advance_checked, advance_ensemble, Gauss, TaoExtended)
    ok = ([0.5, 0.5], [1.0, 0.0])
    @test_throws ArgumentError SplittingMethod([0.5, 0.5 + 1e-6], [1.0, 0.0]; order = 2, palindromic = false, source = "Yoshida1990", name = "bad sum")
    @test_throws ArgumentError SplittingMethod([0.5, 0.5], [1.0 + 1e-6, 0.0]; order = 2, palindromic = false, source = "Yoshida1990", name = "bad sum d")
    @test_throws ArgumentError SplittingMethod([0.25, 0.25, 0.3, 0.2], [0.5, 0.25, 0.25, 0.0]; order = 2, palindromic = true, source = "Yoshida1990", name = "c2 != c3")
    @test_throws ArgumentError SplittingMethod(Float64[], Float64[]; order = 1, palindromic = false, source = "Yoshida1990", name = "no stages")
    good = SplittingMethod(ok...; order = 2, palindromic = true, source = "Yoshida1990", name = "leapfrog by hand")
    z = SVector(1.0, 0.0); flow = HarmonicOscillator([1.0])
    @test advance(good, flow, z, 0.5) == advance(Leapfrog(Float64, :DKD), flow, z, 0.5)
    # non-finite input: non-finite output without an error, for every explicit method on the four separable flows
    for M in IC_ALLX, kind in (:f1_3, :f4, :f2, :t2)
        m = M.mk(Float64); fl = ic_flow(Float64, kind)
        z0 = ic_to(Float64, IC_START[kind]); zi = SVector{length(z0), Float64}(ntuple(i -> i == 1 ? Inf : z0[i], length(z0)))
        out = try
            advance(m, fl, zi, 0.05)
        catch
            nothing
        end
        @test out !== nothing && !all(isfinite, collect(out))
    end
    # implicit and Tao: no throw, failure visible
    for (s, mk) in IC_GAUSS
        zn, st = advance_checked(mk(Float64), ic_flow(Float64, :f4), SVector(NaN, 0.1), 0.1)
        @test !st.ok && !any(isfinite, zn)
    end
    zt = try
        advance(TaoExtended(Float64, 10.0; order = 2), TaoNonseparable(), SVector(Inf, 0.1, Inf, 0.1), 0.01)
    catch
        nothing
    end
    @test zt !== nothing && !all(isfinite, collect(zt))
    # an ensemble with a non-finite orbit: the others are unaffected, nothing throws inside the kernel
    xs = [ic_to(Float64, [1.0, 0.3]), SVector(NaN, 0.0), ic_to(Float64, [0.5, -0.2])]
    res = try
        advance_ensemble(Leapfrog(Float64, :DKD), ic_flow(Float64, :f4), xs, 0.05, 20; backend = :cpu)
    catch
        nothing
    end
    @test res !== nothing && res[1][1] == advance_n(Leapfrog(Float64, :DKD), ic_flow(Float64, :f4), xs[1], 0.05, 20) && !any(isfinite, res[1][2])
end

# upward imports: every module named on a using/import line must be Foundation, Systems or Integrators itself (spec criterion 23)
function ic_upward(lines)
    bad = String[]
    for line in lines
        occursin(r"^\s*(?:using|import)\s", line) || continue
        for m in eachmatch(r"(?:\.\.|Chaoslab\.)(\w+)", line)
            m.captures[1] in ("Foundation", "Systems", "Integrators") || push!(bad, m.captures[1])
        end
    end
    return bad
end

@testset "AC23 layering of the Integrators sources" begin
    _ = (Chaoslab,)
    root = joinpath(@__DIR__, "..", "..", "src")
    dir = joinpath(root, "Integrators")
    @test isdir(dir)
    files = isdir(dir) ? [joinpath(d, f) for (d, _, fs) in walkdir(dir) for f in fs if endswith(f, ".jl")] : String[]
    @test !isempty(files)
    for f in files
        @test isempty(ic_upward(readlines(f)))
    end
    # Systems must not import Integrators (the flow side is new in Systems but stays below)
    for (d, _, fs) in walkdir(joinpath(root, "Systems")), f in fs
        endswith(f, ".jl") && @test !("Integrators" in ic_upward(readlines(joinpath(d, f))))
    end
    # position in the package file: Systems, then Integrators, then Tangent
    top = read(joinpath(root, "Chaoslab.jl"), String)
    a = findfirst("Systems/Systems.jl", top); b = findfirst("Integrators/Integrators.jl", top); c = findfirst("Tangent/Tangent.jl", top)
    @test a !== nothing && b !== nothing && c !== nothing && first(a) < first(b) < first(c)
    # the checker itself catches an upward import (negative control)
    @test ic_upward(["using ..Foundation, ..Systems", "import ..Tangent: propagate", "using ..Lyapunov"]) == ["Tangent", "Lyapunov"]
end

@testset "API: hand values of the first-order and second-order maps, IntegratorMap, Tangent consumption" begin
    _ = (SymplecticEuler, Leapfrog, IntegratorMap, advance, forward, inverse, dimension, nominal_order, method_name, propagate)
    flow = HarmonicOscillator([1.0]); z = SVector(1.0, 0.0)
    # H = (p^2 + q^2)/2, tau = 0.5, all dyadic: hand values
    @test advance(Leapfrog(Float64, :DKD), flow, z, 0.5) == SVector(0.875, -0.5)            # D(1/4) K(1/2) D(1/4)
    @test advance(Leapfrog(Float64, :KDK), flow, z, 0.5) == SVector(0.875, -0.46875)        # K(1/4) D(1/2) K(1/4)
    @test advance(SymplecticEuler(Float64, :DK), flow, z, 0.5) == SVector(1.0, -0.5)        # drift then kick
    @test advance(SymplecticEuler(Float64, :KD), flow, z, 0.5) == SVector(0.75, -0.5)       # kick then drift
    @test nominal_order(Leapfrog(Float64)) == 2 && nominal_order(SymplecticEuler(Float64)) == 1
    @test !isempty(method_name(Leapfrog(Float64)))
    # the one-step map as a Systems map
    m = Leapfrog(Float64, :DKD); fl = ic_flow(Float64, :f4); im = IntegratorMap(m, fl, 0.05)
    @test im isa AbstractMap && dimension(im) == 2
    rng = StableRNG(ic_seed(24, 1))
    for _ in 1:20
        x = ic_to(Float64, ic_draw(rng, :f4))
        @test forward(im, x) == advance(m, fl, x, 0.05)
        @test maximum(abs, collect(inverse(im, forward(im, x))) .- collect(x)) <= 16 * ic_eps(Float64) * ic_mmax(x)
        @test maximum(abs, Matrix(inverse_tangent(im, forward(im, x))) * Matrix(tangent(im, x)) - I) <= 16 * ic_eps(Float64) * (1 + maximum(abs, Matrix(tangent(im, x)))^2)
    end
    # the Tangent layer consumes it like any map: 8 orbits, 20 steps, bitwise the serial iteration of forward
    x0 = [ic_to(Float64, ic_draw(rng, :f4)) for _ in 1:8]
    res = propagate(im, x0, 20; k = 2, seed = ic_seed(24, 2))
    for j in 1:8
        x = x0[j]
        for _ in 1:20
            x = forward(im, x)
        end
        @test res.x[j] == x
    end
    # propagate with the symmetric default: a negative step is the inverse (I8)
    x = x0[1]
    @test maximum(abs, collect(advance(m, fl, advance(m, fl, x, 0.05), -0.05)) .- collect(x)) <= 16 * ic_eps(Float64) * ic_mmax(x)
end
