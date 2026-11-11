# AC10 to AC12: portrait data (chaotic fractions, the printed orbits R, S, C), the fractions across the three energies, escape above E = 1/6.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).
# Package sizes: 24 x 24 grid at E = 1/8 (full 200 x 200), 14 x 14 for the energy scan (full 100 x 100); the paper's 500 x 500 is the preregistered size of the study.

hh_exit_centres() = (pi / 2, 7pi / 6, 11pi / 6)                    # Equation 14: sectors k = 0, 1, 2
function hh_sector(z)                                              # nearest exit by the polar angle of q (test-side, independent of the module)
    th = atan(z[2], z[1])
    d = [abs(mod(th - c + pi, 2pi) - pi) for c in hh_exit_centres()]
    return argmin(d) - 1
end

@testset "AC10 portrait reproduction at E = 1/8 (AntonopoulosManosSkokos2010 Sect. 3.1, O5)" begin
    _ = (energy_grid, energy_start)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); E = 1 / 8
    p = hh_portrait(m, E, hh_size(24, 200), 2000.0)
    @test p.n > 0
    @test abs(p.fraction - 0.55) <= 0.10                                                          # (a) about 55 percent, 10 percentage points (proposal)
    zR, sR = energy_start(flow, E, 0.1, 0.0); zS, sS = energy_start(flow, E, -0.0135, 0.0); zC, sC = energy_start(flow, E, -0.015, 0.0)
    @test sR === :ok && sS === :ok && sC === :ok
    # the printed p1 (spec Equation 9) reproduces H = 1/8 only to about 2e-7, so the rebuilt p1 differs from it by about 4e-7
    @test abs(zR[3] - 0.490578) <= 1.0e-6 && abs(zS[3] - 0.499816) <= 1.0e-6 && abs(zC[3] - 0.499773) <= 1.0e-6
    @test 0.5 <= hh_sali(m, flow, zR, 2000.0) <= 2                                                # (b) ordered orbit R, printed about 1.25
    @test hh_sali(m, flow, zC, 2000.0) < 1.0e-8                                                   # (c) chaotic orbit C
    @test minimum(hh_sali(m, flow, zC, T) for T in (1000.0, 2000.0, 4000.0, 10000.0)) <= 1.0e-12
    @info "AC10 reported: chaotic fraction, sticky orbit S (printed: turns chaotic after about 5000 steps)" fraction = p.fraction S2000 = hh_sali(m, flow, zS, 2000.0) S10000 = hh_sali(m, flow, zS, 10000.0)
end

@testset "AC11 chaotic fractions across the three energies and the integrator effect" begin
    _ = (energy_grid,)
    n = hh_size(14, 100); f8 = Float64[]; f4 = Float64[]
    for E in (1 / 12, 1 / 8, 1 / 6)
        push!(f8, hh_portrait(KahanLi(Float64, 8), E, n, 2000.0).fraction)
        push!(f4, hh_portrait(ForestRuth(Float64), E, n, 2000.0).fraction)
    end
    for f in (f8, f4)
        @test f[1] + 0.05 <= f[2] && f[2] + 0.05 <= f[3]                                          # f(1/12) < f(1/8) < f(1/6), gaps at least 5 percentage points
    end
    @test all(abs.(f8 .- f4) .<= 0.03)                                                           # order-8 against order-4 table at tau = 1/20
    @info "AC11 reported: fractions at E = 1/12, 1/8, 1/6" f8 f4
end

@testset "AC12a escape: declared at r >= R_esc, no return, sector by the polar angle" begin
    _ = (random_energy_start, escape_orbit)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; R = 10.0; et = hh_eps(Float64); n = round(Int, 200 / tau)
    nst = hh_size(100, 1000); nescaped = 0; nret = 0; nbadr = 0; nbads = 0; nbadt = 0
    for (k, E) in enumerate((0.2, 0.25))
        rng = StableRNG(hh_seed(12, k))
        for _ in 1:nst
            z0, st = random_energy_start(flow, E, (-1.0, 1.0), (-1.0, 1.0), rng)
            res = escape_orbit(m, flow, z0, tau, n; R_esc = R)
            res.status === :escaped || continue
            nescaped += 1
            hh_rad(res.z) >= R * (1 - 16 * et) || (nbadr += 1)
            res.sector == hh_sector(res.z) || (nbads += 1)
            0 < res.t <= n * tau || (nbadt += 1)
            z = res.z; rmin = Inf                                                                # continued to r = 2 R_esc or blow-up
            for _ in 1:4000
                z = advance(m, flow, z, tau)
                (all(isfinite, z) && hh_rad(z) <= 2R) || break
                rmin = min(rmin, hh_rad(z))
            end
            rmin >= R / 2 || (nret += 1)
        end
    end
    @test nescaped > 0                                                                           # some starts at E > 1/6 escape within T = 200
    @test nbadr == 0 && nbads == 0 && nbadt == 0 && nret == 0
    @test isequal(escape_orbit(m, flow, SVector(0.0, 0.1, 0.3, 0.1), tau, n), escape_orbit(m, flow, SVector(0.0, 0.1, 0.3, 0.1), tau, n; R_esc = 10.0))   # default R_esc = 10 (P5)
end

@testset "AC12b escape time against the reference" begin
    _ = (random_energy_start, escape_orbit)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; R = 10.0; n = round(Int, 200 / tau)
    nref = hh_size(8, 100); done = 0; worst = 0.0; nmis = 0
    for (k, E) in enumerate((0.2, 0.25))
        rng = StableRNG(hh_seed(12, 10 + k))
        while done < nref * k
            z0, _ = random_energy_start(flow, E, (-1.0, 1.0), (-1.0, 1.0), rng)
            res = escape_orbit(m, flow, z0, tau, n; R_esc = R)
            (res.status === :escaped && res.t < 100) || continue
            ref = hh_ref_escape(z0, R, 100.0)
            ref === nothing ? (nmis += 1) : (worst = max(worst, abs(res.t - ref)))
            done += 1
        end
    end
    @test nmis == 0 && worst <= 1.0e-6
end

@testset "AC12c exit label equivariance: reflection bitwise, rotation shifts the label" begin
    _ = (random_energy_start, escape_orbit)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; R = 10.0; n = round(Int, 200 / tau)
    nst = hh_size(100, 1000); nrefl = 0; nrot_mis = 0; nearly = 0; nrot_esc = 0
    setprecision(BigFloat, HH_BITS) do
        G = hh_rot_big(1)
        for (k, E) in enumerate((0.2, 0.25))
            rng = StableRNG(hh_seed(12, 20 + k))
            for _ in 1:nst
                z0, _ = random_energy_start(flow, E, (-1.0, 1.0), (-1.0, 1.0), rng)
                res = escape_orbit(m, flow, z0, tau, n; R_esc = R)
                rr = escape_orbit(m, flow, hh_reflect(z0), tau, n; R_esc = R)
                if res.status === :escaped
                    (rr.status === :escaped && rr.t == res.t && rr.sector == mod(-res.sector, 3)) || (nrefl += 1)   # labels 1 and 2 swapped, 0 fixed
                else
                    rr.status === res.status || (nrefl += 1)
                end
                if res.status === :escaped && res.t < 50
                    nearly += 1
                    rot = escape_orbit(m, flow, hh_apply(G, z0), tau, n; R_esc = R)
                    if rot.status === :escaped
                        nrot_esc += 1
                        rot.sector == mod(res.sector + 1, 3) || (nrot_mis += 1)
                    end
                end
            end
        end
    end
    @test nrefl == 0
    @test nearly > 0 && nrot_mis == 0 && nrot_esc >= 0.9 * nearly                                # the rotation is rounded: almost every early escape escapes again, with the shifted label
end

@testset "AC12d non-finite and overflowing starts end with a status, never an exception" begin
    _ = (escape_orbit, section_map, random_energy_start)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05
    @test escape_orbit(m, flow, SVector(NaN, 0.0, 0.0, 0.0), tau, 100).status === :nonfinite
    @test escape_orbit(m, flow, SVector(0.0, Inf, 0.0, 0.0), tau, 100).status === :nonfinite
    @test escape_orbit(m, flow, SVector(1.0e100, 1.0, 0.0, 0.0), tau, 100; R_esc = Inf).status === :nonfinite   # the step overflows
    @test section_map(m, flow, SVector(NaN, 0.0, 0.0, 0.0), tau, 100).status === :nonfinite
    @test section_map(m, flow, SVector(-1.0e100, 1.0, 0.0, 0.0), tau, 100).status === :nonfinite
    z, st = random_energy_start(flow, NaN, (-1.0, 1.0), (-1.0, 1.0), StableRNG(hh_seed(12, 30)))
    @test z === nothing && st !== :ok
end
