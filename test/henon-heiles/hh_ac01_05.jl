# AC1 to AC5: the system (critical points, force and Hessian, D3 and time-reversal symmetry), energy-surface starts and the bound below the escape energy.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError). Testsets whose names end in "(existing code)" use only
# names that exist today (HenonHeilesFlow, force, force_jacobian, hamiltonian, advance_n, propagate) and are expected to pass already; they are listed in the test plan.

@testset "AC1a critical points and levels from the force and Hamiltonian of the flow (existing code)" begin
    flow = HenonHeilesFlow()
    setprecision(BigFloat, HH_BITS) do
        s3 = sqrt(BigFloat(3)); half = BigFloat(1) / 2; tol = BigFloat(2)^-200
        pts = ((BigFloat(0), BigFloat(0)), (BigFloat(0), BigFloat(1)), (s3 / 2, -half), (-s3 / 2, -half))   # O, S0, S1, S2 (Equation 4)
        vals = (BigFloat(0), BigFloat(1) / 6, BigFloat(1) / 6, BigFloat(1) / 6)
        eigs = ((BigFloat(1), BigFloat(1)), (BigFloat(3), BigFloat(-1)), (BigFloat(3), BigFloat(-1)), (BigFloat(3), BigFloat(-1)))
        for k in 1:4
            q = SVector{2, BigFloat}(pts[k]...)
            @test maximum(abs, force(flow, q)) <= tol
            @test abs(hamiltonian(flow, SVector{4, BigFloat}(q[1], q[2], 0, 0)) - vals[k]) <= tol
            Hs = force_jacobian(flow, q); tr = Hs[1, 1] + Hs[2, 2]; dt = Hs[1, 1] * Hs[2, 2] - Hs[1, 2] * Hs[2, 1]
            disc = sqrt(tr^2 - 4dt)
            @test abs((tr + disc) / 2 - eigs[k][1]) <= tol && abs((tr - disc) / 2 - eigs[k][2]) <= tol
        end
        # the points (+-sqrt(3)/2, +1/2) printed in benchmarks.md before the D116 correction are not critical: gradient norm at least 1 (it is 2)
        for sgn in (1, -1)
            @test norm(force(flow, SVector{2, BigFloat}(sgn * s3 / 2, half))) >= 1
        end
    end
end

@testset "AC1b critical_points and escape_energy" begin
    _ = (critical_points, escape_energy)
    flow = HenonHeilesFlow()
    setprecision(BigFloat, HH_BITS) do
        s3 = sqrt(BigFloat(3)); half = BigFloat(1) / 2; tol = BigFloat(2)^-200
        exact = ((BigFloat(0), BigFloat(0)), (BigFloat(0), BigFloat(1)), (s3 / 2, -half), (-s3 / 2, -half))
        cp = critical_points(flow, BigFloat)
        @test length(cp) == 4 && all(c -> c isa SVector{2, BigFloat}, cp)
        for k in 1:4
            @test maximum(abs, cp[k] .- BigFloat.(collect(exact[k]))) <= tol
            @test maximum(abs, force(flow, cp[k])) <= tol                               # the module's own points are critical
        end
        @test abs(escape_energy(flow, BigFloat) - BigFloat(1) / 6) <= tol
        for T in (Float64, Double64)
            c = critical_points(flow, T)
            @test all(x -> x isa SVector{2, T}, c) && escape_energy(flow, T) isa T
            for k in 1:4, i in 1:2
                @test abs(hh_big(c[k][i]) - exact[k][i]) <= hh_eps(T)                   # correctly rounded data, error at most eps
            end
            @test abs(hh_big(escape_energy(flow, T)) - BigFloat(1) / 6) <= hh_eps(T)
        end
    end
end

@testset "AC2 force and Hessian against ForwardDiff and complex-step (existing code)" begin
    flow = HenonHeilesFlow()
    for (T, seed) in ((Float64, 1), (Double64, 2))
        rng = StableRNG(hh_seed(2, seed)); et = hh_eps(T)
        worst = 0.0
        for _ in 1:1000
            x = 3 .* rand(rng, 2) .- 1.5
            q = SVector{2, T}(T(x[1]), T(x[2]))
            gad = ForwardDiff.gradient(w -> hh_V(w[1], w[2]), T[q[1], q[2]])
            had = ForwardDiff.hessian(w -> hh_V(w[1], w[2]), T[q[1], q[2]])
            g = force(flow, q); Hs = force_jacobian(flow, q)
            for i in 1:2
                worst = max(worst, Float64(abs(g[i] - gad[i])) / (et * max(1, abs(Float64(gad[i])))))
                for j in 1:2
                    worst = max(worst, Float64(abs(Hs[i, j] - had[i, j])) / (et * max(1, abs(Float64(had[i, j])))))
                end
            end
            if T === Float64                                                               # complex step: gradient from V, Hessian column j from the force
                h = 1.0e-20
                for j in 1:2
                    e = SVector(j == 1 ? 1.0 : 0.0, j == 2 ? 1.0 : 0.0)
                    gcs = imag(hh_V(q[1] + im * h * e[1], q[2] + im * h * e[2])) / h
                    fcs = imag.(force(flow, SVector{2, ComplexF64}(q[1] + im * h * e[1], q[2] + im * h * e[2]))) ./ h
                    worst = max(worst, abs(g[j] - gcs) / (et * max(1, abs(gcs))))
                    for i in 1:2
                        worst = max(worst, abs(Hs[i, j] - fcs[i]) / (et * max(1, abs(fcs[i]))))
                    end
                end
            end
        end
        @test worst <= 8
    end
end

@testset "AC3a D3 elements and the invariance of V" begin
    _ = (d3_elements,)
    flow = HenonHeilesFlow()
    setprecision(BigFloat, HH_BITS) do
        F = BigFloat[-1 0; 0 1]
        expected = (BigFloat[1 0; 0 1], hh_rot_big(1), hh_rot_big(2), F, hh_rot_big(1) * F, hh_rot_big(2) * F)   # identity, rotations by 2 pi/3 and 4 pi/3, reflection q1 -> -q1, and its rotations
        for T in (Float64, Double64, BigFloat)
            G = d3_elements(T)
            @test length(G) == 6 && all(g -> g isa SMatrix{2, 2, T, 4}, G)
            for k in 1:6, i in 1:2, j in 1:2
                @test abs(hh_big(G[k][i, j]) - expected[k][i, j]) <= hh_eps(T)
            end
        end
        # (a) |V(G q) - V(q)| <= 8 eps max(1, |V|) at 1000 seeded points, 256 bits (the matrices are the module's, the points exact Float64 data)
        rng = StableRNG(hh_seed(3, 1)); Gb = d3_elements(BigFloat); et = hh_eps(BigFloat)
        for _ in 1:1000
            x = 3 .* rand(rng, 2) .- 1.5; q = SVector{2, BigFloat}(x[1], x[2])
            V = hamiltonian(flow, SVector{4, BigFloat}(q[1], q[2], 0, 0))
            for g in Gb
                gq = g * q
                @test abs(hamiltonian(flow, SVector{4, BigFloat}(gq[1], gq[2], 0, 0)) - V) <= 8 * et * max(1, abs(V))
            end
        end
        # the same in binary64: the rounding of the matrix and of the product enters as 4 eps |q|_1 |grad V|_1 (first order, derived)
        rng = StableRNG(hh_seed(3, 2)); G64 = d3_elements(Float64); e64 = hh_eps(Float64); bad = 0
        for _ in 1:1000
            x = 3 .* rand(rng, 2) .- 1.5; q = SVector(x[1], x[2])
            V = hamiltonian(flow, SVector(q[1], q[2], 0.0, 0.0)); g1 = force(flow, q)
            for g in G64
                gq = g * q
                abs(hamiltonian(flow, SVector(gq[1], gq[2], 0.0, 0.0)) - V) <= 8 * e64 * max(1, abs(V)) + 4 * e64 * sum(abs, q) * sum(abs, g1) || (bad += 1)
            end
        end
        @test bad == 0
    end
end

@testset "AC3b one step commutes with the reflection q1 -> -q1, bitwise (existing code)" begin
    flow = HenonHeilesFlow(); rng = StableRNG(hh_seed(3, 3)); m = KahanLi(Float64, 8)
    bad = 0
    for _ in 1:1000
        z = 2 .* rand(rng, 4) .- 1; z = SVector{4, Float64}(z)
        advance(m, flow, hh_reflect(z), 0.05) == hh_reflect(advance(m, flow, z, 0.05)) || (bad += 1)
    end
    @test bad == 0
end

@testset "AC3c rotation by 2 pi/3 commutes with the flow up to the amplified round-off" begin
    _ = (random_energy_start, section_interval)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; n = 1000
    rng = StableRNG(hh_seed(3, 4)); E = 1 / 8; a, b = section_interval(flow, E); pm = sqrt(2E)
    setprecision(BigFloat, HH_BITS) do
        G = hh_rot_big(1)
        for _ in 1:20
            z0, st = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
            @test st === :ok
            zn = advance_n(m, flow, z0, tau, n)
            zr = advance_n(m, flow, hh_apply(G, z0), tau, n)
            err = maximum(abs, zr .- hh_apply(G, zn))
            lam = max(hh_lambda1(m, flow, z0, tau, n * tau), 0.0)
            @test err <= 16 * hh_B(Float64, max(1.0, maximum(abs, zn), maximum(abs, z0)), n, tau) * exp(lam * n * tau)
        end
    end
end

@testset "AC3d time reversal: n steps, p -> -p, n steps return to the start (existing code)" begin
    flow = HenonHeilesFlow(); rng = StableRNG(hh_seed(3, 5)); tau = 0.05; n = 10^4; E = 1 / 20   # regular orbits: B(n) is the bound for the regular test orbits of the preamble
    for m in (KahanLi(Float64, 8), ForestRuth(Float64), Leapfrog(Float64, :KDK))
        for _ in 1:5
            a, b = hh_interval(E); pm = sqrt(2E)
            z0 = nothing
            while z0 === nothing
                z0 = hh_start(E, a + (b - a) * rand(rng), pm * (2 * rand(rng) - 1))
            end
            zn = advance_n(m, flow, z0, tau, n)
            zb = hh_flipp(advance_n(m, flow, hh_flipp(zn), tau, n))
            @test maximum(abs, zb .- z0) <= hh_B(Float64, 1.0, n, tau)
        end
    end
end

@testset "AC4a energy-surface starts: energy, sign of p1, section domain" begin
    _ = (random_energy_start, section_interval, energy_start)
    flow = HenonHeilesFlow(); et = hh_eps(Float64)
    for (k, E) in enumerate((1 / 12, 1 / 8, 1 / 6))
        rng = StableRNG(hh_seed(4, k)); a, b = section_interval(flow, E); pm = sqrt(2E)
        nbad_e = 0; nbad_s = 0; nbad_p = 0; nbad_d = 0; nbad_q = 0
        for _ in 1:10^5
            z, st = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
            st === :ok || (nbad_s += 1; continue)
            z isa SVector{4, Float64} || (nbad_s += 1; continue)
            abs(hh_H(z) - E) <= 4 * et * max(1, E) || (nbad_e += 1)
            z[3] > 0 || (nbad_p += 1)
            (z[1] == 0 && a <= z[2] <= b && z[4]^2 / 2 + hh_u(z[2]) <= E + 4 * et) || (nbad_d += 1)
        end
        @test nbad_s == 0 && nbad_e == 0 && nbad_p == 0 && nbad_d == 0
        # the deterministic constructor of Equation 8 agrees with the test-side formula
        for (q2, p2) in ((0.0, 0.0), (0.1, 0.1), (-0.02, 0.2))
            hh_start(E, q2, p2) === nothing && continue
            z, st = energy_start(flow, E, q2, p2)
            @test st === :ok && z isa SVector{4, Float64} && abs(z[3] - hh_start(E, q2, p2)[3]) <= 2 * et
        end
    end
end

@testset "AC4b section interval and the area of the section domain" begin
    _ = (section_interval, energy_grid)
    flow = HenonHeilesFlow(); et = hh_eps(Float64)
    for E in (1 / 12, 1 / 8)
        a, b = section_interval(flow, E); ra, rb = hh_interval(E)
        @test abs(a - ra) <= 1.0e-14 && abs(b - rb) <= 1.0e-14
    end
    a, b = section_interval(flow, 1 / 6)
    @test abs(a + 1 / 2) <= 4 * et && abs(b - 1) <= 4 * et                              # Equation 7 at the escape energy
    # harmonic limit: A(E) / (2 pi E) -> 1 (Equation 8), by the test-side quadrature; the defect shrinks with E
    d1 = abs(hh_area(1.0e-3) / (2pi * 1.0e-3) - 1); d2 = abs(hh_area(1.0e-4) / (2pi * 1.0e-4) - 1)
    @test d1 < 0.05 && d2 < 0.5 * d1
    # empirical area fraction of the grid against the quadrature, relative 1e-3 (a 1000 x 1000 grid), and the properties of every grid start
    for E in (1 / 12, 1 / 8, 1 / 6)
        zs = energy_grid(flow, E, 1000); a, b = section_interval(flow, E); pm = sqrt(2E)
        frac = length(zs) * (b - a) * 2pm / 1000^2 / hh_area(E)
        @test abs(frac - 1) <= 1.0e-3
        @test all(z -> z[1] == 0 && z[3] > 0 && abs(hh_H(z) - E) <= 4 * et * max(1, E), zs)
    end
    zs = energy_grid(flow, 1 / 8, 10)                                                     # the grid is deterministic
    @test zs == energy_grid(flow, 1 / 8, 10) && eltype(zs) == SVector{4, Float64}
end

@testset "AC4c infeasible requests return a status, never a value" begin
    _ = (random_energy_start, energy_start)
    flow = HenonHeilesFlow(); rng = StableRNG(hh_seed(4, 5))
    for E in (1 / 12, 1 / 8, 1 / 6, 0.2)
        nbad = 0
        for _ in 1:1000                                                                  # p2 in (2, 3): p2^2/2 >= 2 exceeds E for every energy used here
            z, st = random_energy_start(flow, E, (-0.2, 0.2), (2.0, 3.0), rng; maxtries = 100)
            (z === nothing && st !== :ok) || (nbad += 1)
        end
        @test nbad == 0
        z, st = energy_start(flow, E, 0.1, 2.0)
        @test z === nothing && st === :outside
    end
    z, st = energy_start(flow, NaN, 0.1, 0.1); @test z === nothing && st === :nonfinite
    z, st = energy_start(flow, 1 / 8, NaN, 0.1); @test z === nothing && st === :nonfinite
    # above the escape energy the start is accepted only inside r <= 1 (Algorithm 2): with q1 = 0, |q2| <= 1
    for E in (0.2, 0.25)
        for _ in 1:200
            z, st = random_energy_start(flow, E, (-3.0, 3.0), (-1.0, 1.0), rng)
            @test st === :ok && abs(z[2]) <= 1 && abs(hh_H(z) - E) <= 4 * hh_eps(Float64) * max(1, E)
        end
    end
    @test energy_start(flow, Double64(1) / 8, Double64(1) / 10, Double64(0))[1] isa SVector{4, Double64}   # generic in the element type
end

@testset "AC5 bound below the escape energy: r <= r_b(E), |p|^2 <= 2E, no escape" begin
    _ = (random_energy_start, section_interval, escape_orbit)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; T = 10^3
    nstart = hh_size(40, 200)
    for (k, E) in enumerate((1 / 12, 1 / 8, 1 / 6))
        rng = StableRNG(hh_seed(5, k)); a, b = section_interval(flow, E); pm = sqrt(2E); rb = hh_rb(E)
        rmax = 0.0; pmax = 0.0; nesc = 0
        for _ in 1:nstart
            z, st = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
            st === :ok || continue
            escape_orbit(m, flow, z, tau, round(Int, T / tau)).status === :bound || (nesc += 1)
            for _ in 1:round(Int, T / tau)
                z = advance(m, flow, z, tau)
                rmax = max(rmax, hh_rad(z)); pmax = max(pmax, z[3]^2 + z[4]^2)
            end
        end
        @test rmax <= rb * (1 + 1.0e-6)
        @test pmax <= 2E * (1 + 1.0e-6)
        @test nesc == 0
    end
end
