# AC6 to AC9: energy conservation along the orbit, the surface-of-section map (crossings, grazing, Jacobian, symplecticity of the sub-step).
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).
# AC7 and AC8 compare with a Double64 reference (package oracle O3, hh_ref_*) only up to T_AGREE: a chaotic orbit amplifies the 1e-14 level differences of two
# integrators by exp(lambda t), so a pointwise 1e-8 agreement is meaningful for a bounded time only (test plan, spec point U1); beyond it only the invariants are checked.

const HH_TAGREE = 60.0

@testset "AC6 energy conservation of the order-8 and order-4 tables" begin
    _ = (random_energy_start, section_interval)
    flow = HenonHeilesFlow(); E = 1 / 8; a, b = section_interval(flow, E); pm = sqrt(2E)
    rng = StableRNG(hh_seed(6, 1)); nstart = hh_size(10, 200); T = 10^3
    starts = [random_energy_start(flow, E, (a, b), (-pm, pm), rng)[1] for _ in 1:nstart]
    @test all(z -> z isa SVector{4, Float64}, starts)
    for (name, m, bound) in (("order-8", KahanLi(Float64, 8), 1.0e-8), ("order-4", ForestRuth(Float64), 1.0e-4))
        worst = Dict{Float64, Float64}()
        for tau in (1 / 20, 1 / 40, 1 / 80)
            w = 0.0
            for z0 in starts
                H0 = hh_H(z0); z = z0
                for _ in 1:round(Int, T / tau)
                    z = advance(m, flow, z, tau)
                    w = max(w, abs(hh_H(z) - H0))
                end
            end
            worst[tau] = w
        end
        @test worst[1 / 20] <= bound                                                       # no secular drift: the bound holds over the whole run
        @info "AC6 reported: max |H - H0| by step size (scaling with tau^p)" name worst
    end
end

@testset "AC7a crossings of the plane q1 = 0: invariants, count and location against the reference" begin
    _ = (random_energy_start, section_interval, section_map)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; et = hh_eps(Float64)
    Trun = hh_size(100, 1000); nst = hh_size(4, 20)
    for (k, E) in enumerate((1 / 12, 1 / 8, 1 / 6))
        rng = StableRNG(hh_seed(7, k)); a, b = section_interval(flow, E); pm = sqrt(2E)
        for _ in 1:nst
            z0, st = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
            @test st === :ok
            rec = section_map(m, flow, z0, tau, round(Int, Trun / tau))
            @test rec.status === :ok
            @test length(rec.t) >= 2 && rec.t[1] == 0 && rec.z[1] == z0 && issorted(rec.t; lt = <=)   # the start is the first section point; times strictly increase
            @test all(z -> abs(z[1]) <= 16 * et, rec.z[2:end])                              # I6: |q1| at the root
            @test all(z -> z[3] > 0, rec.z)                                                  # I6: p1 > 0
            @test all(z -> abs(hh_H(z) - E) <= 1.0e-9, rec.z)                                # I6: in D(E) up to the energy error
            @test length(rec.t) == length(rec.z)
            ref = hh_ref_section(z0, HH_TAGREE)
            idx = findall(<=(HH_TAGREE), rec.t)
            @test length(idx) == length(ref.t)                                               # no missed and no extra crossing
            if length(idx) == length(ref.t)
                @test maximum(abs.(rec.t[idx] .- ref.t)) <= 1.0e-8
                @test maximum(abs.([rec.z[i][2] for i in idx] .- ref.q2)) <= 1.0e-8
                @test maximum(abs.([rec.z[i][4] for i in idx] .- ref.p2)) <= 1.0e-8
            end
        end
    end
end

@testset "AC7b grazing is flagged, not recorded; one crossing solve" begin
    _ = (section_map, locate_crossing)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; et = hh_eps(Float64)
    # q1 = -1e-40, p1 = 0: q1 returns to 0 with p1 = 1e-40 < 8 eps at the root (a grazing bracket)
    rec = section_map(m, flow, SVector(-1.0e-40, 0.3, 0.0, 0.0), tau, 60)
    @test rec.grazing >= 1 && isempty(rec.t) && rec.status === :ok
    # a normal bracket: the solved state has |q1| <= 16 eps, 0 < s <= tau, and is the step of size s from the bracket start
    rng = StableRNG(hh_seed(7, 20)); E = 1 / 8; a, b = section_interval(flow, E); pm = sqrt(2E)
    for _ in 1:5
        z0, _ = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
        zb = hh_bracket(m, flow, z0, tau)
        zc, s = locate_crossing(m, flow, zb, tau)
        @test 0 < s <= tau && abs(zc[1]) <= 16 * et && zc[3] > 0
        @test zc == advance(m, flow, zb, s)
    end
end

@testset "AC8 section Jacobian: reference finite differences, determinant" begin
    _ = (random_energy_start, section_interval, section_map)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; et = hh_eps(Float64)
    dtol = hh_defaults()["det_residual_tol"]; im = IntegratorMap(m, flow, tau)
    for (k, E) in enumerate((1 / 12, 1 / 8, 1 / 6))
        rng = StableRNG(hh_seed(8, k)); a, b = section_interval(flow, E); pm = sqrt(2E)
        for _ in 1:hh_size(2, 10)
            z0, _ = random_energy_start(flow, E, (a, b), (-pm, pm), rng); E0 = hh_H(z0)
            rec = section_map(m, flow, z0, tau, round(Int, 30 / tau); jacobian = true)
            @test length(rec.J) == length(rec.t) && rec.J[1] == one(SMatrix{2, 2, Float64})  # the start leg is the identity
            for i in 2:min(length(rec.t), 4)
                J = rec.J[i]; zp = rec.z[i - 1]
                Jref = hh_ref_jacobian(E0, zp[2], zp[4])                                      # leg from the previous section point, reference return map
                @test maximum(abs.(J .- Jref) ./ max.(1, abs.(Jref))) <= 1.0e-6
                nsteps = round(Int, (rec.t[i] - rec.t[i - 1]) / tau) + 1
                M = max(1.0, maximum(abs, tangent(im, zp)), maximum(abs, J))                  # Amendment D117: M_max includes the leg's own product (|J| up to about 23 at E = 1/6)
                @test abs(det(J) - 1) <= dtol * et * (1 + M) * nsteps + 4 * et
            end
        end
    end
end

@testset "AC9 symplecticity of the section map with the sub-step: growth of the determinant defect" begin
    _ = (energy_start, section_map)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05
    z0, st = energy_start(flow, 1 / 8, -0.015, 0.0)                                           # orbit C (chaotic, Equation 9, p1 rebuilt)
    @test st === :ok
    rec = section_map(m, flow, z0, tau, round(Int, 10^4 / tau); jacobian = true)
    @test length(rec.J) >= 1001
    d = [det(J) - 1 for J in rec.J[2:1001]]
    c = abs.(cumsum(d)); idx = [i for i in 10:1000 if c[i] > 0]                                # signed running sum: 0.5 (random walk) to 1 (bias)
    s_signed = hh_slope(log.(idx), log.(c[idx]))
    ca = cumsum(abs.(d)); s_abs = hh_slope(log.(10:1000), log.(ca[10:1000]))
    @info "AC9 reported: slope of the cumulative determinant defect, signed and absolute sums; mean |det J - 1|" s_signed s_abs hh_mean(abs.(d))
    @test s_signed <= 1.2
    @test s_abs <= 1.2
end
