# AC9 to AC12: structure (symplecticity, reversibility, tangent consistency, tangent against reference).
# The one-step map of a method is IntegratorMap(m, flow, tau) (a Systems map: forward, inverse, tangent, inverse_tangent, dimension); see the test plan.

const IC_SYMFL = ((:f1_3, 0.05), (:f4, 0.05), (:f2, 0.01))         # (flow, tau) pairs for the structure tests
ic_imap(m, flow, tau) = IntegratorMap(m, flow, tau)
ic_jacobian(f, z) = ForwardDiff.jacobian(w -> Vector(f(SVector{length(z)}(w))), collect(z))
# central difference Jacobian with step h (BigFloat tests)
function ic_fdjac(f, z, h)
    d = length(z); J = Matrix{eltype(z)}(undef, d, d)
    for j in 1:d
        e = zeros(eltype(z), d); e[j] = h
        J[:, j] = (collect(f(SVector{d}(collect(z) .+ e))) .- collect(f(SVector{d}(collect(z) .- e)))) ./ (2h)
    end
    return J
end

@testset "AC9 symplecticity of the one-step Jacobian" begin
    _ = (IntegratorMap, tangent, symplectic_residual, Gauss, TaoExtended)
    rng = StableRNG(ic_seed(9, 1))
    for T in (Float64, Double64), M in IC_ALLX, (kind, tau) in IC_SYMFL
        m = M.mk(T); flow = ic_flow(T, kind); d = ic_dim(kind); worst = 0.0
        for _ in 1:ic_size(10, 1000)                                                    # D38: 10 points per case (spec 1000)
            z = ic_to(T, ic_draw(rng, kind))
            Mx = Matrix(tangent(ic_imap(m, flow, T(tau)), z))
            r = Float64(maximum(abs, Mx' * Matrix(symplectic_form(T, d)) * Mx - symplectic_form(T, d)))
            worst = max(worst, r / (8 * d * ic_eps(T) * Float64(maximum(abs, Mx))^2))
        end
        @test worst <= 1
    end
    for T in (Float64, Double64), (s, mk) in IC_GAUSS, kind in (:f1_3, :f4)           # implicit: plus the solver contribution 64 eps d
        m = mk(T); flow = ic_flow(T, kind); d = ic_dim(kind)
        for _ in 1:ic_size(5, 100)
            z = ic_to(T, ic_draw(rng, kind))
            Mx = Matrix(tangent(ic_imap(m, flow, T(0.1)), z))
            r = Float64(maximum(abs, Mx' * Matrix(symplectic_form(T, d)) * Mx - symplectic_form(T, d)))
            @test r <= 8 * d * ic_eps(T) * Float64(maximum(abs, Mx))^2 + 64 * ic_eps(T) * d
        end
    end
    # Tao: symplectic in the extended space with the extended form (Tao2016, Sec. I); the physical block residual is reported
    for T in (Float64, Double64), order in (2, 4), omega in (10.0, 100.0)
        m = TaoExtended(T, T(omega); order = order); flow = ic_flow(T, :t1); Om = ic_tao_omega(T, 1)
        for _ in 1:ic_size(5, 100)
            z = extend_state(ic_to(T, ic_draw(rng, :t1)))
            Mx = Matrix(tangent(ic_imap(m, flow, T(0.01)), z))
            r = Float64(maximum(abs, Mx' * Om * Mx - Om))
            @test r <= 8 * 4 * ic_eps(T) * Float64(maximum(abs, Mx))^2
        end
    end
end

@testset "AC10 reversibility, adjoint property" begin
    _ = (advance_n, adjoint_method, is_palindromic, Gauss, TaoExtended)
    rng = StableRNG(ic_seed(10, 1))
    for M in IC_ALLX, (kind, tau) in IC_SYMFL
        m = M.mk(Float64); flow = ic_flow(Float64, kind); nH = round(Int, 10 / tau)
        for _ in 1:ic_size(3, 100)
            z = ic_to(Float64, ic_draw(rng, kind)); Mn = ic_mmax(z)
            if is_palindromic(m)
                @test M.sym !== false                                                          # a table declared non-symmetric must not report palindromic
                for n in (1, nH)
                    back = advance_n(m, flow, advance_n(m, flow, z, tau, n), -tau, n)
                    @test maximum(abs, back .- z) <= 16 * ic_B(Float64, Mn, n, tau) * 4      # the orbit stays in a compact set; M_n of the start times 4
                end
                # S* = S as maps for a symmetric table
                a = advance(adjoint_method(m), flow, z, tau); b = advance(m, flow, z, tau)
                @test maximum(abs, collect(a) .- collect(b)) <= 16 * ic_eps(Float64) * Mn^2
            else
                # adjoint property: S*(tau) S(-tau) = id for a general table (I8; S1: the adjoint of A then B is B then A)
                back = advance(adjoint_method(m), flow, advance(m, flow, z, -tau), tau)
                @test maximum(abs, back .- z) <= 16 * ic_eps(Float64) * Mn^2
                @test M.sym !== true                                                           # a table declared palindromic must report palindromic
            end
        end
    end
    z = ic_to(Float64, IC_START[:f4])
    a = advance(SymplecticEuler(Float64, :KD), ic_flow(Float64, :f4), z, 0.1); b = advance(adjoint_method(SymplecticEuler(Float64, :DK)), ic_flow(Float64, :f4), z, 0.1)
    @test a == b                                                                          # the adjoint of S1 (A then B) is S1 (B then A)
    for (s, mk) in IC_GAUSS, (kind, tau) in IC_SYMFL[1:2]                                 # implicit: to the solver tolerance
        m = mk(Float64); flow = ic_flow(Float64, kind); nH = round(Int, 10 / tau)
        for _ in 1:ic_size(2, 50)
            z = ic_to(Float64, ic_draw(rng, kind)); Mn = ic_mmax(z)
            for n in (1, nH)
                back = advance_n(m, flow, advance_n(m, flow, z, tau, n), -tau, n)
                @test maximum(abs, back .- z) <= 16 * ic_B(Float64, Mn, n, tau) * 4 + 16 * 8 * ic_eps(Float64) * Mn * n
            end
        end
    end
    # Tao Strang: symmetric in the extended space
    flow = ic_flow(Float64, :t1); m = TaoExtended(Float64, 10.0; order = 2)
    for _ in 1:ic_size(5, 100)
        z = extend_state(ic_to(Float64, ic_draw(rng, :t1)))
        back = advance_n(m, flow, advance_n(m, flow, z, 0.01, 1000), -0.01, 1000)
        @test maximum(abs, back .- z) <= 16 * ic_B(Float64, 2.0, 1000, 0.01) * 4
    end
end

# D112 (runtime, D38): ForwardDiff compiles a new closure and tag per method and flow (0.3 s each, 88 s of the 3 min 12 s of the suite for the 41 explicit methods). The package run uses one
# representative of every code path of the tangent (drift-first and kick-first tables, non-palindromic S1 and S4, merged kicks S2KDK, the longest table S13, S5 with negative weights,
# S9 with kick first, the generated SABA and SBAB, the two corrector kinds, the generalized order table); IC_FULL_RUN (verify/integrator-catalog.jl) runs all 41 with the same assertions.
const IC_AC11_IDS = ("S1DK", "S1KD", "S2DKD", "S2KDK", "S3", "S4", "S5", "S9", "S12-6", "S13", "SABA3", "SBAB3", "SABAC2", "SBABC2", "ABA864")
ic_ac11_set() = IC_FULL ? IC_ALLX : filter(M -> M.id in IC_AC11_IDS, IC_ALLX)

@testset "AC11 tangent step equals the derivative of the one-step map (ForwardDiff, 256-bit finite difference)" begin
    _ = (IntegratorMap, tangent, tangent_advance, stage_count, advance, advance_n)
    rng = StableRNG(ic_seed(11, 1)); set11 = ic_ac11_set()
    # binary64 and double-double: ForwardDiff; the bound 16 eps (1 + |M|) is scaled by the stage count k (rounding of k elementary Jacobian products), see the plan
    for T in (Float64, Double64), M in set11, (kind, tau) in IC_SYMFL[2:3]
        m = M.mk(T); flow = ic_flow(T, kind); d = ic_dim(kind); k = max(1, stage_count(m))
        for _ in 1:ic_size(3, 1000)
            z = ic_to(T, ic_draw(rng, kind)); h = T(tau)
            A = Matrix(tangent(ic_imap(m, flow, h), z)); B = ic_jacobian(w -> advance(m, flow, w, h), z)
            @test Float64(maximum(abs, A .- B)) <= 16 * ic_eps(T) * (1 + Float64(maximum(abs, A))) * k
        end
    end
    # n steps: the product of tangent steps equals the tangent of the n-step map
    for M in set11, (n, kind, tau) in ((1, :f4, 0.05), (10, :f4, 0.05), (10, :f2, 0.01))
        m = M.mk(Float64); flow = ic_flow(Float64, kind); d = ic_dim(kind); k = max(1, stage_count(m))
        z = ic_to(Float64, ic_draw(rng, kind)); W = SMatrix{d, d, Float64}(I); zz = z
        for _ in 1:n
            zz, W = tangent_advance(m, flow, zz, W, tau)
        end
        B = ic_jacobian(w -> advance_n(m, flow, w, tau, n), z)
        @test maximum(abs, Matrix(W) .- B) <= 16 * ic_eps(Float64) * n * (1 + maximum(abs, W))^2 * k
    end
    # 256-bit: central difference with step 2^-100; bound 2^-130 (1 + |M|) (rounding eps/h = 2^-155, truncation 2^-200); implicit methods by the same difference
    setprecision(BigFloat, IC_BITS) do
        h = BigFloat(2)^(-100)
        ms = [(M.id, M.mk(BigFloat)) for M in set11]
        append!(ms, [("G$s", mk(BigFloat)) for (s, mk) in IC_GAUSS])
        for (id, m) in ms, (kind, tau) in (IC_SYMFL[2], IC_SYMFL[3])
            flow = ic_flow(BigFloat, kind); z = ic_to(BigFloat, ic_draw(rng, kind)); tb = BigFloat(tau)
            A = Matrix(tangent(ic_imap(m, flow, tb), z)); B = ic_fdjac(w -> advance(m, flow, w, tb), z, h)
            @test Float64(maximum(abs, A .- B)) <= 2.0^(-130) * (1 + Float64(maximum(abs, A)))
        end
    end
end

@testset "AC12a tangent against the exact matrix exponential of the linearized flow (F1)" begin
    _ = (tangent_advance, advance)
    kind = :f1_3; tau = 0.05; n = 200; Tend = 10.0
    exact = setprecision(BigFloat, IC_BITS) do
        hcat([ic_f1_exact(kind, [i == j ? 1.0 : 0.0 for i in 1:6], BigFloat(Tend)) for j in 1:6]...)
    end
    for M in IC_ALLX
        m = M.mk(Float64); flow = ic_flow(Float64, kind)
        z = SVector{6, Float64}(zeros(6)); W = SMatrix{6, 6, Float64}(I)
        for _ in 1:n
            z, W = tangent_advance(m, flow, z, W, tau)
        end
        # truncation of the method itself: the same product in double-double against the exact matrix (the truncation is far above 2^-100 here)
        trunc = begin
            mb = M.mk(Double64); fb = ic_flow(Double64, kind)
            zb = SVector{6, Double64}(zeros(Double64, 6)); Wb = SMatrix{6, 6, Double64}(Matrix{Double64}(I, 6, 6))
            for _ in 1:n
                zb, Wb = tangent_advance(mb, fb, zb, Wb, Double64(tau))
            end
            Float64(maximum(abs, ic_big.(Matrix(Wb)) .- exact))
        end
        err = Float64(maximum(abs, ic_big.(Matrix(W)) .- exact))
        @test err <= ic_B(Float64, 1.0 + maximum(abs, W), n, tau) + 10 * trunc
    end
end

# reference tangent columns by central finite difference of Gauss s = 8 orbits (step 1/64) at `bits` bits; relative to the largest entry
function ic_fd_columns(z0, Tend, cols, delta, bits)
    return setprecision(BigFloat, bits) do
        flow = ic_flow(BigFloat, :f2); g = Gauss(BigFloat, 8; maxiter = 500, tol = 1e-45); out = Dict{Int, Vector{BigFloat}}()
        for j in cols
            e = zeros(BigFloat, 4); e[j] = delta
            zp = advance_n(g, flow, ic_to(BigFloat, z0) .+ SVector{4}(e), BigFloat(1) / 64, round(Int, Tend * 64))
            zm = advance_n(g, flow, ic_to(BigFloat, z0) .- SVector{4}(e), BigFloat(1) / 64, round(Int, Tend * 64))
            out[j] = Vector(zp .- zm) ./ (2 * delta)
        end
        out
    end
end

@testset "AC12b tangent of S2 and of a 6th order table against finite differences of reference orbits (F2, e = 0.5, T = 10)" begin
    _ = (tangent_advance, advance_n, Gauss)
    z0 = IC_START[:f2]; Tend = 10.0
    cols = ic_size([1, 4], [1, 2, 3, 4])                                                   # D38: two columns (two position/momentum pairs of orbits in 192 bits)
    fd = ic_fd_columns(z0, Tend, cols, BigFloat(2)^(-40), 192)
    scale = Float64(maximum(maximum(abs, v) for v in values(fd)))
    frame(mk, tau) = begin
        m = mk(Float64); flow = ic_flow(Float64, :f2); z = ic_to(Float64, z0); W = SMatrix{4, 4, Float64}(I)
        for _ in 1:round(Int, Tend / tau)
            z, W = tangent_advance(m, flow, z, W, tau)
        end
        W
    end
    relerr(W) = maximum(Float64(maximum(abs, ic_big.(collect(W[:, j])) .- fd[j])) for j in cols) / scale
    # D112: the tolerance follows the order. The tangent of a method of order p differs from the reference tangent by its truncation C tau^p (S2 at 2^-16: 3.9e-9 measured, tau^2 = 2.3e-10 times the
    # growth of the tangent over T = 10, far above 2^-30). For S2 (symmetric, order 2, even error expansion) the order law is asserted (slope 2 between 2^-15 and 2^-16) and the Richardson
    # extrapolate (4 W(tau/2) - W(tau))/3, whose error is O(tau^4), meets 2^-30; the 6th order table meets 2^-30 directly.
    W1 = frame(T -> Leapfrog(T, :DKD), ldexp(1.0, -15)); W2 = frame(T -> Leapfrog(T, :DKD), ldexp(1.0, -16))
    r1 = relerr(W1); r2 = relerr(W2); rx = relerr((4 .* W2 .- W1) ./ 3)
    @info "AC12b reported: S2 relative difference to the reference tangent, tau = 2^-15, 2^-16, extrapolated" r1 r2 rx log2(r1 / r2)
    @test abs(log2(r1 / r2) - 2) <= 0.3
    @test rx <= 2.0^(-30)
    rel6 = relerr(frame(T -> KahanLi(T, 6), 1 / 128))
    @info "AC12b reported: S12-6 relative difference to the reference tangent" rel6
    @test rel6 <= 2.0^(-30)
end

@testset "AC12c frame over n steps against the two-trajectory finite difference of the method's own orbits (256-bit, perturbation 2^-40)" begin
    _ = (tangent_advance, advance_n)
    z0 = IC_START[:f2]; Tend = 10.0; tau = 1 / 64; n = 640
    for (id, mk) in (("S2", T -> Leapfrog(T, :DKD)), ("S12-6", T -> KahanLi(T, 6)), ("S5", T -> ForestRuth(T)))
        m = mk(Float64); flow = ic_flow(Float64, :f2)
        z = ic_to(Float64, z0); W = SMatrix{4, 4, Float64}(I)
        for _ in 1:n
            z, W = tangent_advance(m, flow, z, W, tau)
        end
        rel = setprecision(BigFloat, IC_BITS) do
            mb = mk(BigFloat); fb = ic_flow(BigFloat, :f2); delta = BigFloat(2)^(-40); err = 0.0; sc = 0.0
            for j in 1:4
                e = zeros(BigFloat, 4); e[j] = delta
                zp = advance_n(mb, fb, ic_to(BigFloat, z0) .+ SVector{4}(e), BigFloat(tau), n)
                zm = advance_n(mb, fb, ic_to(BigFloat, z0) .- SVector{4}(e), BigFloat(tau), n)
                col = Vector(zp .- zm) ./ (2 * delta)
                err = max(err, Float64(maximum(abs, ic_big.(collect(W[:, j])) .- col))); sc = max(sc, Float64(maximum(abs, col)))
            end
            err / sc
        end
        @test rel <= 2.0^(-30)
    end
end
