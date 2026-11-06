# FL: the Hamiltonian-flow side added to Systems by this module (open point P1): F1 (HarmonicOscillator), F4 (Pendulum), F2 (Kepler), the test-only
# T1 (TaoNonseparable) and T2 (HenonHeilesFlow). Not an acceptance criterion by itself: the integrators cannot be tested without these flows (spec Scope).
# Each testset begins by naming the new symbols it uses, so a missing module errors that testset (UndefVarError) and no top-level code calls the module.

const IC_SEP = (:f1_1, :f1_3, :f4, :f2, :t2)      # separable test flows (T = |p|^2/2)

# test-side partial flows, written out by hand (independent of the module)
function ic_drift(kind, z, t)
    h = length(z) ÷ 2
    return vcat(z[1:h] .+ t .* z[(h + 1):end], z[(h + 1):end])
end
function ic_gradV(kind, q; eps = 1.0)
    kind === :f4 && return [eps * sin(q[1])]
    kind === :f2 && return q ./ sqrt(q[1]^2 + q[2]^2)^3
    kind === :t2 && return [q[1] + 2 * q[1] * q[2], q[2] + q[1]^2 - q[2]^2]
    ws = ic_f1w(kind); T = eltype(q)                          # w^2 in the element type of q (D112)
    return [T(ws[j])^2 * q[j] for j in eachindex(ws)]
end
function ic_kick(kind, z, t; eps = 1.0)
    h = length(z) ÷ 2
    return vcat(z[1:h], z[(h + 1):end] .- t .* ic_gradV(kind, z[1:h]; eps = eps))
end

@testset "FL1 dimension, element type and Hamiltonian" begin
    _ = (HarmonicOscillator, Pendulum, Kepler, TaoNonseparable, HenonHeilesFlow, hamiltonian, dimension)
    rng = StableRNG(ic_seed(0, 1))
    for kind in (IC_SEP..., :t1), T in (Float64, Double64, BigFloat)
        flow = ic_flow(T, kind; eps = 0.25)
        @test dimension(flow) == ic_dim(kind)
        for _ in 1:20
            z0 = ic_draw(rng, kind)
            z = ic_to(T, z0)
            H = hamiltonian(flow, z)
            @test H isa T                                         # element type preserved
            ref = ic_H(kind, ic_to(BigFloat, z0); eps = 0.25)
            @test ic_err([H], [ref]) <= 64 * ic_eps(T) * ic_mmax(z0)^3
        end
    end
end

@testset "FL2 partial flows are exact: drift and kick" begin
    _ = (drift, kick, HarmonicOscillator)
    rng = StableRNG(ic_seed(0, 2))
    for kind in IC_SEP, T in (Float64, Double64)
        flow = ic_flow(T, kind; eps = 0.25)
        for _ in 1:20
            z0 = ic_draw(rng, kind); z = ic_to(T, z0); t = T(0.3)
            M = ic_mmax(z0)
            @test ic_err(drift(flow, z, t), ic_drift(kind, ic_to(BigFloat, z0), BigFloat(0.3))) <= 8 * ic_eps(T) * M
            @test ic_err(kick(flow, z, t), ic_kick(kind, ic_to(BigFloat, z0), BigFloat(0.3); eps = 0.25)) <= 16 * ic_eps(T) * M^3
            # group property of an exact flow, either sign
            @test maximum(abs, drift(flow, drift(flow, z, t), -t) .- z) <= 8 * ic_eps(T) * M
            @test maximum(abs, kick(flow, kick(flow, z, t), -t) .- z) <= 16 * ic_eps(T) * M^3
            @test maximum(abs, drift(flow, drift(flow, z, t), t) .- drift(flow, z, 2t)) <= 8 * ic_eps(T) * M
        end
    end
end

@testset "FL3 elementary Jacobians equal the derivative of the partial flows" begin
    _ = (drift_jacobian, kick_jacobian, drift, kick)
    rng = StableRNG(ic_seed(0, 3))
    for kind in IC_SEP
        flow = ic_flow(Float64, kind; eps = 0.25)
        for _ in 1:20
            z = ic_to(Float64, ic_draw(rng, kind)); t = 0.3
            Jd = ForwardDiff.jacobian(w -> Vector(drift(flow, SVector{length(z)}(w), t)), collect(z))
            Jk = ForwardDiff.jacobian(w -> Vector(kick(flow, SVector{length(z)}(w), t)), collect(z))
            @test maximum(abs, Matrix(drift_jacobian(flow, z, t)) .- Jd) <= 16 * ic_eps(Float64) * (1 + maximum(abs, Jd))
            @test maximum(abs, Matrix(kick_jacobian(flow, z, t)) .- Jk) <= 16 * ic_eps(Float64) * (1 + maximum(abs, Jk))
        end
    end
end

@testset "FL4 gradient and vector field of H, including the nonseparable T1" begin
    _ = (hamiltonian_gradient, vector_field)
    rng = StableRNG(ic_seed(0, 4))
    for kind in (IC_SEP..., :t1)
        flow = ic_flow(Float64, kind; eps = 0.25)
        for _ in 1:20
            z0 = ic_draw(rng, kind); z = ic_to(Float64, z0)
            g = ForwardDiff.gradient(w -> ic_H(kind, w; eps = 0.25), z0)
            @test maximum(abs, collect(hamiltonian_gradient(flow, z)) .- g) <= 32 * ic_eps(Float64) * (1 + maximum(abs, g))
            v = ic_vf(kind, z0; eps = 0.25)
            @test maximum(abs, collect(vector_field(flow, z)) .- v) <= 32 * ic_eps(Float64) * (1 + maximum(abs, v))
        end
    end
end

@testset "FL5 F1 exact solution and its matrix" begin
    _ = (exact_flow, exact_flow_matrix, HarmonicOscillator)
    rng = StableRNG(ic_seed(0, 5))
    for kind in (:f1_1, :f1_3), T in (Float64, Double64, BigFloat)
        flow = ic_flow(T, kind)
        for _ in 1:10
            z0 = ic_draw(rng, kind); t = 2.5
            ref = ic_f1_exact(kind, z0, BigFloat(t))
            @test ic_err(exact_flow(flow, ic_to(T, z0), T(t)), ref) <= 16 * ic_eps(T) * ic_mmax(z0) * (1 + t)
            Mx = Matrix(exact_flow_matrix(flow, T(t)))
            @test maximum(abs, Mx * ic_big.(Float64.(z0)) .- ref) <= 16 * ic_eps(T) * ic_mmax(z0) * (1 + t) # also the matrix acts as the flow
        end
    end
end

@testset "FL6 F2 invariants, perturbation scale, corrector force, force counter" begin
    _ = (angular_momentum, runge_lenz, perturbation_scale, corrector_force, corrector_hessian, CountingFlow, kick_count, reset_counts!)
    rng = StableRNG(ic_seed(0, 6))
    K = ic_flow(Float64, :f2)
    for _ in 1:20
        z0 = ic_draw(rng, :f2); z = ic_to(Float64, z0)
        @test abs(angular_momentum(K, z) - ic_L(z0)) <= 8 * ic_eps(Float64) * ic_mmax(z0)^2
        # Runge-Lenz vector A = p x L - q/|q| (planar: (p2 L, -p1 L) - q/r), a first integral of the continuous flow: magnitude e
        r = sqrt(z0[1]^2 + z0[2]^2); L = ic_L(z0)
        A = [z0[4] * L - z0[1] / r, -z0[3] * L - z0[2] / r]
        @test maximum(abs, collect(runge_lenz(K, z)) .- A) <= 64 * ic_eps(Float64) * ic_mmax(z0)^3
        # each elementary flow preserves L (spec I6, derived)
        @test abs(angular_momentum(K, drift(K, z, 0.3)) - angular_momentum(K, z)) <= 16 * ic_eps(Float64) * ic_mmax(z0)^2
        @test abs(angular_momentum(K, kick(K, z, 0.3)) - angular_momentum(K, z)) <= 16 * ic_eps(Float64) * ic_mmax(z0)^2
    end
    P = Pendulum(0.25)
    @test perturbation_scale(P) == 0.25 && perturbation_scale(K) == 1.0
    # corrector potential W = |grad V|^2 with V the unscaled perturbation (V = -cos q): W = sin^2 q, grad W = sin 2q, Hessian 2 cos 2q (LaskarRobutel2001 Eq. (62) use)
    for q in (-1.3, 0.2, 0.9, 2.1)
        z = SVector(q, 0.1)
        @test only(corrector_force(P, z)) ≈ sin(2q) rtol = 1e-14
        @test only(corrector_hessian(P, z)) ≈ 2cos(2q) rtol = 1e-13
    end
    # F2: W = 1/r^4; the gradient follows from ForwardDiff of the test-side potential
    for _ in 1:10
        z0 = ic_draw(rng, :f2); q0 = z0[1:2]
        g = ForwardDiff.gradient(w -> 1 / (w[1]^2 + w[2]^2)^2, q0)
        @test maximum(abs, collect(corrector_force(K, ic_to(Float64, z0))) .- g) <= 64 * ic_eps(Float64) * (1 + maximum(abs, g))
    end
    # force counter: wraps a flow, counts the kick evaluations, forwards everything else
    cf = CountingFlow(P)
    @test kick_count(cf) == 0
    z = SVector(1.0, 0.3)
    @test kick(cf, z, 0.1) == kick(P, z, 0.1)
    kick(cf, z, 0.1)
    @test kick_count(cf) == 2
    reset_counts!(cf)
    @test kick_count(cf) == 0
    @test dimension(cf) == 2 && drift(cf, z, 0.1) == drift(P, z, 0.1) && kick_count(cf) == 0
end

@testset "FL7 T1 and T2 values" begin
    _ = (TaoNonseparable, HenonHeilesFlow, hamiltonian)
    # T1: H = (Q^2 + 1)(P^2 + 1)/2 (Tao2016 Sec. II A): hand values
    T1 = TaoNonseparable()
    @test hamiltonian(T1, SVector(0.0, 0.0)) == 0.5
    @test hamiltonian(T1, SVector(1.0, 1.0)) == 2.0
    @test hamiltonian(T1, SVector(1.0, 0.5)) == 1.25
    # T2: Henon-Heiles (spec item 18): H at the Hairer2008 start q = (0, 0.3), p2 = 0.2, p1 from H0 = 1/8 (hand value of p1^2 = 2 (1/8 - V - p2^2/2))
    T2 = HenonHeilesFlow()
    V = 0.3^2 / 2 - 0.3^3 / 3
    p1 = sqrt(2 * (1 / 8 - V - 0.2^2 / 2))
    @test hamiltonian(T2, SVector(0.0, 0.3, p1, 0.2)) ≈ 1 / 8 rtol = 1e-14
end
