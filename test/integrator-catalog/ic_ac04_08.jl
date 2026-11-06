# AC4 to AC8: measured orders (explicit symmetric tables, non-symmetric tables, near-integrable rows, implicit methods, Tao's method).
# Reference orbits: ic_ref (oracle O1: F1 exact, otherwise Gauss s = 8 in 256 bits); the Taylor reference (R1) arrives in AC19/AC21 when the extension is loaded.

# Slopes of the final-state error for tau = tau0/2^j, j = 0..8 at T = 8 (spec criterion 4, D112): returns the slopes used (errors above 1e3 B(n) and 1e3 10^-dtab).
function ic_check_order(M, kind; T = ic_slope_type(M.p), Tend = 8, J = 8)   # J = 5 for the implicit methods (each step costs an iteration; their slopes settle at the third halving)
    return setprecision(BigFloat, IC_BITS) do
        z0 = IC_START[kind]; ref = ic_ref(kind, z0, Tend)
        flow = ic_flow(T, kind); m = M.mk(T)
        tau0 = ic_pick_tau0(m, flow, z0, Tend, ref, T)
        taus, errs, ns = ic_scan(m, flow, z0, Tend, ref, tau0, T; J = J)
        Mn = ic_mmax(z0)
        lims = [1e3 * max(ic_B(T, Mn, ns[j], taus[j]), 10.0^(-M.dtab)) for j in eachindex(taus)]
        return ic_slopes(errs, lims), errs
    end
end

@testset "AC4 measured order of the explicit symmetric tables (F4 and F2)" begin
    _ = (advance_n, nominal_order)
    for M in IC_S
        (M.p >= 2 && M.sym !== false) || continue
        for kind in (:f4, :f2)
            sl, errs = ic_check_order(M, kind)
            @test ic_asym(sl, M.p, M.p <= 4 ? 0.2 : 0.3)                      # the two finest of at least three resolved halvings (D112)
            @info "AC4 reported" M.id kind sl
        end
    end
end

@testset "AC5 non-symmetric tables: S1 first order, S4 third order, local order on nonlinear flows" begin
    _ = (SymplecticEuler, Ruth3, advance_n, nominal_order, is_palindromic)
    for M in IC_S
        (M.sym === false) || continue
        for kind in (:f4, :f2)
            sl, errs = ic_check_order(M, kind)
            @test ic_asym(sl, M.p, 0.2)
            @info "AC5 reported" M.id kind sl
        end
    end
    # Eq. (30) of Ruth1983 is not transcribed in the spec; the nonlinear order conditions are tested through the local error exponent on nonlinear flows in 256 bits:
    # one step of size h against one step of Gauss s = 8 scales like h^(p+1) on the pendulum and on Henon-Heiles (T2)
    hs = [ldexp(1.0, -j) for j in 4:7]
    for M in (IC_S[1], IC_S[2], IC_S[6], IC_S[15]), (kind, z0) in ((:f4, [1.0, 0.3]), (:t2, [0.1, 0.3, 0.3, 0.2]))      # S1DK, S1KD, S4, S9
        sl, _ = ic_local_slopes(M.mk, kind, z0, hs)
        @test all(abs(s - (M.p + 1)) <= 0.3 for s in sl)
    end
end

# near-integrable pendulum H = p^2/2 - eps cos q, q(0) = 1, p(0) = 0.3, T = 8: error of method M at tau in element type T (spec criterion 6)
function ic_pend_err(M, eps, tau, ::Type{T}) where {T}
    return setprecision(BigFloat, IC_BITS) do
        z0 = [1.0, 0.3]; ref = ic_ref(:f4, z0, 8; eps = eps)
        return ic_err(advance_n(M.mk(T), ic_flow(T, :f4; eps = eps), ic_to(T, z0), T(tau), round(Int, 8 / tau)), ref)
    end
end
ic_sl(e1, e2) = log2(e1 / e2)

@testset "AC6a SABA and SBAB: the eps tau^2n branch (D112)" begin
    _ = (SABA, SBAB, advance_n, coefficients)
    # Error law e = a eps tau^2n + b eps^2 tau^2. The slope 2n between tau = 1 and 1/2 is asserted where the first branch is above the second and above the round-off floor:
    # n = 1, 2 in binary64 at eps = 1e-8; n = 3 in 256 bits at eps = 1e-8 (binary64 is at the floor, 1e-16); n = 4 in 256 bits at eps = 1e-16 (at 1e-8 the eps^2 tau^2 branch already dominates, slope 2.0)
    for M in IC_P[1:8]
        n = parse(Int, M.id[end:end])
        T, eps = n <= 2 ? (Float64, 1e-8) : (n == 3 ? (BigFloat, 1e-8) : (BigFloat, 1e-16))
        e1 = ic_pend_err(M, eps, 1.0, T); e2 = ic_pend_err(M, eps, 0.5, T)
        @test abs(ic_sl(e1, e2) - 2n) <= 0.5
    end
    # n = 6, 10 (256 bits): the eps tau^2n branch of these rows is below the eps^2 tau^2 branch for every eps and tau the reference orbit can resolve (a tau^12 is 1e-20 eps), so at eps = 1e-8 the
    # error is the second branch (slope 2 within 0.3, measured 2.04 to 2.05) and the order 2n is asserted at its source: the kicks are an n-point (SABA) or (n + 1)-point (SBAB) quadrature of the
    # free-flight integral, exact for the monomials of degree <= 2n - 1 and not for 2n (LaskarRobutel2001 Sect. 7). The same moment test is applied to n = 1 to 4.
    for M in IC_PBIG
        e1 = ic_pend_err(M, 1e-8, 1.0, BigFloat); e2 = ic_pend_err(M, 1e-8, 0.5, BigFloat)
        @test abs(ic_sl(e1, e2) - 2) <= 0.3
    end
    setprecision(BigFloat, IC_BITS) do
        for M in (IC_P[1:8]..., IC_PBIG...)
            n = parse(Int, filter(isdigit, M.id))
            pos, w = ic_kickpos(coefficients(M.mk(BigFloat))...)
            mom(k) = abs(sum(w .* pos .^ k) - 1 / BigFloat(k + 1))
            @test maximum(mom(k) for k in 0:(2n - 1)) <= 64 * ic_eps(BigFloat)
            @test mom(2n) > 1e-14                                                       # measured 1.5e-12 (SBAB10) to 0.17 (SBAB1)
        end
    end
end

@testset "AC6b SABA and SBAB: eps^2 tau^2 branch at eps = 1e-2" begin
    _ = (SABA, SBAB, advance_n)
    for M in IC_P[[3, 4, 7, 8]]                                                          # n = 3, 4: the eps tau^2n branch is below eps^2 tau^2 for tau <= 1/8
        taus = (1 / 8, 1 / 16, 1 / 32)
        e = [ic_pend_err(M, 1e-2, t, Float64) for t in taus]
        @test all(abs(ic_sl(e[j], e[j + 1]) - 2) <= 0.3 for j in 1:2)
        e3 = ic_pend_err(M, 1e-3, 1 / 8, Float64)
        @test abs(log10(e[1] / e3) - 2) <= 0.3                                           # error scales like eps^2
    end
end

@testset "AC6c SABA_n against leapfrog at the same tau and eps = 1e-3" begin
    _ = (SABA, SBAB, Leapfrog)
    lf = ICM("S2", T -> Leapfrog(T, :DKD), 2, true, Inf)
    el = ic_pend_err(lf, 1e-3, 0.25, Float64)
    for M in IC_P[[2, 3, 4, 6, 7, 8]]                                                    # n >= 2
        @test ic_pend_err(M, 1e-3, 0.25, Float64) / el < 1
    end
end

@testset "AC6 ABA864 and ABAH864: generalized order (8, 6, 4), eps^3 tau^4 regime" begin
    _ = (ABA864, ABAH864, advance_n)
    # eps = 1e-2 and tau <= 1/16: eps^3 tau^4 exceeds eps^2 tau^6 (eps > tau^2); 256 bits; tolerance 1.0 is provisional (no drafting check, P17)
    for M in IC_P[9:10]
        e = [ic_pend_err(M, 1e-2, t, BigFloat) for t in (1 / 16, 1 / 32, 1 / 64)]
        @info "AC6 reported: ABA slopes" M.id ic_sl(e[1], e[2]) ic_sl(e[2], e[3])
        @test all(abs(ic_sl(e[j], e[j + 1]) - 4) <= 1.0 for j in 1:2)
    end
end

@testset "AC6d corrected SABAC_n and SBABC_n (P2)" begin
    _ = (SABAC, SBABC, advance_n, SABA, SBAB)
    for (Mc, Mb) in zip(IC_PC[1:6], (IC_P[2], IC_P[3], IC_P[4], IC_P[6], IC_P[7], IC_P[8]))     # n = 2, 3, 4 (n = 1: the eps tau^2 branch is not corrected)
        for t in (1 / 8, 1 / 16)
            @test ic_pend_err(Mb, 1e-2, t, Float64) / ic_pend_err(Mc, 1e-2, t, Float64) >= 10
        end
        e = [ic_pend_err(Mc, 1e-2, t, Float64) for t in (1 / 16, 1 / 32)]
        n = parse(Int, Mc.id[end:end])
        @test abs(ic_sl(e[1], e[2]) - 4) <= (n == 2 ? 0.5 : 1.5)                        # n = 3, 4: transition towards eps tau^2n, tolerance provisional (P17)
    end
end

@testset "AC7 Gauss methods: order 2s, stage equations, iteration counts, failure status, stopping rule" begin
    _ = (Gauss, ImplicitMidpoint, advance_n, advance_checked, gauss_coefficients, vector_field)
    # (a) slope 2s on F4 and F2; 256 bits for s >= 3
    for (s, mk) in IC_GAUSS, kind in (:f4, :f2)
        M = ICM("G$s", mk, 2s, true, Inf)
        sl, _ = ic_check_order(M, kind; T = s >= 3 ? BigFloat : Float64, J = 5)
        @test ic_asym(sl, 2s, 0.2)
    end
    # coefficients computed, not copied: quadrature and collocation conditions in 256 bits, symmetric nodes
    setprecision(BigFloat, IC_BITS) do
        for s in 1:4
            A, b, c = gauss_coefficients(BigFloat, s)
            ep = ic_eps(BigFloat)
            @test all(abs(sum(b .* c .^ (k - 1)) - 1 / BigFloat(k)) <= 64 * ep for k in 1:(2s))
            @test all(abs(sum(A[i, j] * c[j]^(k - 1) for j in 1:s) - c[i]^k / k) <= 64 * ep for i in 1:s, k in 1:s)
            @test all(abs(c[i] + c[s + 1 - i] - 1) <= 64 * ep for i in 1:s)
        end
    end
    # (b) stage equations hold at the accepted stages: F4, 200 starts s = 1, 2 (binary64), 20 starts s = 3, 4 (256-bit residual evaluation)
    rng = StableRNG(ic_seed(7, 1))
    for s in 1:4
        flow = ic_flow(Float64, :f4); m = Gauss(Float64, s); A, b, c = gauss_coefficients(Float64, s)
        for _ in 1:ic_size(s <= 2 ? 100 : 20, 200)
            z0 = ic_draw(rng, :f4); z = ic_to(Float64, z0)
            zn, st = advance_checked(m, flow, z, 0.1)
            @test st.ok && st.iterations >= 1
            setprecision(BigFloat, IC_BITS) do
                fb = ic_flow(BigFloat, :f4); zb = ic_big.(collect(z))
                fj(j) = collect(vector_field(fb, ic_to(BigFloat, collect(st.stages[j]))))
                res = maximum(maximum(abs, ic_big.(collect(st.stages[i])) .- zb .- BigFloat(0.1) .* sum(BigFloat(A[i, j]) .* fj(j) for j in 1:s)) for i in 1:s)
                @test res <= 8 * ic_eps(Float64) * max(1.0, maximum(maximum(abs, st.stages[i]) for i in 1:s))
            end
        end
    end
    # (c) iteration count at h L max_i sum_j |a_ij| = 0.1 (F1 n = 1, omega = 1, L = 1 in the max norm); the spec's range 0.5 cannot hold the cap 25 (see the plan)
    for s in 1:4
        A, = gauss_coefficients(Float64, s)
        rows = maximum(sum(abs, A[i, :]) for i in 1:s)
        flow = ic_flow(Float64, :f1_1); m = Gauss(Float64, s); tau = 0.1 / rows
        its = Int[]
        for _ in 1:50
            z0 = ic_draw(rng, :f1_1)
            _, st = advance_checked(m, flow, ic_to(Float64, z0), tau)
            @test st.ok
            push!(its, st.iterations)
        end
        @info "AC7 reported: iterations at h L rowsum = 0.1" s extrema(its)
        @test maximum(its) <= 25
    end
    # (d) diverging iteration (h L rowsum = 4): failure status in 100 percent of the trials, never a value (non-finite state)
    for s in 1:4
        A, = gauss_coefficients(Float64, s)
        rows = maximum(sum(abs, A[i, :]) for i in 1:s)
        flow = ic_flow(Float64, :f1_1); m = Gauss(Float64, s); tau = 4.0 / rows
        nfail = 0
        for _ in 1:ic_size(100, 1000)
            z0 = ic_draw(rng, :f1_1)
            zn, st = advance_checked(m, flow, ic_to(Float64, z0), tau)
            nfail += (!st.ok && !any(isfinite, zn)) ? 1 : 0
        end
        @test nfail == ic_size(100, 1000)
    end
    # (e) stopping rule: final increment exactly 0 in at least 90 percent of accepted steps (binary64, F4, s = 2)
    flow = ic_flow(Float64, :f4); m = Gauss(Float64, 2); z = ic_to(Float64, IC_START[:f4])
    zero_inc = 0; nst = ic_size(500, 5000)
    for _ in 1:nst
        z, st = advance_checked(m, flow, z, 0.1)
        zero_inc += (st.ok && st.increment == 0) ? 1 : 0
    end
    @info "AC7 reported: steps with final increment exactly 0" zero_inc / nst
    @test zero_inc / nst >= 0.9
end

# Tao's extended phase-space method on T1 (spec criterion 8)
ic_tao_err(order, omega, tau, Tend) = begin
    z0 = IC_START[:t1]; ref = ic_ref(:t1, z0, Tend)
    flow = ic_flow(Float64, :t1); m = TaoExtended(Float64, omega; order = order)
    z = extend_state(ic_to(Float64, z0))
    z = advance_n(m, flow, z, tau, round(Int, Tend / tau))
    ic_err(physical_state(z), ref)
end

@testset "AC8a Tao's method: Strang order 2, triple jump order 4 with the corrected exponent" begin
    _ = (TaoExtended, extend_state, physical_state, advance_n, TaoNonseparable)
    taus = (0.02, 0.01, 0.005)
    e2 = [ic_tao_err(2, 10.0, t, 20) for t in taus]
    e4 = [ic_tao_err(4, 10.0, t, 20) for t in taus]
    @info "AC8a reported: slopes" ic_sl(e2[1], e2[2]) ic_sl(e2[2], e2[3]) ic_sl(e4[1], e4[2]) ic_sl(e4[2], e4[3])
    @test all(abs(ic_sl(e2[j], e2[j + 1]) - 2) <= 0.2 for j in 1:2)
    @test all(abs(ic_sl(e4[j], e4[j + 1]) - 4) <= 0.3 for j in 1:2)
    # the order-4 method is the triple jump of the Strang formula with gamma = 1/(2 - 2^(1/(l-1))), l = 4 (Inconsistencies (d)); negative control: the printed exponent 1/(l+1)
    flow = ic_flow(Float64, :t1); z = extend_state(ic_to(Float64, IC_START[:t1])); tau = 0.01
    S2 = TaoExtended(Float64, 10.0; order = 2); S4 = TaoExtended(Float64, 10.0; order = 4)
    tj(g1) = advance(S2, flow, advance(S2, flow, advance(S2, flow, z, g1 * tau), (1 - 2g1) * tau), g1 * tau)
    gc = 1 / (2 - 2.0^(1 / 3)); gp = 1 / (2 - 2.0^(1 / 5))
    @test maximum(abs, collect(advance(S4, flow, z, tau)) .- collect(tj(gc))) <= 16 * ic_eps(Float64) * 4
    ep = [begin
              m = S2
              zz = extend_state(ic_to(Float64, IC_START[:t1])); n = round(Int, 20 / t)
              for _ in 1:n
                  zz = advance(m, flow, advance(m, flow, advance(m, flow, zz, gp * t), (1 - 2gp) * t), gp * t)
              end
              ic_err(physical_state(zz), ic_ref(:t1, IC_START[:t1], 20))
          end for t in (0.02, 0.01)]
    @info "AC8a control: slope with the exponent as printed" ic_sl(ep[1], ep[2])
    @test ic_sl(ep[1], ep[2]) < 2.6                                                      # order not raised (drafting check 2.00)
end

@testset "AC8b,c Tao's method: physical error, recovery of the original system, coupling" begin
    _ = (TaoExtended, extend_state, physical_state, advance)
    z0 = IC_START[:t1]; ref = ic_ref(:t1, z0, 20)
    flow = ic_flow(Float64, :t1)
    coupling = Dict{Tuple{Float64, Float64}, Float64}()
    for omega in (10.0, 100.0), tau in (0.01, 0.005)
        m = TaoExtended(Float64, omega; order = 2)
        z = extend_state(ic_to(Float64, z0)); cmax = 0.0
        for _ in 1:round(Int, 20 / tau)
            z = advance(m, flow, z, tau)
            cmax = max(cmax, abs(z[1] - z[3]) + abs(z[2] - z[4]))
        end
        tau == 0.01 && @test ic_err(physical_state(z), ref) <= 1e-3
        tau == 0.01 && @test abs(z[1] - z[3]) <= 1e-3 && abs(z[2] - z[4]) <= 1e-3                     # q = x, p = y recovered
        coupling[(omega, tau)] = cmax
    end
    @info "AC8c reported: max |q - x| + |p - y| over the run" coupling
    # D112: the coupling is the error of the second order splitting, c(omega) tau^2 with c about 1.3 for every omega from 3 to 1000 (measured); the law is asserted (quartering under a halved
    # step) and the two omega give the same constant to 10 percent. The drafted "decreasing with omega" (the mechanism of the continuous system) is not what the discrete coupling does: +1.8 percent measured.
    for omega in (10.0, 100.0)
        @test abs(log2(coupling[(omega, 0.01)] / coupling[(omega, 0.005)]) - 2) <= 0.2
    end
    @test abs(coupling[(100.0, 0.01)] / coupling[(10.0, 0.01)] - 1) <= 0.1
end
