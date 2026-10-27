# Recorded lines (no threshold, K18): the loose sanity values of the spec's Oracle section (R8), reported next to the drafting values. The class separation (R1), the sensitivity (R2), the
# weak-chaos test (R3, deferred to the coupled maps), the L13 distributions (R4), the realized ratios (R5, reported by the @info of AC6, AC13), the cost (R6), the printed properties (R7)
# are stage 04 or later (the test plan maps them). Nothing here is asserted beyond finiteness, the range rule and the conservation of the counts.

@testset "R8 sanity values of the drafting check (reported)" begin
    _ = (StretchingNumbers, stretching_range)
    # Standard map, Float64, wrapped, orthonormal start pair (0.6, 0.8) and (-0.8, 0.6), B = 1000 over the default range, cumulative counts, one start each (not the L1 S7 panel), n = 1e5:
    # D at Z = 1e3, 1e4, 1e5 and D times Z, the mean of z. The spec's drafting values (not sources, not citable): K = 0.5 start (2, 0): D 0.878, 0.100, 0.0100; start (3, 0): 1.06, 0.257, 0.102;
    # K = 2 pi start (3, 0.7): 0.0363, 0.00363, 0.00036; K = 5 start (1, 0.5): 0.0311, 0.00311, 0.00031. Reported, never pass / fail.
    T = Float64
    cases = ((0.5, (2.0, 0.0)), (0.5, (3.0, 0.0)), (2pi, (3.0, 0.7)), (5.0, (1.0, 0.5)))
    W0 = SMatrix{2, 2, T}(0.6, 0.8, -0.8, 0.6)
    for (K, (x, p)) in cases
        sys = WrappedStandardMap(K)
        r = sn_red(sys, 1000, 1)
        f = sn_run((r,), sys, [torus_state(SVector{2, T}(x, p))], [W0], 100_000, 1; steps = [1000, 10_000, 100_000])[1][1].snaps
        @test all(g -> isfinite(g.D) && g.D >= 0 && g.below == [0, 0] && g.above == [0, 0] && sum(sn_cnt(g, 1)) == g.handoffs, f)
        Ds, DZs = [g.D for g in f], [g.D * g.handoffs for g in f]
        @info "R8 sanity (K, start): D at Z = 1e3, 1e4, 1e5; D times Z; mean z at 1e5, column 1" K x p D = string(round.(Ds; sigdigits = 4)) DZ = string(round.(DZs; sigdigits = 4)) mean_z = round(f[end].mean[1]; sigdigits = 4)
    end
end
