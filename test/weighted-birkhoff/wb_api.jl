# API tests: the pure function weighted_birkhoff, the index birkhoff_index, the observable kinds, the reducer constructor and provenance, and the reducer protocol on hand arithmetic. Every testset
# starts with a line naming the new symbols it uses (UndefVarError when missing). The names are proposals of the plan (stage 03 may rename them; the spec is WHAT-only).

@testset "API pure function weighted_birkhoff" begin
    _ = (weighted_birkhoff,)
    # weighted_birkhoff(series, window; order = 1, weights = nothing) -> (avg1, avg2, S, B): the series is a vector (one observable) or an n x K matrix (a column per observable), `window` is N or a
    # tuple of N (horizons); avg1, avg2, S are K x H matrices of the element type, B a vector of H; the pair is [0, N), [N, 2N) of the first 2N samples; S the largest |f| of the 2N samples (the
    # samples of weight zero count). Hand values: f_j = j gives N/2 and 3N/2 (the weighted mean of jj is N/2 by symmetry); equal weights give (N - 1)/2 and (3N - 1)/2.
    for T in WB_TYPES
        epsT = wb_eps(T)
        f = T.(0:31)
        r = weighted_birkhoff(f, 16)
        @test r.avg1 isa Matrix{T} && r.avg2 isa Matrix{T} && r.S isa Matrix{T} && r.B isa Vector{T} && size(r.avg1) == (1, 1) && length(r.B) == 1
        @test abs(r.avg1[1] - 8) <= 4 * epsT * 31 && abs(r.avg2[1] - 24) <= 4 * epsT * 31 && r.S[1] == 31
        @test wb_err(r.B[1], sum(wb_wref(16, 1))) <= 4 * epsT * Float64(sum(wb_wref(16, 1)))
        @test isequal(r, weighted_birkhoff(f, 16; order = 1)) && isequal(r, weighted_birkhoff(reshape(f, :, 1), 16))        # defaults, vector = n x 1 matrix
        r2 = weighted_birkhoff(f, 16; order = 2)
        @test abs(r2.avg1[1] - 8) <= 16 * epsT * 31 && abs(r2.avg2[1] - 24) <= 16 * epsT * 31 && r2.B[1] != r.B[1]
        # K = 2 columns and H = 2 horizons: shapes K x H, horizon values (N = 32: the means of jj are 16 and 48), horizon 1 equal to the single-horizon call
        two = [T.(0:63) T(100) .- T.(0:63)]
        rr = weighted_birkhoff(two, (16, 32))
        @test size(rr.avg1) == (2, 2) && length(rr.B) == 2
        @test all(abs.(rr.avg1 .- T[8 16; 92 84]) .<= 4 * epsT * 100) && all(abs.(rr.avg2 .- T[24 48; 76 52]) .<= 4 * epsT * 100) && rr.S == T[31 63; 100 100]
        @test isequal(rr.avg1[:, 1], weighted_birkhoff(two, 16).avg1[:, 1]) && isequal(rr.B[1], weighted_birkhoff(two, 16).B[1])
        # samples beyond 2N are ignored
        longer = vcat(f, T[1e6, -1e6])
        @test isequal(weighted_birkhoff(longer, 16), r)
        # equal weights (a weight vector of length N, used as given)
        re = weighted_birkhoff(f, 16; weights = ones(T, 16))
        @test abs(re.avg1[1] - T(7.5)) <= 4 * epsT * 31 && abs(re.avg2[1] - T(23.5)) <= 4 * epsT * 31 && re.B[1] == 16 && re.S[1] == 31
        # the samples of weight zero count for S but not for the averages: a spike at the first sample of each window
        sp = ones(T, 32)
        sp[1], sp[17] = 100, -200
        rs = weighted_birkhoff(sp, 16)
        @test rs.S[1] == 200 && abs(rs.avg1[1] - 1) <= 4 * epsT * 200 && abs(rs.avg2[1] - 1) <= 4 * epsT * 200
        # a constant series returns the constant; BigFloat is a supported element type at the ambient precision
        rb = setprecision(BigFloat, WB_BITS) do
            weighted_birkhoff(fill(BigFloat(3) / 10, 40), 16; order = 2)
        end
        @test rb.avg1 isa Matrix{BigFloat} && abs(rb.avg1[1] - BigFloat(3) / 10) < 1e-70 && rb.avg1[1] == rb.avg2[1]
        # argument errors naming the parameter (the window, order, series are also in AC9)
        named(g, word) = (e = wb_error(g); e isa ArgumentError && occursin(word, e.msg))
        @test named(() -> weighted_birkhoff(f, 16; weights = ones(T, 15)), "weights") && named(() -> weighted_birkhoff(T.(0:127), (16, 32); weights = ones(T, 16)), "weights")
    end
end

@testset "API index birkhoff_index" begin
    _ = (birkhoff_index,)
    # birkhoff_index(avg1, avg2, S; floor_multiplier = 32) -> (delta, absdig, reldig, dig, below_floor, zero), the type of the arguments (WB15): Delta = avg1 - avg2 unclipped; phi = 32 eps S;
    # absdig = -log10 max(|Delta|, phi); reldig = -log10(max(|Delta|, phi)/m), m = (|avg1| + |avg2|)/2, undefined (NaN) where m <= phi; dig = max of the two where reldig is defined else absdig;
    # below_floor iff |Delta| <= phi; S = 0: digits NaN and zero set.
    for T in WB_TYPES
        epsT = wb_eps(T)
        bad = String[]
        one_, tiny = T(1), T(ldexp(1.0, -10))
        cases = ((one_, one_ - tiny, one_), (one_ - tiny, one_, one_), (T(3), T(3) + T(0.25), T(4)), (T(1e-3), T(2e-3), T(1)),
                 (T(epsT), T(2epsT), T(1)), (one_, one_, one_), (T(1e-20), T(1e-20), T(1e-20)), (T(5), T(-5), T(5)), (one_, one_ - T(ldexp(1.0, -40)), T(1e3)))
        for (i, (a1, a2, S)) in enumerate(cases)
            ix = birkhoff_index(a1, a2, S)
            ix.delta isa T && ix.absdig isa T && ix.reldig isa T && ix.dig isa T && ix.below_floor isa Bool && ix.zero isa Bool || push!(bad, "case $i types")
            isequal(ix.delta, a1 - a2) || push!(bad, "case $i delta")
            wb_check_index!(bad, "T=$T case $i", T, a1, a2, S, ix.delta, ix)
        end
        # the first case by hand: Delta = 2^-10, absdig = 10 log10 2, m = 1 - 2^-11
        ix = birkhoff_index(one_, one_ - tiny, one_)
        @test ix.delta == tiny && abs(ix.absdig - 10 * log10(T(2))) <= 4 * epsT * 10 && !ix.below_floor && !ix.zero
        # the floor: |Delta| = phi exactly is flagged, one unit more is not; phi scales with S (here S = 1)
        phi = T(32) * effective_eps(T)
        @test birkhoff_index(one_, one_ - phi, one_).below_floor && !birkhoff_index(one_, one_ - T(2) * phi, one_).below_floor
        fl = birkhoff_index(one_, one_, one_)
        @test fl.delta == 0 && fl.below_floor && abs(fl.absdig + log10(phi)) <= 4 * epsT * 15 && fl.reldig == fl.absdig && fl.dig == fl.absdig
        # m <= phi: reldig undefined, dig = absdig
        un = birkhoff_index(T(epsT), T(2epsT), one_)
        @test isnan(un.reldig) && un.dig == un.absdig && un.below_floor
        # signed Delta, symmetric digits
        @test birkhoff_index(T(2), T(3), T(3)).delta == -1 && isequal(birkhoff_index(T(2), T(3), T(3)).dig, birkhoff_index(T(3), T(2), T(3)).dig)
        # S = 0
        z = birkhoff_index(zero(T), zero(T), zero(T))
        @test z.zero && isnan(z.absdig) && isnan(z.reldig) && isnan(z.dig) && z.delta == 0
        # the keyword and its check
        ix8 = birkhoff_index(one_, one_, one_; floor_multiplier = 8)
        @test abs(ix8.absdig + log10(T(8) * effective_eps(T))) <= 4 * epsT * 16
        e = wb_error(() -> birkhoff_index(one_, one_, one_; floor_multiplier = 0))
        @test e isa ArgumentError && occursin("floor_multiplier", e.msg)
        # broadcast over a matrix of entries (the form the tests use on the K x H results)
        @test size(birkhoff_index.(fill(one_, 2, 3), fill(one_, 2, 3), fill(one_, 2, 3))) == (2, 3)
        @test wb_report(bad, "API index $T")
    end
end

@testset "API observable kinds" begin
    _ = (CosObservable, SinObservable, LiftedObservable, observable_name)
    # CosObservable(i, k = 1)(x) = cos(k x[i]), SinObservable(i, k = 1)(x) = sin(k x[i]) (on the wrapped coordinate; k an integer frequency, exact on the circle), LiftedObservable(iq, im)(x) = x[iq] +
    # 2 pi_T x[im], the lifted coordinate rebuilt from the winding count, bitwise the `lift` of the wrapped map (D28); isbits; observable_name gives the identifier recorded in the provenance.
    for T in WB_TYPES
        epsT = wb_eps(T)
        x = SVector{4, T}(T(0.25), T(0.5), T(2), T(3))
        @test abs(CosObservable(1)(x) - cos(T(0.25))) <= 2 * epsT && abs(CosObservable(2, 3)(x) - cos(3 * T(0.5))) <= 2 * epsT && abs(SinObservable(1, 2)(x) - sin(2 * T(0.25))) <= 2 * epsT
        @test CosObservable(1) == CosObservable(1, 1) && SinObservable(2) == SinObservable(2, 1)
        @test LiftedObservable(2, 4)(x) === lift(x)[2] && LiftedObservable(1, 3)(x) === lift(x)[1]
        @test CosObservable(1)(x) isa T && LiftedObservable(2, 4)(x) isa T
    end
    @test isbitstype(CosObservable) && isbitstype(SinObservable) && isbitstype(LiftedObservable)
    n1, n2, n3 = observable_name(CosObservable(1, 2)), observable_name(SinObservable(1, 2)), observable_name(LiftedObservable(2, 4))
    @test n1 isa String && n2 isa String && n3 isa String && occursin("cos", n1) && occursin("sin", n2) && occursin("lift", n3) && n1 != observable_name(CosObservable(1, 1)) && n1 != n2
    @test observable_name(wb_Coord(1)) isa String
end

@testset "API reducer constructor and provenance" begin
    _ = (WeightedBirkhoff, CosObservable, SinObservable, LiftedObservable, reducer_provenance, reducer_workspace, reducer_bytes)
    # WeightedBirkhoff(; window, observables, birkhoff_order = 1): window is N or a tuple of N (horizons of one run), observables a tuple of K isbits callables of the orbit point (any number type is
    # converted to the element type); the reducer is isbits, a reader (any k >= 1), with no workspace; the provenance is merged into the propagate record.
    obs = (CosObservable(1), SinObservable(1, 2), LiftedObservable(2, 4))
    r = WeightedBirkhoff(; window = 100, observables = obs)
    @test isbitstype(typeof(r))
    pr = reducer_provenance(r)
    @test pr.indicator == "weighted-birkhoff" && pr.weight == "wb_exp" && pr.birkhoff_order == 1 && pr.window == (100,) && pr.handoffs == 200 && pr.floor_multiplier == 32 && pr.min_window == 16
    @test pr.observables isa Tuple && length(pr.observables) == 3 && all(s -> s isa String, pr.observables) && occursin("cos", pr.observables[1]) && occursin("sin", pr.observables[2])
    pl = reducer_provenance(WeightedBirkhoff(; window = (50, 100), observables = (wb_Coord(1),), birkhoff_order = 2))
    @test pl.window == (50, 100) && pl.birkhoff_order == 2 && pl.handoffs == 200 && length(pl.observables) == 1
    for T in WB_TYPES, k in (1, 2, 3)
        @test reducer_workspace(r, T, Val(2), Val(k)) == 0
        st = reducer_init(r, T, Val(2), Val(k))
        @test isbitstype(typeof(st)) && !reducer_failed(r, st)
        f0 = reducer_final(r, st)
        @test f0.handoffs == 0 && f0.n == 0 && !f0.complete && all(isnan, f0.avg1) && all(isnan, f0.avg2) && all(isnan, f0.delta) && all(isnan, f0.dig) && size(f0.avg1) == (3, 1)
    end
end

@testset "API reducer protocol on hand arithmetic" begin
    _ = (WeightedBirkhoff, reducer_step)
    # N = 16 (the minimum window), K = 2 observables (the two coordinates of the point), samples x_j = (j, 100 - j), j = 0..: the first handoff is the sample 0 (weight 0 but counted in S); avg_1 appears
    # at the handoff N (8 and 92 by symmetry), avg_2, Delta and the digits at 2N (24 and 76, Delta = -16 and 16); the running S is 15 and 100 at the handoff 16 and 31 and 100 at 32; W is returned
    # unchanged by every step; after the horizon is complete further handoffs change nothing; two horizons share the run.
    N = 16
    for T in WB_TYPES
        epsT = wb_eps(T)
        r = WeightedBirkhoff(; window = N, observables = (wb_Coord(1), wb_Coord(2)))
        W = SMatrix{2, 1, T}(1, 0)
        xs = [SVector{2, T}(T(j), T(100 - j)) for j in 0:(4N)]
        st = reducer_init(r, T, Val(2), Val(1))
        for n in 1:(2N)
            st, W2 = reducer_step(r, st, n, xs[n], W)
            @test W2 === W
            f = reducer_final(r, st)
            @test f.handoffs == n && f.n == n && f.complete == (n == 2N)
            n < N && @test all(isnan, f.avg1)
            n < 2N && @test all(isnan, f.avg2) && all(isnan, f.delta) && all(isnan, f.absdig)
            n == 1 && @test f.S[1] == 0 && f.S[2] == 100                                   # the sample 0 counts for S
            n == N && @test abs(f.avg1[1] - 8) <= 4 * epsT * 31 && abs(f.avg1[2] - 92) <= 4 * epsT * 100 && f.S[1] == 15 && f.S[2] == 100
        end
        f = reducer_final(r, st)
        @test !reducer_failed(r, st) && f.complete && size(f.avg1) == (2, 1) && size(f.B) == (1,)
        @test abs(f.avg2[1] - 24) <= 4 * epsT * 31 && abs(f.avg2[2] - 76) <= 4 * epsT * 100 && f.S[1] == 31 && f.S[2] == 100
        @test abs(f.delta[1] + 16) <= 8 * epsT * 31 && abs(f.delta[2] - 16) <= 8 * epsT * 100 && f.delta[1] == f.avg1[1] - f.avg2[1]
        @test wb_err(f.absdig[1], -log10(BigFloat(16))) <= 16 * epsT && !f.below_floor[1] && !f.zero[1]
        g = weighted_birkhoff(hcat(T.(0:(2N - 1)), T(100) .- T.(0:(2N - 1))), N)
        @test isequal(Matrix(f.avg1), g.avg1) && isequal(Matrix(f.avg2), g.avg2) && isequal(Matrix(f.S), g.S) && isequal(Vector(f.B), g.B)
        # one more handoff after the horizon is complete changes nothing but the counters
        st2, _ = reducer_step(r, st, 2N + 1, xs[2N + 1], W)
        f2 = reducer_final(r, st2)
        @test f2.handoffs == 2N + 1 && all(fld -> isequal(getfield(f2, fld), getfield(f, fld)), (:avg1, :avg2, :S, :B, :delta, :absdig, :reldig, :dig, :below_floor, :zero))
        # two horizons (16, 32) share the run of 64 handoffs: the first completes at 32 and equals the one-horizon reducer, `complete` waits for the largest
        r2 = WeightedBirkhoff(; window = (N, 2N), observables = (wb_Coord(1), wb_Coord(2)))
        s2 = reducer_init(r2, T, Val(2), Val(1))
        for n in 1:(4N)
            s2, _ = reducer_step(r2, s2, n, xs[n], W)
            n == 2N && @test !reducer_final(r2, s2).complete && all(!isnan, reducer_final(r2, s2).delta[:, 1]) && all(isnan, reducer_final(r2, s2).delta[:, 2])
        end
        h = reducer_final(r2, s2)
        @test h.complete && size(h.avg1) == (2, 2) && isequal(h.avg1[:, 1], f.avg1[:, 1]) && isequal(h.avg2[:, 1], f.avg2[:, 1]) && isequal(h.B[1], f.B[1])
        @test abs(h.avg1[1, 2] - 16) <= 4 * epsT * 63 && abs(h.avg2[1, 2] - 48) <= 4 * epsT * 63
    end
end
