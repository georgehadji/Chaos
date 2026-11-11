# AC16 to AC18: execution invariance and allocation, the external reference (package stand-in), the CPU part of the CPU-GPU comparison.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).

@testset "AC16a execution invariance: chunk sizes and the serial body give bitwise identical series" begin
    _ = (OFLI,)
    flow = HenonHeilesFlow(); im = of_im(Float64); n = of_size(48, 4096); N = of_size(1000, 20_000); seed = of_seed(16, 1)
    a, b = section_interval(flow, 1 / 6); w = sqrt(2 / 6); rng = StableRNG(of_seed(16, 2))
    xs = [random_energy_start(flow, 1 / 6, (a, b), (-w, w), rng)[1] for _ in 1:n]
    @info "AC16 reported" threads = Threads.nthreads() starts = n
    for tau in (1, 20)
        serial = [propagate_orbit(im, xs[j], initial_deviations(seed, j, Val(4), Val(4), Float64), N, tau, (OFLI(flow),); orbit = j) for j in 1:n]
        for chunk in (of_size(7, 64), of_size(64, 1000), nothing)
            res = propagate(im, xs, N; k = 4, tau = tau, seed = seed, reducers = (OFLI(flow),), chunk = chunk)
            @test all(j -> isequal(res.finals[j][1], reducer_final(OFLI(flow), serial[j][3][1])) && res.x[j] == serial[j][1] && res.W[j] == serial[j][2], 1:n)
        end
    end
end

@testset "AC16b the hot loop allocates nothing after warm-up" begin
    _ = (OFLI,)
    flow = HenonHeilesFlow(); z0 = of_orbit_C()
    for T in (Float64, Double64), tau in (1, 20)
        im = of_im(T); x0 = SVector{4, T}(map(T, z0)); W0 = of_frame(of_seed(16, 3), 1, 4, 4, T)
        @test of_alloc_bytes(im, x0, W0, 2000, tau, (OFLI(flow),)) == 0
        W2 = SMatrix{4, 2, T}(W0[:, 1:2])
        @test of_alloc_bytes(im, x0, W2, 2000, tau, (OFLI(flow; replace = false), SALI())) == 0      # the chain with SALI needs k = 2
    end
end

# ------------------------------------------------------------------ AC17: external reference at the package level (Double64 order-10 table at h/4 stands in for TaylorIntegration)
# G_B(t) of the exact flow at the unit-time handoffs t = 0 .. Tmax (the s = 0 term included): Double64 order-10 table at 1/80 with its tangent.
function of_ref_hh(z0, W0, Tmax::Integer)
    D = Double64; m = SofroniouSpaletta10(D); flow = HenonHeilesFlow(); h = D(1) / 80
    z = SVector{4, D}(map(D, z0)); W = SMatrix{4, 4, D}(map(D, W0)); G = -Inf; out = zeros(Tmax + 1)
    for t in 0:Tmax
        f = vector_field(flow, z); ff = dot(f, f)
        for i in 1:4
            w = W[:, i]; b = norm(w - (dot(f, w) / ff) * f); G = max(G, Float64(log(b)))
        end
        out[t + 1] = G
        if t < Tmax
            for _ in 1:80
                z, W = tangent_advance(m, flow, z, W, h)
            end
        end
    end
    return out
end

@testset "AC17 external reference, package stand-in: order-8 table at h = 1/20 against the Double64 order-10 reference (the Taylor form is stage 04)" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); Tmax = 60
    for (name, z0) in ((:R, of_orbit_R()), (:C, of_orbit_C()))
        W0 = of_frame(of_seed(17, name === :R ? 1 : 2), 1, 4, 4, Float64); ref = of_ref_hh(z0, W0, Tmax)
        for key in (:o8, :o4)
            im = of_im(Float64, key); snaps = of_snaps(OFLI(flow), im, z0, W0, 20 * Tmax, 20); g0 = ofli_initial(flow, z0, W0)
            errs = [abs(max(maximum(g0), snaps[t].ofli) - ref[t + 1]) for t in (10, 30, 60)]
            if key === :o8
                @test all(<=(1.0e-6), errs)                                              # delta_8: calibrated at stage 02 (spec: order 1e-6)
            else
                @info "AC17 reported: order-4 table, error at t = 10, 30, 60" orbit = name errors = errs
            end
        end
    end
end

# ------------------------------------------------------------------ AC18: CPU part of the CPU-GPU comparison (sali Amendment D68 form); the device runs are stage 04
@testset "AC18 CPU part: each backend against its own reference, the gap rule, the fallback" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); T100 = of_size(100, 400); N = 20 * T100; m = min(OF_NCLASS, of_size(2, 1024))
    nbad = 0; ngap = 0; nagree = 0
    for ks in ([(:reg, i) for i in 1:m], [(:cha, i) for i in 1:m])
        cs = Dict(T => [of_case(k, T, :o8, N) for k in ks] for T in (Float64, Double64))
        fin = Dict{Any, Any}(); rowsT = Dict{Any, Any}()
        for T in (Float64, Double64)
            c = cs[T]
            res = propagate(c[1].im, [x.x0 for x in c], N; k = 4, tau = 20, seed = 1, reducers = (OFLI(flow),), W0 = [x.W0 for x in c], backend = :cpu)
            @test res.provenance.backend === :cpu
            fin[T] = res
            rowsT[T] = [of_rows(c[j], T, 20, [res.finals[j][1]], ofli_initial(flow, c[j].x0, c[j].W0), [N]) for j in eachindex(c)]
            for rows in rowsT[T]
                all(r -> r.err <= r.bound, rows) || (nbad += 1)
            end
        end
        # gap only where both orbits agree (final states within 1e-9 in every coordinate; the B_G gap rule gates the regular class only): |G_a - G_b| <= B_a + B_b, both forms
        for j in 1:m
            xa, xb = fin[Float64].x[j], fin[Double64].x[j]
            maximum(abs.(Float64.(xb) .- Float64.(xa))) <= 1.0e-9 || continue
            nagree += 1
            ks[j][1] === :reg || continue          # chaotic class: the gap includes the orbit divergence term D_orb (D68), reported here, gated at stage 04
            for (form, field) in ((1, :ofli), (2, :ofli1))
                ra, rb = rowsT[Float64][j][form], rowsT[Double64][j][form]
                ga = form == 1 ? max(maximum(ofli_initial(flow, cs[Float64][j].x0, cs[Float64][j].W0)), Float64(fin[Float64].finals[j][1].ofli)) :
                                 max(ofli_initial(flow, cs[Float64][j].x0, cs[Float64][j].W0)[1], Float64(fin[Float64].finals[j][1].ofli1))
                gb = form == 1 ? max(maximum(ofli_initial(flow, cs[Double64][j].x0, cs[Double64][j].W0)), Float64(fin[Double64].finals[j][1].ofli)) :
                                 max(ofli_initial(flow, cs[Double64][j].x0, cs[Double64][j].W0)[1], Float64(fin[Double64].finals[j][1].ofli1))
                abs(ga - gb) <= ra.bound + rb.bound || (ngap += 1)
            end
        end
    end
    @test nbad == 0 && ngap == 0
    @info "AC18 reported" orbits_in_agreement = nagree
    # a request for the device with no device installed falls back to the CPU and records the reason; the result is the CPU one
    if Chaoslab.Foundation.DEVICE[] isa Chaoslab.Foundation.NoDevice
        c = of_case((:reg, 1), Float64, :o8, N)
        r1 = propagate(c.im, [c.x0], N; k = 4, tau = 20, seed = 1, reducers = (OFLI(flow),), W0 = [c.W0], backend = :gpu)
        r0 = propagate(c.im, [c.x0], N; k = 4, tau = 20, seed = 1, reducers = (OFLI(flow),), W0 = [c.W0], backend = :cpu)
        @test r1.provenance.backend === :cpu && !isempty(string(r1.provenance.backend_fallback_reason)) && isequal(r1.finals, r0.finals)
    end
end
