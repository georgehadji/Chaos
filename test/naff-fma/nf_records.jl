# Recorded lines (no threshold, never a gate; the rli D77 form): the cheap, deterministic parts of R2 (estimator error laws) and R8 (the N at which E1 of each order reaches its rounding floor), written with
# @info next to the printed orders. R1 (class separation on the 1000-start panel), R3 (observables), R4 (sticky orbits), R5 (cost), R6 (a 2D map with an elliptic point: Y2 does not exist yet), the
# F_1 / F_2 test function of Laskar2003 Eq. (79) in R2(b) and the full R7 / R8 are stage 04 or later (the test plan maps them). R7 (Float64 against Double64) is reported inside AC13. Nothing here is asserted beyond
# finiteness and that the estimator returns.

@testset "R2 estimator error laws (reported)" begin
    _ = (frequency_map, weighted_phase_frequency)
    # E1, window order p in 0..5, the two-line family of the drafting check (nu_1 = 0.1234567891 with amplitude 1, a second line of amplitude 0.1 at a separation of 0.0113), N in {256, 512, 1024, 2048}:
    # error of the first window against nu_1 and the log-log slope, next to the printed theory 2p + 2 (Laskar2003 Eq. (63), p.11; the fitted a_0 of Table 1, p.15, for a different function: 1.98, 3.94, 5.82,
    # 7.64, 8.89, 9.15 for p = 0..5). Float64.
    nu1 = nf_dec("0.1234567891")
    Ns = (256, 512, 1024, 2048)
    tab = Dict{Int, Vector{Float64}}()
    for p in 0:5
        errs = Float64[]
        for N in Ns
            z = nf_series(Float64, nf_lines_big([(nu1, BigFloat(1), BigFloat(0)), (nu1 + nf_dec("0.0113"), nf_dec("0.1"), BigFloat(0))], 2N))
            r = frequency_map(z, N; order = p)
            push!(errs, nf_cderr(r.w1[1].nu, nu1))
        end
        tab[p] = errs
        slope = -sum((log10.(Ns) .- sum(log10.(Ns)) / 4) .* (log10.(max.(errs, 1e-17)) .- sum(log10.(max.(errs, 1e-17))) / 4)) / sum((log10.(Ns) .- sum(log10.(Ns)) / 4) .^ 2)
        @info "R2(a) E1 two-line error vs N, order p=$p (theory 2p+2 = $(2p + 2))" errors = string(round.(errs; sigdigits = 3)) slope = round(slope; digits = 2)
        @test all(isfinite, errs)
    end
    # E2 on NF15 with the exponential, sin^2 and equal weights (Das2017 p.11: faster than any power / N^-3 / N^-1 as printed), N in {250, 500, 1000, 2000}
    incr, rho = nf_nf15_big(4000)
    for (name, wf) in (("exponential", N -> nothing), ("sin^2", N -> (w = [sin(Float64(pi) * j / N)^2 for j in 0:(N - 1)]; w ./ sum(w))), ("equal", N -> fill(1 / N, N)))
        errs = Float64[]
        for N in (250, 500, 1000, 2000)
            inT = Float64[Float64(x) for x in incr[1:(2N)]]
            w = wf(N)
            nu = w === nothing ? weighted_phase_frequency(inT, N) : weighted_phase_frequency(inT, N; weights = w)
            push!(errs, Float64(abs(nf_big(nu[1]) - rho)))
        end
        @info "R2 E2 on NF15, $name weight, error vs N in (250, 500, 1000, 2000)" errors = string(round.(errs; sigdigits = 3))
        @test all(isfinite, errs)
    end
end
