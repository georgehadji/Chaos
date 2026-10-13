# Criterion 14 (oracle O2): Float64 and Double64 orbits vs 256-bit orbits from identical exactly representable data.

# Largest relative error of sin and cos in Double64 on the sample xs, against 256-bit values (call inside setprecision).
function sm_dd_trig_error(xs)
    m = 0.0
    for x in xs
        b = sm_big(x)
        for (f, g) in ((sin, sin), (cos, cos))
            r = g(b)
            iszero(r) && continue
            m = max(m, Float64(abs(sm_big(f(x)) - r) / abs(r)))
        end
    end
    return m
end

# worst (error / tau) over cases; lifted and torus errors vs the BigFloat orbit, eps = epsT (Float64 or the effective DD value)
function sm_ac14_ratio(::Type{T}, cases, epsT) where {T}
    worst = 0.0
    for (K, z) in cases
        nH = sm_nH(K)
        ref = sm_orbit(StandardMap(BigFloat(K)), sm_cast(BigFloat, z), nH)
        low = sm_orbit(StandardMap(T(K)), sm_cast(T, z), nH)
        tau = sm_tau(epsT, K, ref, nH)
        worst = max(worst, sm_absdiff(low[end], ref[end]) / tau, sm_tdist(low[end], ref[end]) / tau)
    end
    return worst
end

@testset "AC14 element-type agreement" begin
    # |orbit_T - orbit_256| <= tau(n_H) with eps of the lower precision; 200 orbits per listed K; 256 bits via setprecision(do)
    cases = sm_cases(StableRNG(SM_SEED + 14001), 200; uniform = false)
    setprecision(BigFloat, SM_BITS) do
        @test sm_ac14_ratio(Float64, cases, sm_eps(Float64)) <= 1
        # sanity: the Chaoslab 256-bit orbit agrees with an independent plain-formula 256-bit orbit
        worst = 0.0
        for (K, z) in cases[1:400]
            nH = sm_nH(K)
            a = sm_orbit(StandardMap(BigFloat(K)), sm_cast(BigFloat, z), nH)
            b = sm_ref_orbit(BigFloat(K), sm_cast(BigFloat, z), nH)
            worst = max(worst, sm_absdiff(a[end], b[end]) / sm_tau(ldexp(1.0, 1 - SM_BITS), K, a, nH))
        end
        @test worst <= 1
        # Double64: effective eps = largest relative sin/cos error on the sample x values, floored at 2^-104
        xs = Double64[]
        for (K, z) in cases
            append!(xs, first.(sm_orbit(StandardMap(Double64(K)), sm_cast(Double64, z), sm_nH(K))))
        end
        eff = max(SM_EPS_DD_FLOOR, sm_dd_trig_error(xs))
        @test sm_ac14_ratio(Double64, cases, eff) <= 1
        @test eff <= ldexp(1.0, -100)          # otherwise report to the owner (spec criterion 14)
    end
end
