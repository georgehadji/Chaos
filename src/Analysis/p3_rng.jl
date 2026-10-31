# Study p3-abstention (Amendment D93): named streams (E55). One stream per purpose and unit, seeded from the first 8 bytes of SHA-256 over
# "chaoslab-analysis-1|<test_seed>|p3|<purpose>|<unit>"; the generator is xoshiro256++ seeded by SplitMix64; a uniform is in the open unit interval; an integer in 1 to m is
# drawn by rejection; a standard normal variate is the inverse distribution function (AS241) of one uniform. The draws of a unit are a function of its stream only.

"The four state words of xoshiro256++ from a seed: the first four outputs of SplitMix64."
function p3_splitmix_state(seed::UInt64)
    s = seed
    out = ntuple(_ -> UInt64(0), 4)
    for i in 1:4
        s += 0x9e3779b97f4a7c15
        z = s
        z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
        z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
        out = Base.setindex(out, z ⊻ (z >> 31), i)
    end
    return out
end

"The generator xoshiro256++ (its four state words)."
mutable struct P3Gen
    s0::UInt64
    s1::UInt64
    s2::UInt64
    s3::UInt64
end

p3_generator(state::NTuple{4, UInt64}) = P3Gen(state...)

"The next 64-bit output of xoshiro256++."
function p3_u64!(g::P3Gen)
    r = bitrotate(g.s0 + g.s3, 23) + g.s0
    t = g.s1 << 17
    g.s2 ⊻= g.s0; g.s3 ⊻= g.s1; g.s1 ⊻= g.s2; g.s0 ⊻= g.s3
    g.s2 ⊻= t
    g.s3 = bitrotate(g.s3, 45)
    return r
end

_p3_unit(u::Tuple) = join(u, ",")
_p3_unit(u) = string(u)

"The `UInt64` seed of the stream of a unit (E55); a pair unit is written `a,b`. It is the rule of `stream_seed` with the label `p3|<purpose>|<unit>`."
p3_stream_seed(test_seed::Integer, purpose::AbstractString, unit) = stream_seed("p3|$purpose|$(_p3_unit(unit))", test_seed)

"The generator of the stream of a unit."
p3_stream(test_seed::Integer, purpose::AbstractString, unit) = p3_generator(p3_splitmix_state(p3_stream_seed(test_seed, purpose, unit)))

"A uniform `(k + 1/2) / 2^53` with `k` the top 53 bits of an output; it lies in the open unit interval."
p3_uniform!(g::P3Gen) = (Float64(p3_u64!(g) >> 11) + 0.5) / 2.0^53

"An integer in 1 to `m`, unbiased: `x mod m + 1` for the first output `x <= 2^64 - 1 - (2^64 mod m)`."
function p3_randint!(g::P3Gen, m::Integer)
    mm = UInt64(m)
    lim = typemax(UInt64) - (typemax(UInt64) % mm + 1) % mm
    while true
        x = p3_u64!(g)
        x <= lim && return Int(x % mm) + 1
    end
end

# AS241 (PPND16): rational approximations in three regions, relative accuracy of about 1e-16.
const _AS241_A = (3.3871328727963666080, 1.3314166789178437745e2, 1.9715909503065514427e3, 1.3731693765509461125e4, 4.5921953931549871457e4, 6.7265770927008700853e4,
                  3.3430575583588128105e4, 2.5090809287301226727e3)
const _AS241_B = (1.0, 4.2313330701600911252e1, 6.8718700749205790830e2, 5.3941960214247511077e3, 2.1213794301586595867e4, 3.9307895800092710610e4,
                  2.8729085735721942674e4, 5.2264952788528545610e3)
const _AS241_C = (1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550, 3.64784832476320460504, 1.27045825245236838258, 2.41780725177450611770e-1,
                  2.27238449892691845833e-2, 7.74545014278341407640e-4)
const _AS241_D = (1.0, 2.05319162663775882187, 1.67638483018380384940, 6.89767334985100004550e-1, 1.48103976427480074590e-1, 1.51986665636164571966e-2,
                  5.47593808499534494600e-4, 1.05075007164441684324e-9)
const _AS241_E = (6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580, 2.96560571828504891230e-1, 2.65321895265761230930e-2,
                  1.24266094738807843860e-3, 2.71155556874348757815e-5, 2.01033439929228813265e-7)
const _AS241_F = (1.0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2, 7.86869131145613259100e-4, 1.84631831751005468180e-5,
                  1.42151175831644588870e-7, 2.04426310338993978564e-15)

_poly(c, x) = evalpoly(x, c)

"The quantile of the standard normal distribution at `u` in the open unit interval (algorithm AS241); an `ArgumentError` outside it."
function p3_normal_quantile(u::Real)
    (0 < u < 1) || throw(ArgumentError("the normal quantile needs a probability in the open unit interval, got $u"))
    q = u - 0.5
    if abs(q) <= 0.425
        r = 0.180625 - q * q
        return q * _poly(_AS241_A, r) / _poly(_AS241_B, r)
    end
    r = sqrt(-log(q < 0 ? u : 1 - u))
    val = if r <= 5.0
        r -= 1.6
        _poly(_AS241_C, r) / _poly(_AS241_D, r)
    else
        r -= 5.0
        _poly(_AS241_E, r) / _poly(_AS241_F, r)
    end
    return q < 0 ? -val : val
end

"A standard normal variate: the inverse distribution function of one uniform of the stream."
p3_normal!(g::P3Gen) = p3_normal_quantile(p3_uniform!(g))
