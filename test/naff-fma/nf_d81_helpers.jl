# Helpers for nf_d81.jl (Amendment D81: D81a, D81b and two clauses of criterion 11(b)). Definitions only; none of them calls a name of the module under test.
# Written from the statements of the amendment (NF27, NF30) and checked against the held text data/literature/naff-fma/CERN-SL-95-84_Bartolini1996.txt (Eq. (24), p.9; Eqs. (6), (7), (A.3), (A.5), pp.5, 6, 17).

# Bartolini1996 Eq. (24): the modulus of the transform (1/N) sum_{n=0}^{N-1} 2 sin^2(pi n/N) exp(2 pi i nu_0 n) exp(-2 pi i nu_j n) at a frequency nu_j, with delta = nu_0 - nu_j:
# sin^2(pi/N) |sin(pi N delta) cos(pi delta)| / (N |sin(pi delta) sin(pi delta + pi/N) sin(pi delta - pi/N)|). BigFloat.
function nf_bart_eq24(N::Integer, delta::BigFloat)
    return setprecision(BigFloat, NF_BITS) do
        a = nf_pi() * delta
        h = nf_pi() / N
        sin(h)^2 * abs(sin(N * a) * cos(a)) / (N * abs(sin(a) * sin(a + h) * sin(a - h)))
    end
end

# The two-sine signal of D81a (Bartolini1996 Eq. (A.1)): z_n = exp(2 pi i nu_0 n) + b_0 exp(2 pi i nu_1 n), nu_1 = nu_0 - Delta, n = 0..2N (2N + 1 samples, 2N increments), BigFloat.
nf_d81_two_sine_big(nu0::BigFloat, b0::BigFloat, Delta::BigFloat, N::Integer) = setprecision(BigFloat, NF_BITS) do
    nf_lines_big([(nu0, BigFloat(1), BigFloat(0)), (nu0 - Delta, b0, BigFloat(0))], 2N + 1)
end

# The printed bound with N increments (Eq. (7) with C_APA <= 1/2 and the arctangent bound (A.5), N - 1 read as N increments), plus 16 eps of rounding.
nf_d81a_bound(N::Integer, epsT::Real) = 1 / (2N) + 16 * epsT
