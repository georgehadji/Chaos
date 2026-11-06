# Near-integrable rows P1 to P3 (H = A + eps B with A the free drift): SABA_n and SBAB_n (LaskarRobutel2001 Table I), ABA864 and ABAH864 (Blanes2013 Tables 3 and 4),
# and the corrected SABAC_n and SBABC_n (LaskarRobutel2001 Eq. (62), Table II).
#
# SABA_n and SBAB_n are quadrature formulas (LaskarRobutel2001 Sect. 7, Eq. (47): "the Gauss integration formula"): SABA_n kicks at the n Gauss-Legendre nodes
# with the weights, SBAB_n at the n + 1 Gauss-Lobatto nodes (both on [0, 1], weights summing to 1). They are computed in 320-bit arithmetic and rounded into the element
# type, not copied; the printed decimals of Table I (n = 6 to 10) and the closed forms for the small n agree with them to the 36 digits printed (checked, stage 03 notes).

# Newton iteration for the roots of Legendre P_n (Gauss) or of P_n' (Lobatto, interior nodes); returns nodes and weights on [-1, 1].
function _legendre_nodes(n::Integer, lobatto::Bool)
    function PP(x)                                                  # P_n, P_n'
        a = one(x); b = x
        for k in 2:n
            a, b = b, ((2k - 1) * x * b - (k - 1) * a) / k
        end
        return b, n * (a - x * b) / (1 - x * x)                    # P_n' = n (P_(n-1) - x P_n)/(1 - x^2)
    end
    if lobatto
        xs = [-one(BigFloat)]
        for i in (n - 1):-1:1
            x = cospi(BigFloat(i) / n)
            for _ in 1:60
                p, dp = PP(x); ddp = (2x * dp - n * (n + 1) * p) / (1 - x * x)
                dx = dp / ddp; x -= dx
                abs(dx) < ldexp(BigFloat(1), -300) && break
            end
            push!(xs, x)
        end
        push!(xs, one(BigFloat)); xs = sort(xs)
        w = [2 / (n * (n + 1) * PP(x)[1]^2) for x in xs]
        return xs, w
    end
    xs = BigFloat[]
    for i in 1:n
        x = cospi((i - BigFloat(1) / 4) / (n + BigFloat(1) / 2))
        for _ in 1:60
            p, dp = PP(x); dx = p / dp; x -= dx
            abs(dx) < ldexp(BigFloat(1), -300) && break
        end
        push!(xs, x)
    end
    xs = sort(xs)
    return xs, [2 / ((1 - x * x) * PP(x)[2]^2) for x in xs]
end

# The table of SABA_n (kicks at the n Gauss nodes) or SBAB_n (kicks at the n + 1 Lobatto nodes); stages (c_i, d_i) with the first half mirrored (bitwise symmetry).
function _quadrature_table(::Type{T}, n::Integer, lobatto::Bool) where {T}
    n >= 1 || throw(ArgumentError("SABA/SBAB need n >= 1, got $n"))
    prec = precision(BigFloat)
    return setprecision(BigFloat, max(320, prec + 64)) do
        xs, ws = _legendre_nodes(n, lobatto)
        xs = (xs .+ 1) ./ 2; ws = ws ./ 2
        m = length(xs)
        cb = BigFloat[xs[1]; [xs[j] - xs[j - 1] for j in 2:m]; 1 - xs[m]]      # drifts before each kick and the last one: m + 1 values
        conv(x) = T === BigFloat ? BigFloat(x; precision = prec) : T(x)
        k = m + 1
        c = [conv(cb[i <= (k + 1) ÷ 2 ? i : k + 1 - i]) for i in 1:k]
        d = [i == k ? zero(T) : conv(ws[i <= (m + 1) ÷ 2 ? i : m + 1 - i]) for i in 1:k]
        return c, d
    end
end

"`SABA(T, n)`: the symmetric splitting A B A ... with kicks at the n Gauss-Legendre nodes (LaskarRobutel2001 Table I): order (2n, 2) in (tau, eps)."
function SABA(::Type{T}, n::Integer) where {T}
    c, d = _quadrature_table(T, n, false)
    return SplittingMethod(c, d; order = 2, palindromic = true, source = "LaskarRobutel2001", name = string("SABA", n), digits = 40)
end

"`SBAB(T, n)`: the symmetric splitting B A B ... with kicks at the n + 1 Gauss-Lobatto nodes (LaskarRobutel2001 Table I)."
function SBAB(::Type{T}, n::Integer) where {T}
    c, d = _quadrature_table(T, n, true)                         # kicks at 0 and 1: the first and the last drift are exactly 0
    return SplittingMethod(c, d; order = 2, palindromic = true, source = "LaskarRobutel2001", name = string("SBAB", n), digits = 40)
end

# Blanes2013 Table 3 (ABA864, 7 stages, generalized order (8, 6, 4)) and Table 4 (ABAH864, 8 stages): 40 digits; structure a1 b1 a2 b2 ... a_k b_k [a_(k+1)] b_k ... b1 a1 (Eq. (19)).
const ABA864_AB = (("0.0711334264982231177779387300061549964174", "0.241153427956640098736487795326289649618", "0.521411761772814789212136078067994229991", "-0.333698616227678005726562603400438876027"),
                   ("0.183083687472197221961703757166430291072", "0.310782859898574869507522291054262796375", "-0.0265646185119588006972121379164987592663", "0.0653961422823734184559721793911134363710"))
const ABAH864_AB = (("0.06810235651658372084723976682061164571212", "0.2511360387221033233072829580455350680082", "-0.07507264957216562516006821767601620052338",
                     "-0.009544719701745007811488218957217113269121", "0.5307579480704471776340674235341732001443"),
                    ("0.1684432593618954534310382697756917558148", "0.4243177173742677224300351657407231801453", "-0.5858109694681756812309015355404036521923", "0.4930499927320125053698281000239887162321"))

"`ABA864(T)`: Blanes2013 Table 3, 7 stages, generalized order (8, 6, 4)."
function ABA864(::Type{T}) where {T}
    a = [_P(T, s) for s in ABA864_AB[1]]; b = [_P(T, s) for s in ABA864_AB[2]]
    # The printed b_i of Table 3 sum to 1 only to 5.0e-31 (the a_i and every other table to 1e-39 or better); the consistency condition is restored through the central b4
    # (a change of 5.0e-31 in a coefficient of size 0.065, below the accuracy of the other conditions at double-double); recorded in the stage 03 notes.
    b[4] = 1 - 2 * (b[1] + b[2] + b[3])
    # a1 b1 a2 b2 a3 b3 a4 b4 a4 b3 a3 b2 a2 b1 a1: stages (a1,b1),(a2,b2),(a3,b3),(a4,b4),(a4,b3),(a3,b2),(a2,b1),(a1,0)
    c = [a[1], a[2], a[3], a[4], a[4], a[3], a[2], a[1]]; d = [b[1], b[2], b[3], b[4], b[3], b[2], b[1], zero(T)]
    return SplittingMethod(c, d; order = 4, palindromic = true, source = "Blanes2013", name = "ABA864", digits = 40)
end

"`ABAH864(T)`: Blanes2013 Table 4, 8 stages, generalized order (8, 6, 4), for use with an exact B flow here."
function ABAH864(::Type{T}) where {T}
    a = [_P(T, s) for s in ABAH864_AB[1]]; b = [_P(T, s) for s in ABAH864_AB[2]]
    # a1 b1 a2 b2 a3 b3 a4 b4 a5 b4 a4 b3 a3 b2 a2 b1 a1
    c = [a[1], a[2], a[3], a[4], a[5], a[4], a[3], a[2], a[1]]; d = [b[1], b[2], b[3], b[4], b[4], b[3], b[2], b[1], zero(T)]
    return SplittingMethod(c, d; order = 4, palindromic = true, source = "Blanes2013", name = "ABAH864", digits = 40)
end

# Corrector coefficients c of LaskarRobutel2001 Table II (the coefficient of {{A,B},B} in the modified Hamiltonian), 40 digits.
const CORR_SABA = Dict(4 => "0.003396775048208601331532157783492144", 5 => "0.002270543121419264819434955050039130", 6 => "0.001624459841624282521452258512463608",
                       7 => "0.001219643912760418472579211822331645", 8 => "0.000949308177745602234792177503535054", 9 => "0.000759846022860436646358196674176815",
                       10 => "0.000621934331486166426497049845358646")
const CORR_SBAB = Dict(5 => "0.002381486672953634187470386232181453", 6 => "0.001681346512091906326563693215296434", 7 => "0.001251765616039400003072516100251191",
                       8 => "0.000968797968073688571654684208462982", 9 => "0.000772349023999952078227686810260323", 10 => "0.000630320044163167840798638762665112")
# n = 1 to 3 (both) and SBAB n = 4 are the closed forms of Table II: SABA 1/12, (2 - sqrt 3)/24, (54 - 13 sqrt 15)/648; SBAB -1/24, 1/72, (13 - 5 sqrt 5)/288, (3861 - 791 sqrt 21)/64800.
function _corrector_c(::Type{T}, n::Integer, sbab::Bool) where {T}
    1 <= n <= 10 || throw(ArgumentError("the corrector coefficient is tabulated for n = 1 to 10 (LaskarRobutel2001 Table II), got $n"))
    sbab ? (n == 1 ? -T(1) / 24 : n == 2 ? T(1) / 72 : n == 3 ? (13 - 5 * sqrt(T(5))) / 288 : n == 4 ? (3861 - 791 * sqrt(T(21))) / 64800 : _P(T, CORR_SBAB[n])) :
           (n == 1 ? T(1) / 12 : n == 2 ? (2 - sqrt(T(3))) / 24 : n == 3 ? (54 - 13 * sqrt(T(15))) / 648 : _P(T, CORR_SABA[n]))
end

"""
    CorrectedMethod(base, c)

The corrected formula S_C = exp(-kappa L_F) S exp(-kappa L_F), F = {{A,B},B} = |grad V|^2, kappa = tau^3 eps^2 c/2 (LaskarRobutel2001 Eq. (62)); `base` a symmetric
table, `c` the coefficient of Table II. Built by [`SABAC`](@ref) and [`SBABC`](@ref).
"""
struct CorrectedMethod{T, M <: SplittingMethod{T}} <: AbstractIntegrator
    base::M
    c::T
    order::Int
end

"`SABAC(T, n)`: SABA_n with the corrector kick before and after; order 4 in eps for n >= 2 (LaskarRobutel2001 Sect. 10)."
SABAC(::Type{T}, n::Integer) where {T} = CorrectedMethod(SABA(T, n), _corrector_c(T, n, false), n >= 2 ? 4 : 2)
"`SBABC(T, n)`: SBAB_n with the corrector."
SBABC(::Type{T}, n::Integer) where {T} = CorrectedMethod(SBAB(T, n), _corrector_c(T, n, true), n >= 2 ? 4 : 2)

for f in (:SABA, :SBAB, :SABAC, :SBABC)
    @eval $f(n::Integer) = $f(Float64, n)
end
for f in (:ABA864, :ABAH864)
    @eval $f() = $f(Float64)
end
