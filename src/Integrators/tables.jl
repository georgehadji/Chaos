# Explicit tables of the catalog rows S1 to S9 (catalog integrators.md). Decimal constants are copied from the cited source as strings and parsed in the
# element type (`parse(T, s)`: no Float64 detour); algebraic constants are computed in the element type. Every constructor is also callable without
# the element type (Float64).

_P(::Type{T}, s::AbstractString) where {T} = parse(T, s)

# S1: first order, one stage: Yoshida1990 Eq. (2.2) (drift then kick, :DK); Ruth1983 Eqs. (10) to (14) (kick, then drift, :KD). The two orderings are adjoint.
function SymplecticEuler(::Type{T}, order::Symbol = :DK) where {T}
    order in (:DK, :KD) || throw(ArgumentError("SymplecticEuler order must be :DK or :KD, got :$order"))
    c, d = order === :DK ? ([one(T)], [one(T)]) : ([zero(T), one(T)], [one(T), zero(T)])
    return SplittingMethod(c, d; order = 1, palindromic = false, source = "Yoshida1990,Ruth1983", name = string("symplectic Euler ", order))
end

# S2: leapfrog: Yoshida1990 Eq. (2.3) (drift-kick-drift), Ruth1983 Eqs. (15) to (18) with a = b = 1/2, ForestRuth1990 Eq. (4.2) (kick-drift-kick).
function Leapfrog(::Type{T}, order::Symbol = :DKD) where {T}
    h = inv(T(2))
    c, d = order === :DKD ? ([h, h], [one(T), zero(T)]) : order === :KDK ? ([zero(T), one(T), zero(T)], [h, h, zero(T)]) :
           throw(ArgumentError("Leapfrog order must be :DKD or :KDK, got :$order"))
    return SplittingMethod(c, d; order = 2, palindromic = true, source = "Yoshida1990,Ruth1983,ForestRuth1990", name = string("leapfrog ", order))
end

# S3: McLachlan1995 Eq. (4.2): exp(z tau A) exp(tau B/2) exp((1 - 2z) tau A) exp(tau B/2) exp(z tau A), z = (y^2 + 6y - 2)/(12 y), y = (2 sqrt(326) - 36)^(1/3)
# (Sect. 4, Order 2; Table 2, order 2, S, m = 2).
function McLachlan2(::Type{T}) where {T}
    y = (2 * sqrt(T(326)) - 36)^inv(T(3)); z = (y^2 + 6y - 2) / (12y); h = inv(T(2))
    return SplittingMethod([z, 1 - 2z, z], [h, h, zero(T)]; order = 2, palindromic = true, source = "McLachlan1995", name = "McLachlan optimized second order")
end

# S4: Ruth1983 Eq. (29) (kick with weight c_i, then drift with weight d_i), solution Eq. (31): c = (7/24, 3/4, -1/24), d = (2/3, -2/3, 1); the same scheme
# and conditions are ForestRuth1990 Eqs. (4.3), (4.4). Stored kick first: stage i is (drift d_(i-1), kick c_i).
function Ruth3(::Type{T}) where {T}
    return SplittingMethod([zero(T), T(2) / 3, -T(2) / 3, one(T)], [T(7) / 24, T(3) / 4, -T(1) / 24, zero(T)]; order = 3, palindromic = false,
                           source = "Ruth1983,ForestRuth1990", name = "Ruth third order")
end

# S5: Forest-Ruth, the triple jump of the leapfrog with theta = 1/(2 - 2^(1/3)): Omelyan2002 Eq. (8); ForestRuth1990 Eq. (5.32) (the same map, Sect. 4 solution Eq. (4.8)
# to (4.12)); closed form Yoshida1990 Eq. (2.11).
function ForestRuth(::Type{T}) where {T}
    m = triple_jump(Leapfrog(T, :DKD), 4)
    return SplittingMethod(collect(m.c), collect(m.d); order = 4, palindromic = true, source = "ForestRuth1990,Omelyan2002,Yoshida1990", name = "Forest-Ruth")
end

# S8: PEFRL (position version), Omelyan2002 Eq. (22) with the parameters of Eq. (20): drifts xi, chi, 1 - 2(chi + xi), chi, xi; kicks (1 - 2 lambda)/2, lambda, lambda, (1 - 2 lambda)/2.
function PEFRL(::Type{T}) where {T}
    xi = _P(T, "0.1786178958448091"); lam = _P(T, "-0.2123418310626054"); chi = _P(T, "-0.06626458266981849")
    kd = (1 - 2lam) / 2
    return SplittingMethod([xi, chi, 1 - 2 * (chi + xi), chi, xi], [kd, lam, lam, kd, zero(T)]; order = 4, palindromic = true, source = "Omelyan2002",
                           name = "PEFRL", digits = 16)
end

# S9: CalvoSanzSerna1993a Eq. (2.8) (five-stage, fourth order, first-same-as-last RKN, 18 digits): gamma_i and b_i; as a splitting: kick b_i then drift gamma_(i+1) - gamma_i,
# stored with the drift in front: c_i = gamma_i - gamma_(i-1), c_1 = gamma_1 = 0 (Eqs. (2.2), (2.3), (2.4), (2.5): alpha_ij = b_j (gamma_i - gamma_j)).
function CalvoSanzSerna4(::Type{T}) where {T}
    g = [_P(T, s) for s in ("0", "0.205177661542286386", "0.608198943146500973", "0.487278066807586965", "1")]
    b = [_P(T, s) for s in ("0.061758858135626325", "0.338978026553643355", "0.614791307175577566", "-0.140548014659373380", "0.125019822794526133")]
    c = [g[1]; [g[i] - g[i - 1] for i in 2:5]]
    return SplittingMethod(c, b; order = 4, palindromic = false, source = "CalvoSanzSerna1993a", name = "Calvo-Sanz-Serna fourth order", digits = 18)
end

# S11: Yoshida1990 Tables 1 and 2: the weights w_1..w_m of the symmetric composition (5.1) of leapfrogs, w_0 = 1 - 2 sum w_i (Eq. (5.18)). 6th order solution A (Table 1);
# 8th order solutions A and D (Table 2); 15 digits as printed.
const YOSHIDA_W = Dict(
    (6, :A) => ("-0.117767998417887E1", "0.235573213359357E0", "0.784513610477560E0"),
    (8, :A) => ("-0.161582374150097E1", "-0.244699182370524E1", "-0.716989419708120E-2", "0.244002732616735E1", "0.157739928123617E0", "0.182020630970714E1", "0.104242620869991E1"),
    (8, :D) => ("0.102799849391985E0", "-0.196061023297549E1", "0.193813913762276E1", "-0.158240635368243E0", "-0.144485223686048E1", "0.253693336566229E0", "0.914844246229740E0"))

function Yoshida(::Type{T}, order::Integer, solution::Symbol = :A) where {T}
    haskey(YOSHIDA_W, (Int(order), solution)) || throw(ArgumentError("Yoshida: no solution $solution of order $order (6 :A, 8 :A, 8 :D)"))
    w = [_P(T, s) for s in YOSHIDA_W[(Int(order), solution)]]
    half = vcat(reverse(w), 1 - 2 * sum(w))          # w_m .. w_1, w_0 (the middle piece last)
    return _palindromic_leapfrog(T, half, Int(order), string("Yoshida ", order, solution), "Yoshida1990", 15)
end

# S12: KahanLi1997 Appendix A, the gamma_j of the palindromic compositions of leapfrogs (Eq. (4)): s7odr6 (equal to Yoshida1990 Table 1 solution A), s9odr6a, s17odr8a; 20 digits.
const KAHANLI_G = Dict(
    :s7odr6 => (6, ("0.78451361047755726382", "0.23557321335935813368", "-1.1776799841788710069", "1.3151863206839112189")),
    :s9odr6a => (6, ("0.39216144400731413928", "0.33259913678935943860", "-0.70624617255763935981", "0.082213596293550800230", "0.79854399093482996340")),
    :s17odr8a => (8, ("0.13020248308889008088", "0.56116298177510838456", "-0.38947496264484728641", "0.15884190655515560090", "-0.39590389413323757734",
                      "0.18453964097831570709", "0.25837438768632204729", "0.29501172360931029887", "-0.60550853383003451170")))

function KahanLi(::Type{T}, scheme::Union{Symbol, Integer} = 6) where {T}
    key = scheme isa Symbol ? scheme : (scheme == 6 ? :s9odr6a : scheme == 8 ? :s17odr8a : throw(ArgumentError("KahanLi: scheme 6, 8 or :s7odr6")))
    haskey(KAHANLI_G, key) || throw(ArgumentError("KahanLi: unknown scheme :$key"))
    order, g = KAHANLI_G[key]
    return _palindromic_leapfrog(T, [_P(T, s) for s in g], order, string("Kahan-Li ", key), "KahanLi1997", 20)
end

# S13: SofroniouSpaletta2005 Eq. (17), the best 35-stage tenth order method (p. 605): gamma_1..gamma_18 (gamma_i = gamma_(36-i)), 32 digits as printed.
const SOFSPA_G = ("0.078795722521686419263907679337684", "0.31309610341510852776481247192647", "0.027918383235078066109520273275299",
                  "-0.22959284159390709415121339679655", "0.13096206107716486317465685927961", "-0.26973340565451071434460973222411",
                  "0.074973343155891435666137105641410", "0.11199342399981020488957508073640", "0.36613344954622675119314812353150",
                  "-0.39910563013603589787862981058340", "0.10308739852747107731580277001372", "0.41143087395589023782070411897608",
                  "-0.0048663605831352617621956593099771", "-0.39203335370863990644808193642610", "0.051942502962449647037182904015976",
                  "0.050665090759924496335874344156866", "0.049674370639729879054568800279461", "0.049317735759594537917680008339338")
SofroniouSpaletta10(::Type{T}) where {T} = _palindromic_leapfrog(T, [_P(T, s) for s in SOFSPA_G], 10, "Sofroniou-Spaletta 10", "SofroniouSpaletta2005", 32)

SymplecticEuler(order::Symbol) = SymplecticEuler(Float64, order)      # callable without the element type (Float64)
Leapfrog(order::Symbol) = Leapfrog(Float64, order)
Yoshida(order::Integer, solution::Symbol = :A) = Yoshida(Float64, order, solution)
KahanLi(scheme::Union{Symbol, Integer} = 6) = KahanLi(Float64, scheme)
for f in (:SymplecticEuler, :Leapfrog, :McLachlan2, :Ruth3, :ForestRuth, :PEFRL, :CalvoSanzSerna4, :SofroniouSpaletta10)
    @eval $f() = $f(Float64)
end
