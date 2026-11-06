using BenchmarkTools, Chaoslab, StaticArrays, LinearAlgebra

const SUITE = BenchmarkGroup()

# One step of the table methods on the pendulum (F4) and the planar Kepler problem (F2), binary64 and double-double: leapfrog, Forest-Ruth, PEFRL, the 8th order Kahan-Li table,
# the 10th order Sofroniou-Spaletta table, SABA3 and the corrected SABAC3.
const PEND = Pendulum(1.0)
const KEPL = Kepler()
const zp = SVector(1.0, 0.3)
const zk = SVector(1.0, 0.0, 0.0, 1.0)
for (name, m) in (("leapfrog", Leapfrog()), ("ForestRuth", ForestRuth()), ("PEFRL", PEFRL()), ("KahanLi8", KahanLi(8)), ("SofroniouSpaletta10", SofroniouSpaletta10()),
                  ("SABA3", SABA(3)), ("SABAC3", SABAC(3)))
    SUITE["$name step pendulum Float64"] = @benchmarkable advance($m, $PEND, $(Ref(zp))[], 0.05)
    SUITE["$name step Kepler Float64"] = @benchmarkable advance($m, $KEPL, $(Ref(zk))[], 0.01)
    SUITE["$name tangent step Kepler Float64"] = @benchmarkable tangent_advance($m, $KEPL, $(Ref(zk))[], $(Ref(SMatrix{4, 4, Float64}(1.0 * I(4))))[], 0.01)
end
const md = KahanLi(Double64, 8)
SUITE["KahanLi8 step Kepler Double64"] = @benchmarkable advance($md, $KEPL, $(Ref(SVector{4, Double64}(1, 0, 0, 1)))[], Double64(0.01))

# Long runs: 10^4 steps of the leapfrog (n + 1 force evaluations in the kick-drift-kick form) and of KahanLi8, plain and with compensated summation (Hairer2008 Sect. 2)
const lf = Leapfrog(Float64, :KDK)
SUITE["leapfrog KDK 10^4 steps Kepler"] = @benchmarkable advance_n($lf, $KEPL, $zk, 0.01, 10000)
SUITE["KahanLi8 10^4 steps Kepler"] = @benchmarkable advance_n($(KahanLi(8)), $KEPL, $zk, 0.01, 10000)
SUITE["KahanLi8 10^4 steps Kepler compensated"] = @benchmarkable advance_n($(KahanLi(8)), $KEPL, $zk, 0.01, 10000; compensated = true)

# Implicit: Gauss s = 2 (order 4, fixed-point iteration to round-off) and implicit midpoint, one step on the pendulum
for s in (1, 2, 4)
    g = Gauss(Float64, s)
    SUITE["Gauss s=$s step pendulum Float64"] = @benchmarkable advance($g, $PEND, $(Ref(zp))[], 0.1)
end

# Tao's extended phase-space method on T1, order 2 and 4
const TAO = TaoNonseparable()
for order in (2, 4)
    t = TaoExtended(Float64, 10.0; order = order)
    SUITE["Tao order $order step T1 Float64"] = @benchmarkable advance($t, $TAO, $(Ref(extend_state(SVector(1.0, 0.5))))[], 0.01)
end

# Ensembles of independent orbits on the CPU backend: 1024 pendulum orbits, 100 steps (leapfrog, KahanLi8, Gauss s = 2)
const x0 = [SVector(0.002 * i, 0.3) for i in 1:1024]
for (name, m) in (("leapfrog", Leapfrog()), ("KahanLi8", KahanLi(8)), ("Gauss2", Gauss(Float64, 2)))
    SUITE["advance_ensemble $name 1024 orbits 100 steps pendulum"] = @benchmarkable advance_ensemble($m, $PEND, $x0, 0.05, 100; backend = :cpu)
end
