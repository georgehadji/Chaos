using BenchmarkTools, Chaoslab, StaticArrays, StableRNGs

const SUITE = BenchmarkGroup()

# Henon-Heiles (F3): one order-8 step and tangent step, the crossing and escape solvers, a section map of 10^3 time units (plain and with the Jacobian), an escape orbit,
# the energy-surface starts and a portrait grid. Binary64 unless named.
const HH = HenonHeilesFlow()
const m8 = KahanLi(Float64, 8)
const tau = 1 / 20
const rng = StableRNG(20261012)
const E = 1 / 8
const (a, b) = section_interval(HH, E)
const zC = first(energy_start(HH, E, -0.015, 0.0))              # orbit C of Antonopoulos2010 Sect. 3.1 (p1 rebuilt)
const I4 = one(SMatrix{4, 4, Float64})
SUITE["KahanLi8 step"] = @benchmarkable advance($m8, $HH, $(Ref(zC))[], $tau)
SUITE["KahanLi8 tangent step"] = @benchmarkable tangent_advance($m8, $HH, $(Ref(zC))[], $I4, $tau)
SUITE["KahanLi8 step Double64"] = @benchmarkable advance($(KahanLi(Double64, 8)), $HH, $(Ref(SVector{4, Double64}(zC)))[], $(Double64(1) / 20))

# the bracket of a crossing and of an escape, found once at load time
function bracket(z, test)
    for _ in 1:10^5
        zn = advance(m8, HH, z, tau)
        test(z, zn) && return z
        z = zn
    end
    error("no bracket")
end
const zb = bracket(zC, (z, zn) -> z[1] < 0 <= zn[1] && zn[3] > 0)
const ze = bracket(first(energy_start(HH, 0.25, 0.9, 0.0)), (z, zn) -> hypot(z[1], z[2]) <= 10 < hypot(zn[1], zn[2]))
SUITE["locate_crossing"] = @benchmarkable locate_crossing($m8, $HH, $(Ref(zb))[], $tau)
SUITE["locate_escape"] = @benchmarkable locate_escape($m8, $HH, $(Ref(ze))[], $tau, 10.0)
SUITE["escape_sector"] = @benchmarkable escape_sector($HH, $(Ref(SVector(ze[1], ze[2])))[])

# one orbit: 1000 time units (2 10^4 steps) of the section map, and 200 time units of the escape tracker above the escape energy
SUITE["section_map T=1000"] = @benchmarkable section_map($m8, $HH, $zC, $tau, 20000)
SUITE["section_map T=1000 jacobian"] = @benchmarkable section_map($m8, $HH, $zC, $tau, 20000; jacobian = true)
SUITE["escape_orbit T=200 E=1/4"] = @benchmarkable escape_orbit($m8, $HH, $(first(energy_start(HH, 0.25, 0.9, 0.0))), $tau, 4000)

# starts: 10^4 rejection draws, the 100 x 100 portrait grid, the section interval
SUITE["random_energy_start 10^4 draws"] = @benchmarkable (for _ in 1:10^4; random_energy_start($HH, $E, ($a, $b), (-0.5, 0.5), $rng); end)
SUITE["energy_grid 100 x 100"] = @benchmarkable energy_grid($HH, $E, 100)
SUITE["section_interval"] = @benchmarkable section_interval($HH, $E)
