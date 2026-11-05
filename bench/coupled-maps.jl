using BenchmarkTools, Chaoslab, StaticArrays, LinearAlgebra

const SUITE = BenchmarkGroup()

const FM = FroeschleMap(0.5, 0.1, 1.0e-3)
const WFM = WrappedFroeschleMap(0.5, 0.1, 1.0e-3)
const ZF = SVector(0.5, 0.5, 0.0, 0.0)
const SF = SVector(0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

SUITE["y3 forward"] = @benchmarkable forward($FM, $(Ref(ZF))[])
SUITE["y3 tangent"] = @benchmarkable tangent($FM, $(Ref(ZF))[])
SUITE["y3 inverse"] = @benchmarkable inverse($FM, $(Ref(ZF))[])
SUITE["y3 wrapped forward"] = @benchmarkable forward($WFM, $(Ref(SF))[])
SUITE["y3 forward 1000 steps"] = @benchmarkable begin
    z = $ZF
    for _ in 1:1000
        z = forward($FM, z)
    end
    z
end
SUITE["y3 tangent frame (4 x 4)"] = @benchmarkable tangent($FM, $(Ref(ZF))[]) * $(Ref(SMatrix{4, 4, Float64}(1.0I)))[]

for N in (2, 10, 100)
    sys = StandardMapChain(N, 3.0, 0.6102; boundary = N >= 3 ? :ring : :open)
    wsys = WrappedStandardMapChain(N, 3.0, 0.6102; boundary = N >= 3 ? :ring : :open)
    z = [0.1 * sin(j) + 1.0 for j in 1:(2N)]
    s = vcat(mod.(z, 2pi), zeros(2N))
    out, sout, M = similar(z), similar(s), zeros(2N, 2N)
    W, WO = [cos(i + 3j) for i in 1:(2N), j in 1:4], zeros(2N, 4)
    SUITE["y4 N=$N forward!"] = @benchmarkable forward!($out, $sys, $z)
    SUITE["y4 N=$N wrapped forward!"] = @benchmarkable forward!($sout, $wsys, $s)
    SUITE["y4 N=$N tangent! (dense)"] = @benchmarkable tangent!($M, $sys, $z)
    SUITE["y4 N=$N tangent_action! (4 vectors)"] = @benchmarkable tangent_action!($WO, $sys, $z, $W)
    SUITE["y4 N=$N inverse_tangent_action! (4 vectors)"] = @benchmarkable inverse_tangent_action!($WO, $sys, $z, $W)
    SUITE["y4 N=$N dense tangent times frame"] = @benchmarkable mul!($WO, $M, $W) setup = (tangent!($M, $sys, $z))
end
