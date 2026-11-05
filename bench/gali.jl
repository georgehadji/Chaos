using BenchmarkTools, Chaoslab, StaticArrays, LinearAlgebra

const SUITE = BenchmarkGroup()

const r = GALI()
const x4 = SVector(0.5, 0.5, 0.0, 0.0)

# Reducer step (one handoff: K norms, K logs, nested QR) and the pure evaluation, Y3 (d = 4) and Y4 N = 100 (d = 200), K = 4.
for (name, d, T) in (("y3", 4, Float64), ("y3", 4, Double64), ("y4 N=100", 200, Float64), ("y4 N=100", 200, Double64))
    W = SMatrix{d, 4, T}([cos(i + 3j) + (i == j) for i in 1:d, j in 1:4])
    x = zeros(SVector{d, T})
    st = reducer_init(r, T, Val(d), Val(4))
    SUITE["$name $T step K=4"] = @benchmarkable reducer_step($r, $st, 1, $x, $(Ref(W))[])
    SUITE["$name $T gali_values K=4"] = @benchmarkable gali_values($(Ref(W))[])
end

W44 = SMatrix{4, 4, Float64}([cos(i + 3j) + (i == j) for i in 1:4, j in 1:4])
SUITE["y3 Float64 step K=d=4 (log form)"] = @benchmarkable reducer_step($r, $(reducer_init(r, Float64, Val(4), Val(4))), 1, $x4, $(Ref(W44))[])
SUITE["check_frame Float64 K=4"] = @benchmarkable check_frame($r, $(Ref(SMatrix{4, 4, Float64}(1.0 * I(4))))[])

const FM = FroeschleMap(0.5, 0.1, 1.0e-3)
SUITE["y3 orbit 1000 steps K=4"] = @benchmarkable propagate_orbit($FM, $x4, $(SMatrix{4, 4, Float64}(1.0 * I(4))), 1000, 1, ($r,))
