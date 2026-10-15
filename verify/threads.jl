# Thread-count invariance of the Batch B CPU estimators (lyap-benettin AC10b, lyap-gram-schmidt, lyap-svd-stabilized AC9,
# lyap-householder-qr AC9): every output of `propagate` hashed in child processes run with 1 and 4 threads, compared bitwise.
# lyap-two-trajectory AC11 does the same in its own script. Log: workspaces/module-build/stages/04-verify/output/threads-verify-log.txt
using Chaoslab, StaticArrays, StableRNGs, DoubleFloats

const SEED = 20261012
const CASES = [("lyap-benettin AC10b", "Benettin", () -> (Benettin(),), 1),
               ("lyap-gram-schmidt", "QRSpectrum(RGS)", () -> (QRSpectrum(RGS()),), 2),
               ("lyap-svd-stabilized AC9", "SVDSpectrum", () -> (SVDSpectrum(),), 2),
               ("lyap-householder-qr AC9", "QRSpectrum(HouseholderQR)", () -> (QRSpectrum(HouseholderQR()),), 2),
               ("lyap-householder-qr AC9", "QRSpectrum(GivensQR)", () -> (QRSpectrum(GivensQR()),), 2),
               ("lyap-householder-qr AC9", "QRSpectrum(AugmentedQR)", () -> (QRSpectrum(AugmentedQR()),), 2),
               ("lyap-householder-qr AC9", "QRSpectrum(LapackQR)", () -> (QRSpectrum(LapackQR()),), 2)]

function case_hash(reducers, k, ::Type{T}) where {T}
    rng = StableRNG(SEED + 60000)
    x0 = [SVector{2, T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:2048]
    res = propagate(StandardMap(exact(T, 2pi)), x0, 1000; k = k, tau = 2, seed = SEED, reducers = reducers, backend = :cpu)
    return hash((res.x, res.W, res.finals))
end

types(name) = name == "QRSpectrum(LapackQR)" ? (Float64,) : (Float64, Double64)

if get(ENV, "THREADS_CHILD", "") == "1"
    for (_, name, red, k) in CASES, T in types(name)
        println("HASH $name $T ", case_hash(red(), k, T))
    end
    exit(0)
end

hashes(nt) = Dict(m[1] => m[2] for m in eachmatch(r"HASH (.+) (\d+)\n",
    read(pipeline(addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "THREADS_CHILD" => "1");
                  stderr = devnull), String)))

out = joinpath(@__DIR__, "..", "workspaces", "module-build", "stages", "04-verify", "output", "threads-verify-log.txt")
h1, h4 = hashes(1), hashes(4)
open(out, "w") do io
    println(io, "# thread-count invariance, julia $(VERSION), child processes -t 1 and -t 4, 2048 orbits, K = 2 pi, n = 1e3, tau = 2, base seed $SEED")
    println(io, "# format: AC | quantity | value | bound | verdict")
    nfail = 0
    for (ac, name, _, k) in CASES, T in types(name)
        key = "$name $T"
        ok = haskey(h1, key) && get(h1, key, "") == get(h4, key, "-")
        nfail += !ok
        println(io, "$ac | $name, $T, k = $k: hash of x, W and finals, 1 thread / 4 threads | $(get(h1, key, "missing")) / $(get(h4, key, "missing")) | equal | $(ok ? "PASS" : "FAIL")")
    end
    println(io, "# SUMMARY FAIL=$nfail")
end
print(read(out, String))
