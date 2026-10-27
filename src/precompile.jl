# Precompile workload (decision D39): every test suite is compile-bound (57 to 88 percent of a fresh run is compilation, the cost
# growing with distinct element type x (d, k) x variant specializations), so the hot paths are compiled once at package precompile.
# Public API only, CPU backend only (never CUDA); BigFloat is left out (too many precision variants to cache usefully).
using PrecompileTools
using StaticArrays: SVector, SMatrix
using DoubleFloats: Double64

@setup_workload begin
    @compile_workload begin
        for T in (Float64, Double64)
            sm = StandardMap(exact(T, 2 * pi))
            x0 = [SVector{2, T}(exact(T, 1.0), exact(T, 2.0)), SVector{2, T}(exact(T, 0.5), exact(T, 3.0))]
            qrs = T === Float64 ? (CGS(), MGS(), RGS(), HouseholderQR(), AugmentedQR(), GivensQR(), LapackQR()) :
                                  (CGS(), MGS(), RGS(), HouseholderQR(), AugmentedQR(), GivensQR())
            initial_deviations(1, 1, Val(2), Val(2), T)
            propagate(sm, x0, 8; k = 1, tau = 2, seed = 1, reducers = (Benettin(),), backend = :cpu)
            propagate(sm, x0, 8; k = 2, tau = 2, seed = 1, reducers = (SVDSpectrum(),), backend = :cpu)
            propagate(sm, x0, 8; k = 2, tau = 2, seed = 1, reducers = (SALI(),), backend = :cpu)
            propagate(sm, x0, 8; k = 2, tau = 2, seed = 1, reducers = (FLI(replace = false), SALI()), backend = :cpu)
            propagate(sm, x0, 8; k = 2, tau = 2, seed = 1, reducers = (FLI(),), backend = :cpu)
            propagate(sm, x0, 8; k = 2, tau = 2, seed = 1, reducers = (MEGNO(replace = false), SALI()), backend = :cpu)
            propagate(sm, x0, 8; k = 2, tau = 2, seed = 1, reducers = (MEGNO(),), backend = :cpu)
            # rli: the orbit pair as one orbit of PairMap, block frame (xi_0, 0), (0, xi_0) (D75)
            z0 = [SVector{4, T}(exact(T, 1.0), exact(T, 2.0), exact(T, 1.0 + 2.0^-30), exact(T, 2.0))]
            o, l = one(T), zero(T)
            Wp = [SMatrix{4, 2, T}(o, l, l, l, l, l, o, l)]
            propagate(PairMap(sm), z0, 8; k = 2, tau = 2, seed = 1, reducers = (RLI(),), W0 = Wp, backend = :cpu)
            propagate(PairMap(sm), z0, 8; k = 2, tau = 2, seed = 1, reducers = (RLI(replace = false), MEGNO()), W0 = Wp, backend = :cpu)
            for alg in qrs, k in (1, 2)
                propagate(sm, x0, 8; k = k, tau = 2, seed = 1, reducers = (QRSpectrum(alg),), backend = :cpu)
            end
            # naff-fma: the weighted phase advance reader (2 window + 1 handoffs) and the series path on one window pair
            propagate(sm, x0, 33; k = 1, tau = 1, seed = 1, reducers = (WeightedPhaseAdvance(window = 16, observable = PlanarObservable(1, 2)),), backend = :cpu)
            # weighted-birkhoff: the reader (K = 2 observables, H = 2 horizons, 2 max(window) handoffs) and the pure function with its index
            wbr = WeightedBirkhoff(window = (16, 32), observables = (CosObservable(1), SinObservable(2)))
            propagate(sm, x0, 64; k = 1, tau = 1, seed = 1, reducers = (wbr,), backend = :cpu)
            wbf = weighted_birkhoff(T[sin(exact(T, 0.7) * j + k) for j in 0:63, k in 1:2], (16, 32))
            birkhoff_index.(wbf.avg1, wbf.avg2, wbf.S)
            frequency_map(Complex{T}[Complex(cos(exact(T, 0.7) * j), sin(exact(T, 0.7) * j)) for j in 0:63], 32)
            propagate_shadow(sm, x0, 8; d0 = exact(T, 2.0^-20), tau = 2, seed = 1, backend = :cpu)
            # d = 4 linear map, k = 4 (the (4, 4) panel of the QR tests)
            lm = LinearSymplecticMap(block_diagonal(rotation_block(T, exact(T, 0.7)), hyperbolic_block(T, exact(T, 0.3))))
            y0 = [SVector{4, T}(exact(T, 0.1), exact(T, 0.2), exact(T, 0.3), exact(T, 0.4))]
            for alg in qrs, k in (2, 4)
                propagate(lm, y0, 8; k = k, tau = 2, seed = 1, reducers = (QRSpectrum(alg),), backend = :cpu)
            end
        end
    end
end
