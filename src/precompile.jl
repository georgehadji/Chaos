# Precompile workload (decision D39): every test suite is compile-bound (57 to 88 percent of a fresh run is compilation, the cost
# growing with distinct element type x (d, k) x variant specializations), so the hot paths are compiled once at package precompile.
# Public API only, CPU backend only (never CUDA); BigFloat is left out (too many precision variants to cache usefully).
using PrecompileTools
using StaticArrays: SVector
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
            for alg in qrs, k in (1, 2)
                propagate(sm, x0, 8; k = k, tau = 2, seed = 1, reducers = (QRSpectrum(alg),), backend = :cpu)
            end
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
