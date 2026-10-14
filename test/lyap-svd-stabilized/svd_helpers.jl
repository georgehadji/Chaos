# Helpers of the lyap-svd-stabilized suite (prefix sv_); the shared ones (gs_ prefix) are in lyap-gram-schmidt/gs_helpers.jl.
# Lyapunov names appear only inside function bodies and testsets.

const SV_GUARD_MESSAGE = "singular-value exponents are specified for d = 2 and k = d only; the general-d algorithm is held"

sv_run(sys, x0s, N; k = 2, tau, seed, tol = 8.0, W0 = nothing) =
    gs_serial(sys, x0s, N; k, tau, seed, reducers = (SVDSpectrum(; orthogonality_tol = tol),), W0)
sv_fin(res, i) = res.finals[i][1]

# 256-bit reference exponents (ln sigma_1, ln sigma_2) of a 2 x 2 product Y with det Y = 1 used analytically (spec Oracle):
# s = sqrt(||Y||_F^2 + 2) = sigma_1 + sigma_2, sigma_1 = (s + sqrt(s^2 - 4)) / 2, sigma_2 = 1 / sigma_1.
function sv_ref_logsigma(Y::Matrix{BigFloat})
    return setprecision(BigFloat, 256) do
        s = sqrt(sum(abs2, Y) + 2)
        l1 = log((s + sqrt(s^2 - 4)) / 2)
        (l1, -l1)
    end
end

# V9b: ln sigma_1 of [[1, n], [0, 1]], sigma_1 = (n + sqrt(n^2 + 4)) / 2
sv_V9b_lnsigma1(n) = setprecision(BigFloat, 256) do
    log((n + sqrt(BigFloat(n)^2 + 4)) / 2)
end

# Plain Float64 product M_{n-1} ... M_0 of the tangents at the stored points (the negative control of AC5).
function sv_naive_product(sys, pts, n)
    Y = SMatrix{2,2,Float64}(1, 0, 0, 1)
    for i in 1:n
        Y = tangent(sys, pts[i]) * Y
    end
    return Y
end

# Test map of dimension 4 whose forward and tangent throw: proves the scope guard acts before any work (AC12).
struct sv_Boom4 <: AbstractMap end
Chaoslab.dimension(::sv_Boom4) = 4
Chaoslab.forward(::sv_Boom4, z) = error("sv_Boom4: forward was called")
Chaoslab.tangent(::sv_Boom4, z) = error("sv_Boom4: tangent was called")
