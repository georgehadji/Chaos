# Amplification sum Psi_n of lyap-benettin S3 (Amendment D35), shared by lyap-benettin, lyap-householder-qr and lyap-two-trajectory.
# Psi_n = sum_{j=1..n} norm(M_n ... M_j) norm(w_{j-1}) / norm(w_n), M_i the one-step tangent at realized point i-1, w_i = M_i w_{i-1}.
# Float64 throughout (a scale for a tolerance, no precision needed); scaled recursions, so no overflow for chaotic orbits.
using StaticArrays, LinearAlgebra

const LB_PSI_LOADED = true                     # include guard (other slugs include this file once)

lb_norm2(P::SMatrix{2, 2, Float64}) = (F = sum(abs2, P); d = det(P); sqrt((F + sqrt(max(F^2 - 4 * d^2, 0.0))) / 2))   # matrix 2-norm

# Tangents M_1 .. M_n of the orbit from x0 in the element type of x0 (iterating `forward` is bitwise the stream's orbit), widened to Float64.
function lb_tangents(sys, x0, n::Integer)
    Ms = Vector{SMatrix{2, 2, Float64, 4}}(undef, n)
    z = x0
    for i in 1:n
        Ms[i] = SMatrix{2, 2, Float64, 4}(map(Float64, tangent(sys, z)))
        z = forward(sys, z)
    end
    return Ms
end

# Psi_n for the tangents Ms = [M_1 .. M_n] and the start vector w0.
function lb_Psi(Ms::AbstractVector{<:SMatrix{2, 2}}, w0)
    n = length(Ms)
    lw = zeros(n + 1)                          # lw[k + 1] = ln norm(w_k)
    w = SVector{2, Float64}(Float64(w0[1]), Float64(w0[2]))
    a = norm(w)
    w = w / a
    acc = log(a)
    lw[1] = acc
    for k in 1:n
        w = SMatrix{2, 2, Float64}(Ms[k]) * w
        a = norm(w)
        acc += log(a)
        w = w / a
        lw[k + 1] = acc
    end
    P = one(SMatrix{2, 2, Float64})            # P = exp(lP) * M_n ... M_j, kept at unit Frobenius norm
    lP = 0.0
    psi = 0.0
    for j in n:-1:1
        Q = P * SMatrix{2, 2, Float64}(Ms[j])
        f = norm(Q)
        P = Q / f
        lP += log(f)
        psi += exp(lP + log(lb_norm2(P)) + lw[j] - lw[n + 1])
    end
    return psi
end
lb_Psi(sys, x0, w0, n::Integer) = lb_Psi(lb_tangents(sys, x0, n), w0)
