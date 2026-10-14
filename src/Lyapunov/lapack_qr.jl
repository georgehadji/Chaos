# Library path (spec lyap-householder-qr Compute): LAPACK Householder QR for Float64, CPU only. The state is not isbits
# (a Matrix workspace allocated in reducer_init and reused by every handoff), so the stream resolves the backend to
# the CPU with the reason "reducer is CPU-only" without any code here.

"LAPACK Householder QR through LinearAlgebra (library path, Float64 only, CPU only), provenance \"V-H\" with implementation \"library\"."
struct LapackQR <: Orthonormalizer end

orthonormalizer_name(::LapackQR) = "V-H"

# Buffers of one orbit: the d x k matrix (W, then R and Q), the reflector scalars and the LAPACK work array.
struct LapackWork
    A::Matrix{Float64}
    tau::Vector{Float64}
    work::Vector{Float64}
end

LapackWork(d::Int, k::Int) =
    LapackWork(Matrix{Float64}(undef, d, k), Vector{Float64}(undef, k), Vector{Float64}(undef, max(1, 64 * max(d, k))))

const _BI = LinearAlgebra.BlasInt

# geqrf and orgqr straight through libblastrampoline with the preallocated work array (LinearAlgebra.LAPACK allocates one per call).
function _geqrf!(ws::LapackWork)
    m, n = size(ws.A)
    info = Ref{_BI}()
    ccall((LinearAlgebra.BLAS.@blasfunc(dgeqrf_), LinearAlgebra.BLAS.libblastrampoline), Cvoid,
          (Ref{_BI}, Ref{_BI}, Ptr{Float64}, Ref{_BI}, Ptr{Float64}, Ptr{Float64}, Ref{_BI}, Ptr{_BI}),
          m, n, ws.A, m, ws.tau, ws.work, length(ws.work), info)
    return info[]
end

function _orgqr!(ws::LapackWork)
    m, n = size(ws.A)
    info = Ref{_BI}()
    ccall((LinearAlgebra.BLAS.@blasfunc(dorgqr_), LinearAlgebra.BLAS.libblastrampoline), Cvoid,
          (Ref{_BI}, Ref{_BI}, Ref{_BI}, Ptr{Float64}, Ref{_BI}, Ptr{Float64}, Ptr{Float64}, Ref{_BI}, Ptr{_BI}),
          m, n, n, ws.A, m, ws.tau, ws.work, length(ws.work), info)
    return info[]
end

function _lapack_factor!(ws::LapackWork, W::SMatrix{d, k, Float64}) where {d, k}
    A = ws.A
    for j in 1:k, i in 1:d
        @inbounds A[i, j] = W[i, j]
    end
    _geqrf!(ws)
    R0 = SMatrix{k, k, Float64}(ntuple(l -> (i = mod(l - 1, k) + 1; j = div(l - 1, k) + 1; i > j ? 0.0 : @inbounds(A[i, j])), Val(k * k)))
    _orgqr!(ws)
    Q0 = SMatrix{d, k, Float64}(ntuple(l -> @inbounds(A[l]), Val(d * k)))
    return _finish(Q0, R0, _nref(Val(d), Val(k)))
end

orthonormalize(::LapackQR, W::SMatrix{d, k, Float64}) where {d, k} = _lapack_factor!(LapackWork(d, k), W)

# Non-isbits state: the shared QR accumulators plus the workspace.
struct LapackState{K}
    core::QRState{Float64, K}
    ws::LapackWork
end

function reducer_init(::QRSpectrum{LapackQR}, ::Type{T}, ::Val{d}, ::Val{k}) where {T, d, k}
    T === Float64 || throw(ArgumentError("LapackQR is the Float64 library path, got element type $T"))
    Threads.nthreads() > 1 && LinearAlgebra.BLAS.get_num_threads() != 1 && LinearAlgebra.BLAS.set_num_threads(1)
    sums = SVector{k, CompensatedSum{Float64}}(ntuple(_ -> CompensatedSum{Float64}(), Val(k)))
    return LapackState{k}(QRState{Float64, k}(sums, 0, 0.0, 0.0, 0, Int8(0), 0, 0.0, 0.0, 0), LapackWork(d, k))
end

function reducer_step(r::QRSpectrum{LapackQR}, st::LapackState, n::Int, x::SVector, W::SMatrix)
    st.core.code == 0 || return st, W
    core, Q = _qr_step(r, st.core, n, W, _lapack_factor!(st.ws, W))
    return LapackState(core, st.ws), Q
end

function reducer_update(r::QRSpectrum{LapackQR}, st::LapackState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::QRSpectrum{LapackQR}, st::LapackState) = st.core.code != 0
reducer_error(r::QRSpectrum{LapackQR}, st::LapackState, orbit, backend) = reducer_error(r, st.core, orbit, backend)
reducer_final(r::QRSpectrum{LapackQR}, st::LapackState) = reducer_final(r, st.core)
