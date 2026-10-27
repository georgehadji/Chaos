# Radix-2 FFT for the coarse search of the windowed Fourier frequency (spec naff-fma, Oracle: the transform only locates the starting bins;
# every value that decides a tie or a frequency is a direct sum). Single-threaded, so it adds no thread-count obligation (spec Compute).

# In-place forward transform a_k = sum_j a_j exp(-2 pi i j k / n) of a power-of-two length vector.
function _fft!(a::Vector{ComplexF64})
    n = length(a)
    ispow2(n) || throw(ArgumentError("FFT length must be a power of two, got $n"))
    j = 0
    for i in 1:(n - 1)
        bit = n >> 1
        while j & bit != 0
            j ⊻= bit
            bit >>= 1
        end
        j ⊻= bit
        i < j && ((a[i + 1], a[j + 1]) = (a[j + 1], a[i + 1]))
    end
    tw = [cispi(-2 * k / n) for k in 0:max(n ÷ 2 - 1, 0)]
    len = 2
    while len <= n
        half, stride = len >> 1, n ÷ len
        for s in 0:len:(n - 1), k in 0:(half - 1)
            u = a[s + k + 1]
            v = a[s + k + half + 1] * tw[k * stride + 1]
            a[s + k + 1] = u + v
            a[s + k + half + 1] = u - v
        end
        len <<= 1
    end
    return a
end
