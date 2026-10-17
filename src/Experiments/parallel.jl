# Work distribution over independent items. Results are written by item index, so the outcome never depends on the
# number of tasks or on which task ran an item (spec I8).

# One worker: builds its state once, then takes item indices from the shared counter until none is left.
function _work(f, make, next, n::Integer)
    s = make()
    while true
        i = Threads.atomic_add!(next, 1)
        i > n && break
        f(s, i)
    end
    return nothing
end

"""
    _par(f, make, n, nt)

Call `f(state, i)` for `i = 1:n` on at most `nt` tasks; each task builds its own `state = make()` once.
"""
function _par(f, make, n::Integer, nt::Integer)
    nt = max(1, min(nt, n))
    next = Threads.Atomic{Int}(1)
    if nt == 1
        _work(f, make, next, n)
    else
        @sync for _ in 1:nt
            Threads.@spawn _work(f, make, next, n)
        end
    end
    return nothing
end

_ntasks(nt) = nt === nothing ? Threads.nthreads() : max(1, Int(nt))
