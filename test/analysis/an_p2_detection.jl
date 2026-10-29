# D89h stable verdict, detection time, outcome, censored time, M_X (E28, M4). Spec: analysis-spec.md, Amendment D89, E28, criterion D89h.
# The oracle of the enumeration is the closed characterisation by the maximal constant decided suffix (a different algorithm from p2ref_detection).

# t_det by the suffix characterisation: the last verdict c9 must be decided; the suffix of verdicts equal to c9 starts at k0; stable from max(k0, first k with b_k <= N_c / divisor)
function p2_suffix_detection(V, rho, p)
    nb = length(V)
    last = V[nb]
    t = nothing
    if last !== :undecided
        k0 = nb
        while k0 > 1 && V[k0 - 1] === last
            k0 -= 1
        end
        p.budgets[k0] * p.persist_divisor <= p.N_c && (t = p.budgets[k0])
    end
    out = t === nothing ? :undetected : ((last === :chaotic && rho === :C) || (last === :regular && rho === :R)) ? :correct : :wrong
    return (t_det = t, outcome = out, censored = out === :correct ? t : p.censor_factor * p.N_c)
end

@testset "D89h detection: all 3^10 verdict sequences against the independent characterisation" begin
    p = p2_params()
    syms = (:chaotic, :regular, :undecided)
    n = 0; bad = 0; badref = 0
    counts = Dict(:correct => 0, :wrong => 0, :undetected => 0)
    for idx in Iterators.product(ntuple(_ -> 1:3, 10)...)
        V = Symbol[syms[i] for i in idx]
        for rho in (:C, :R)
            want = p2_suffix_detection(V, rho, p)
            got = AN.p2_detection(V, rho, p)
            n += 1
            bad += (got.t_det !== want.t_det || got.outcome !== want.outcome || got.censored != want.censored)
            if rho === :C
                counts[want.outcome] += 1
                badref += p2ref_detection(V, rho, p) != want              # the reference of the other tests agrees with the characterisation
            end
        end
    end
    @test n == 118_098 && bad == 0 && badref == 0
    # the number of sequences per outcome follows from the characterisation (rho = C): correct needs a chaotic last verdict (3^9 sequences) with a stable run from
    # k <= 7: the suffix of chaotic verdicts must reach index 7 (verdicts 7, 8, 9 chaotic and the sequence before it arbitrary): 3^7 = 2187 sequences
    @test counts[:correct] == 3^7 && counts[:wrong] == 3^7 && counts[:undetected] == 3^10 - 2 * 3^7
end

@testset "D89h detection: hand cases at the persistence limit and the median" begin
    p = p2_params()
    sq(s) = p2_sym(s)
    # b_7 = 5120 = N_c / 4: stability from k = 7 is allowed, from k = 8 (b_8 = 10240) is not
    @test AN.p2_detection(sq("rrrrrrrccc"), :C, p) == (t_det = 5120, outcome = :correct, censored = 5120)
    @test AN.p2_detection(sq("rrrrrrrrcc"), :C, p) == (t_det = nothing, outcome = :undetected, censored = 40960)
    # a flip at k = 8 gives no earlier stability
    @test AN.p2_detection(sq("cccccccrcc"), :C, p) == (t_det = nothing, outcome = :undetected, censored = 40960)
    @test AN.p2_detection(sq("cccccccccc"), :C, p) == (t_det = 40, outcome = :correct, censored = 40)
    # undecided at k = 9 is undetected, whatever came before; undecided inside the run moves the start after it
    @test AN.p2_detection(sq("cccccccccu"), :C, p).outcome === :undetected && AN.p2_detection(sq("cccccccccu"), :C, p).t_det === nothing
    @test AN.p2_detection(sq("uucccccccc"), :C, p) == (t_det = 160, outcome = :correct, censored = 160) && AN.p2_detection(sq("cucccccccc"), :C, p).t_det == 160
    # a stable verdict of the other class is wrong: the censored time is 2 N_c = 40960
    @test AN.p2_detection(sq("rrrrrrrrrr"), :C, p) == (t_det = 40, outcome = :wrong, censored = 40960) && AN.p2_detection(sq("cccrrrrrrr"), :R, p) == (t_det = 320, outcome = :correct, censored = 320)
    @test AN.p2_detection(sq("cccccccccc"), :R, p) == (t_det = 40, outcome = :wrong, censored = 40960)
    # the median of m values: the middle value for odd m, the mean of the two middle values for even m
    @test AN.p2_median([40, 40960, 80]) == 80 && AN.p2_median([40, 80, 160, 40960]) == 120 && AN.p2_median([5.0]) == 5.0 && AN.p2_median([2, 1]) == 1.5
    # M_X over fixture sets, by hand from the censored times of S2 (8 orbits): five with onset 4 (640), three with onset 6 (2560): the 4th and 5th are 640
    @test AN.p2_median(vcat(fill(640, 5), fill(2560, 3))) == 640 && AN.p2_median(vcat(fill(640, 4), fill(2560, 4))) == 1600
end
