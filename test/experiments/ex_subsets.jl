# AC21 subsets (E17), synthetic part. The build-level part (nesting, caps, manifest removed ids) is in ex_panel.jl.
# select_subsets takes members (id, cell, rank, control) with cell the position of the panel cell in canonical order (E20).

function ex_population(; ncells = 27, per_cell = 40, ncontrols = 13, short = Dict{Int, Int}())
    members = NamedTuple{(:id, :cell, :rank, :control), Tuple{Int, Int, Int, Bool}}[]
    id = 0
    for c in 1:ncells, r in 1:get(short, c, per_cell)
        push!(members, (id = (id += 1), cell = c, rank = r, control = false))
    end
    for _ in 1:ncontrols
        push!(members, (id = (id += 1), cell = 0, rank = 0, control = true))
    end
    return members
end
ex_id(ms, c, r) = only(m.id for m in ms if m.cell == c && m.rank == r)

@testset "AC21 subsets (synthetic)" begin
    ms = ex_population()
    ctrl = [m.id for m in ms if m.control]
    @assert length(ms) == 27 * 40 + 13
    # 27 cells x 7 + 13 controls = 202 members before the cap; cap 200 removes the two rank-7 members of the last two cells
    sub = EX.select_subsets(ms; s_D = 7, cap_D = 200, s_L = 4, cap_L = 110)
    @test sort(sub.removed_D) == sort([ex_id(ms, 27, 7), ex_id(ms, 26, 7)])
    @test length(sub.SD) == 200 && issubset(ctrl, sub.SD)               # controls are never removed
    @test isempty(sub.removed_L) && length(sub.SL) == 108               # 4 x 27 = 108 <= 110, all ranks 1..4, no controls
    @test issubset(sub.SL, sub.SD) && isdisjoint(sub.SL, ctrl)
    @test all(m -> m.id in sub.SL, (m for m in ms if !m.control && m.rank <= 4))
    # before the cap every populated cell contributes exactly s_D members (all of them when it has fewer)
    nocap = EX.select_subsets(ms; s_D = 7, cap_D = 10_000, s_L = 4, cap_L = 10_000)
    @test length(nocap.SD) == 202 && isempty(nocap.removed_D) && length(nocap.SL) == 108
    for c in 1:27
        @test count(m -> m.id in nocap.SD && m.cell == c, ms) == 7
        @test count(m -> m.id in nocap.SL && m.cell == c, ms) == 4
    end
    # the ladder cap: S_L is formed from S_D after its cap; 108 members, cap 100 removes the rank-4 members of the last 8 cells
    sub = EX.select_subsets(ms; s_D = 7, cap_D = 200, s_L = 4, cap_L = 100)
    @test length(sub.SL) == 100 && issubset(sub.SL, sub.SD)
    @test sort(sub.removed_L) == sort([ex_id(ms, c, 4) for c in 20:27])
    # the removed rank-7 members of the S_D cap are not in S_L anyway; nesting and caps hold
    @test length(sub.SD) <= 200 && length(sub.SL) <= 100
    # a cell with fewer members than s_D contributes all it has; ties go to the last cell in canonical order
    ms2 = ex_population(; ncells = 5, per_cell = 4, ncontrols = 2, short = Dict(2 => 1, 4 => 2))
    s2 = EX.select_subsets(ms2; s_D = 3, cap_D = 100, s_L = 2, cap_L = 100)
    @test count(m -> m.id in s2.SD && m.cell == 2, ms2) == 1 && count(m -> m.id in s2.SD && m.cell == 4, ms2) == 2
    s3 = EX.select_subsets(ms2; s_D = 3, cap_D = 12, s_L = 2, cap_L = 100)       # 3+1+3+2+3 = 12 members + 2 controls = 14 > 12
    @test length(s3.SD) == 12 && issubset([m.id for m in ms2 if m.control], s3.SD)
    @test sort(s3.removed_D) == sort([ex_id(ms2, 5, 3), ex_id(ms2, 3, 3)])      # largest rank 3 first, last cell first
    # a cap equal to the number of controls leaves only controls
    s4 = EX.select_subsets(ms2; s_D = 3, cap_D = 2, s_L = 2, cap_L = 2)
    @test length(s4.SD) == 2 && isempty(s4.SL)
end
