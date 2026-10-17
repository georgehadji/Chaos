# Nested subsets S_D and S_L with deterministic caps (E17).

# Remove the members above `room` in the order largest rank first, ties to the last cell in canonical order.
function _trim(ms, room::Integer)
    drop = max(0, length(ms) - max(room, 0))
    order = sortperm(ms; by = m -> (-m.rank, -m.cell))
    removed = Set(order[1:drop])
    return [m for (i, m) in enumerate(ms) if !(i in removed)], [ms[i] for i in sort!(collect(removed))]
end

"""
    select_subsets(members; s_D, cap_D, s_L, cap_L) -> (SD, SL, removed_D, removed_L)

Subsets of E17 on members `(id, cell, rank, control)`: `S_D` holds the rank `1..s_D` members of every cell plus all controls, `S_L`
the rank `1..s_L` members that survived the `S_D` cap. A cap removes the largest rank first, ties to the last cell; controls are
never removed. Returns id lists.
"""
function select_subsets(members; s_D::Integer, cap_D::Integer, s_L::Integer, cap_L::Integer)
    ctrl = [m.id for m in members if m.control]
    keepD, remD = _trim([m for m in members if !m.control && m.rank <= s_D], cap_D - length(ctrl))
    keepL, remL = _trim([m for m in keepD if m.rank <= s_L], cap_L)
    return (SD = sort!(vcat([m.id for m in keepD], ctrl)), SL = sort!([m.id for m in keepL]),
            removed_D = sort!([m.id for m in remD]), removed_L = sort!([m.id for m in remL]))
end
