# Murray's law radius propagation

"""
    update_radii!(tree::VascularTree, gamma::Float64)

Propagate radii from terminals to root using Murray's law:
  r_parent^gamma = sum(r_child^gamma)

Terminal radii are kept fixed. Internal segments' radii are recomputed
bottom-up (from deepest generation to root).
"""
function update_radii!(tree::VascularTree, gamma::Float64)
    topo = tree.topology
    seg = tree.segments
    n = seg.n
    n == 0 && return nothing

    # BFS order (root to leaves), then reverse for bottom-up
    # sequential: tree traversal — O(n) instead of O(n * max_gen)
    root_id = Int(tree.root_segment_id)
    root_id < 1 && return nothing
    order = Vector{Int}()
    sizehint!(order, n)
    queue = Int[root_id]
    while !isempty(queue)
        id = popfirst!(queue)
        push!(order, id)
        c1 = Int(topo.child1_id[id])
        c2 = Int(topo.child2_id[id])
        c3 = Int(topo.child3_id[id])
        c1 > 0 && push!(queue, c1)
        c2 > 0 && push!(queue, c2)
        c3 > 0 && push!(queue, c3)
    end

    # Bottom-up: traverse in reverse BFS order (leaves first)
    for k in length(order):-1:1
        i = order[k]
        topo.is_terminal[i] && continue

        r_gamma_sum = 0.0
        c1 = topo.child1_id[i]
        c2 = topo.child2_id[i]
        c3 = topo.child3_id[i]

        c1 > 0 && (r_gamma_sum += seg.radius[c1]^gamma)
        c2 > 0 && (r_gamma_sum += seg.radius[c2]^gamma)
        c3 > 0 && (r_gamma_sum += seg.radius[c3]^gamma)

        r_gamma_sum > 0.0 && (seg.radius[i] = r_gamma_sum^(1.0 / gamma))
    end
    return nothing
end
