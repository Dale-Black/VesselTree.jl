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

    # Find max generation
    max_gen = Int32(0)
    for i in 1:n
        if topo.generation[i] > max_gen
            max_gen = topo.generation[i]
        end
    end

    # Bottom-up: process each generation from deepest non-terminal to root
    for gen in (max_gen - Int32(1)):-Int32(1):Int32(0)
        for i in 1:n
            if topo.generation[i] == gen && !topo.is_terminal[i]
                # Sum of children's r^gamma
                r_gamma_sum = 0.0
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                c3 = topo.child3_id[i]

                if c1 > 0
                    r_gamma_sum += seg.radius[c1]^gamma
                end
                if c2 > 0
                    r_gamma_sum += seg.radius[c2]^gamma
                end
                if c3 > 0
                    r_gamma_sum += seg.radius[c3]^gamma
                end

                seg.radius[i] = r_gamma_sum^(1.0 / gamma)
            end
        end
    end
    return nothing
end
