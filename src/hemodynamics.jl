# Hemodynamic computations: Poiseuille resistance, flow, pressure

"""
    compute_resistances!(tree::VascularTree, viscosity::Float64)

Compute Poiseuille resistance for all active segments using AK kernel:
  R[i] = 8 * mu * L[i] / (pi * r[i]^4)
"""
function compute_resistances!(tree::VascularTree, viscosity::Float64)
    seg = tree.segments
    n = seg.n
    n == 0 && return nothing

    len_view = @view seg.seg_length[1:n]
    rad_view = @view seg.radius[1:n]
    res_view = @view seg.resistance[1:n]

    coeff = 8.0 * viscosity / π

    AK.foreachindex(res_view) do i
        res_view[i] = coeff * len_view[i] / (rad_view[i]^4)
    end
    return nothing
end

"""
    _topo_order(tree::VascularTree) -> Vector{Int}

BFS from root, returning segment IDs in topological order (parent before children).
"""
function _topo_order(tree::VascularTree)
    # sequential: BFS tree traversal
    topo = tree.topology
    n = tree.segments.n
    order = Int[]
    sizehint!(order, n)

    root_id = Int(tree.root_segment_id)
    root_id < 1 && return order

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
    return order
end

"""
    compute_flows!(tree::VascularTree, params::MorphometricParams)

Compute flows through all segments. Uses resistance-weighted flow splitting
and equivalent resistance computation.

Flow is distributed so that Q_root = (P_root - P_term) / R_total_tree,
and at each bifurcation flow splits inversely proportional to subtree resistance.
"""
function compute_flows!(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    n == 0 && return nothing

    order = _topo_order(tree)

    # Step 1: Bottom-up — compute subtree resistance
    subtree_R = zeros(n)
    for k in length(order):-1:1
        i = order[k]
        if topo.is_terminal[i]
            subtree_R[i] = seg.resistance[i]
        else
            conductance_sum = 0.0
            c1 = Int(topo.child1_id[i])
            c2 = Int(topo.child2_id[i])
            c3 = Int(topo.child3_id[i])
            c1 > 0 && (conductance_sum += 1.0 / subtree_R[c1])
            c2 > 0 && (conductance_sum += 1.0 / subtree_R[c2])
            c3 > 0 && (conductance_sum += 1.0 / subtree_R[c3])
            subtree_R[i] = seg.resistance[i] + 1.0 / conductance_sum
        end
    end

    # Step 2: Total flow from root
    root_id = order[1]
    total_flow = (params.root_pressure - params.terminal_pressure) / subtree_R[root_id]
    seg.flow[root_id] = total_flow

    # Step 3: Top-down — distribute flow
    # sequential: tree traversal — parent processed before children
    for k in 1:length(order)
        i = order[k]
        topo.is_terminal[i] && continue

        parent_flow = seg.flow[i]
        c1 = Int(topo.child1_id[i])
        c2 = Int(topo.child2_id[i])
        c3 = Int(topo.child3_id[i])

        conductance_sum = 0.0
        c1 > 0 && (conductance_sum += 1.0 / subtree_R[c1])
        c2 > 0 && (conductance_sum += 1.0 / subtree_R[c2])
        c3 > 0 && (conductance_sum += 1.0 / subtree_R[c3])

        if c1 > 0
            seg.flow[c1] = parent_flow * (1.0 / subtree_R[c1]) / conductance_sum
        end
        if c2 > 0
            seg.flow[c2] = parent_flow * (1.0 / subtree_R[c2]) / conductance_sum
        end
        if c3 > 0
            seg.flow[c3] = parent_flow * (1.0 / subtree_R[c3]) / conductance_sum
        end
    end

    return nothing
end

"""
    compute_pressures!(tree::VascularTree, params::MorphometricParams)

Compute pressure at proximal and distal ends of each segment:
  P_distal = P_proximal - Q * R

Root proximal pressure is set to params.root_pressure.
"""
function compute_pressures!(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    n == 0 && return nothing

    order = _topo_order(tree)

    # sequential: tree traversal — parent processed before children
    root_id = order[1]
    seg.pressure_proximal[root_id] = params.root_pressure
    seg.pressure_distal[root_id] = params.root_pressure - seg.flow[root_id] * seg.resistance[root_id]

    for k in 2:length(order)
        i = order[k]
        parent = Int(topo.parent_id[i])
        seg.pressure_proximal[i] = seg.pressure_distal[parent]
        seg.pressure_distal[i] = seg.pressure_proximal[i] - seg.flow[i] * seg.resistance[i]
    end

    return nothing
end

"""
    validate_hemodynamics(tree::VascularTree, params::MorphometricParams) -> Bool

Check hemodynamic consistency:
- All flows positive
- All pressures non-negative
- Flow conservation at bifurcations
- Pressure monotonically decreases along each segment
"""
function validate_hemodynamics(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    n == 0 && return true

    for i in 1:n
        # All flows must be positive
        seg.flow[i] <= 0.0 && return false

        # All pressures must be non-negative
        seg.pressure_proximal[i] < 0.0 && return false
        seg.pressure_distal[i] < 0.0 && return false

        # Pressure must decrease along segment
        seg.pressure_proximal[i] < seg.pressure_distal[i] - 1e-6 && return false

        # Flow conservation at bifurcations
        if topo.junction_type[i] == :bifurcation
            c1 = topo.child1_id[i]
            c2 = topo.child2_id[i]
            if c1 > 0 && c2 > 0
                flow_sum = seg.flow[c1] + seg.flow[c2]
                if !isapprox(seg.flow[i], flow_sum, rtol=1e-4)
                    return false
                end
            end
        elseif topo.junction_type[i] == :trifurcation
            c1 = topo.child1_id[i]
            c2 = topo.child2_id[i]
            c3 = topo.child3_id[i]
            if c1 > 0 && c2 > 0 && c3 > 0
                flow_sum = seg.flow[c1] + seg.flow[c2] + seg.flow[c3]
                if !isapprox(seg.flow[i], flow_sum, rtol=1e-4)
                    return false
                end
            end
        end
    end

    return true
end
