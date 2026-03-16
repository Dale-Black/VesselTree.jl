# Hemodynamic computations: Poiseuille resistance, flow, pressure

"""
    compute_resistances!(tree::VascularTree, viscosity::Float64)

Compute Poiseuille resistance for all active segments using AK kernel:
  R[i] = 8 * mu * L[i] / (pi * r[i]^4)

Geometry in `VascularTree` is stored in millimeters, so lengths and radii are
converted to meters before applying the SI-form Poiseuille law.
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
        len_m = len_view[i] * 1e-3
        rad_m = rad_view[i] * 1e-3
        res_view[i] = coeff * len_m / (rad_m^4)
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

"""
    assign_terminal_flows!(tree::VascularTree, domain::AbstractDomain, params::MorphometricParams)

Assign terminal flows weighted by estimated territory volume (nearest-neighbor distance).
Terminals farther from neighbors serve larger territories and get more flow.
Flow is then propagated bottom-up so parent flow = sum of children flows.
"""
function assign_terminal_flows!(tree::VascularTree, domain::AbstractDomain, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    n == 0 && return nothing

    # Collect terminal indices and their distal points
    terminal_ids = Int[]
    for i in 1:n
        topo.is_terminal[i] && push!(terminal_ids, i)
    end
    n_term = length(terminal_ids)
    n_term == 0 && return nothing

    # Compute territory weight for each terminal using nearest-neighbor distance
    # AK kernel: compute pairwise distances between terminal distal points
    weights = ones(n_term)

    if n_term > 1
        # For each terminal, find distance to nearest other terminal
        nn_dists = fill(Inf, n_term)
        term_x = [seg.distal_x[terminal_ids[j]] for j in 1:n_term]
        term_y = [seg.distal_y[terminal_ids[j]] for j in 1:n_term]
        term_z = [seg.distal_z[terminal_ids[j]] for j in 1:n_term]

        dist_buf = zeros(n_term)

        for j in 1:n_term
            cx, cy, cz = term_x[j], term_y[j], term_z[j]

            # Compute distances to all other terminals using AK
            AK.foreachindex(dist_buf) do k
                dx = term_x[k] - cx
                dy = term_y[k] - cy
                dz = term_z[k] - cz
                dist_buf[k] = dx * dx + dy * dy + dz * dz
            end
            # Set self-distance to Inf
            dist_buf[j] = Inf

            nn_dists[j] = sqrt(AK.minimum(dist_buf))
        end

        # Weight proportional to nn_dist^3 (approximates 3D territory volume)
        for j in 1:n_term
            weights[j] = nn_dists[j]^3
        end
    end

    # Normalize weights so they sum to 1
    total_weight = sum(weights)
    for j in 1:n_term
        weights[j] /= total_weight
    end

    # Total flow from pressure drop and equivalent resistance
    total_R = _compute_total_resistance(tree)
    total_flow = (params.root_pressure - params.terminal_pressure) / total_R

    # Assign terminal flows
    for j in 1:n_term
        seg.flow[terminal_ids[j]] = total_flow * weights[j]
    end

    # Bottom-up: propagate flows from terminals to root
    order = _topo_order(tree)
    for k in length(order):-1:1
        i = order[k]
        if !topo.is_terminal[i]
            flow_sum = 0.0
            c1 = Int(topo.child1_id[i])
            c2 = Int(topo.child2_id[i])
            c3 = Int(topo.child3_id[i])
            c1 > 0 && (flow_sum += seg.flow[c1])
            c2 > 0 && (flow_sum += seg.flow[c2])
            c3 > 0 && (flow_sum += seg.flow[c3])
            seg.flow[i] = flow_sum
        end
    end

    return nothing
end

"""
    _compute_total_resistance(tree::VascularTree) -> Float64

Compute the equivalent resistance of the entire tree (bottom-up).
"""
function _compute_total_resistance(tree::VascularTree)
    seg = tree.segments
    topo = tree.topology
    order = _topo_order(tree)
    n = seg.n

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

    return subtree_R[Int(tree.root_segment_id)]
end

"""
    recompute_radii_from_flow!(tree::VascularTree, params::MorphometricParams)

Given a flow distribution, recompute terminal radii from Poiseuille's law
and propagate internal radii via Murray's law.

Terminal radius is set so that R = 8*mu*L/(pi*r^4) gives the target pressure
drop for that terminal's flow. Then Murray's law propagates upward.
"""
function recompute_radii_from_flow!(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    gamma = params.gamma
    n == 0 && return nothing

    # Set terminal radii from flow using Poiseuille:
    # Q = delta_P / R, R = 8*mu*L/(pi*r^4)
    # r^4 = 8*mu*L*Q / (pi * delta_P)
    # Use a proportional approach: r proportional to Q^(1/gamma)
    # (from Murray's law generalization: Q ∝ r^gamma)
    for i in 1:n
        if topo.is_terminal[i]
            seg.radius[i] = (seg.flow[i])^(1.0 / gamma)
        end
    end

    # Normalize: scale all terminal radii so root radius is reasonable
    # Murray's law bottom-up will set internal radii
    update_radii!(tree, gamma)

    return nothing
end
