# Core CCO growth loop

"""
    add_bifurcation!(tree, parent_idx, bifurc_t, terminal_point, terminal_radius)

Split segment `parent_idx` at parametric position `bifurc_t`, creating a bifurcation
with two new daughter segments: one continuing to the original distal point,
one going to the new terminal point.

Returns (cont_id, new_id) — IDs of the continuation and new terminal segments.
"""
function add_bifurcation!(
    tree::VascularTree,
    parent_idx::Int,
    bifurc_t::Float64,
    terminal_point::NTuple{3,Float64},
    terminal_radius::Float64,
)
    seg = tree.segments
    topo = tree.topology

    # Original segment endpoints
    px = seg.proximal_x[parent_idx]
    py = seg.proximal_y[parent_idx]
    pz = seg.proximal_z[parent_idx]
    dx = seg.distal_x[parent_idx]
    dy = seg.distal_y[parent_idx]
    dz = seg.distal_z[parent_idx]
    orig_radius = seg.radius[parent_idx]

    # Bifurcation point
    bx = px + bifurc_t * (dx - px)
    by = py + bifurc_t * (dy - py)
    bz = pz + bifurc_t * (dz - pz)

    # Modify parent segment: shorten to proximal → bifurcation point
    seg.distal_x[parent_idx] = bx
    seg.distal_y[parent_idx] = by
    seg.distal_z[parent_idx] = bz
    ex = bx - px; ey = by - py; ez = bz - pz
    seg.seg_length[parent_idx] = sqrt(ex*ex + ey*ey + ez*ez)

    # Save original children of parent (if it was already a bifurcation)
    orig_child1 = topo.child1_id[parent_idx]
    orig_child2 = topo.child2_id[parent_idx]

    # Reset parent children (they will be set by add_segment!)
    topo.child1_id[parent_idx] = Int32(-1)
    topo.child2_id[parent_idx] = Int32(-1)
    topo.is_terminal[parent_idx] = true
    if topo.junction_type[parent_idx] == :bifurcation
        tree.n_bifurcations -= 1
        tree.n_terminals += 1
    end
    topo.junction_type[parent_idx] = :none

    # Continuation segment: bifurcation → original distal
    cont_id = add_segment!(tree, (bx, by, bz), (dx, dy, dz), orig_radius, Int32(parent_idx))

    # Transfer original children to continuation segment
    if orig_child1 > 0
        topo.child1_id[cont_id] = orig_child1
        topo.parent_id[orig_child1] = cont_id
        topo.is_terminal[cont_id] = false
        topo.junction_type[cont_id] = :bifurcation
        tree.n_bifurcations += 1
        tree.n_terminals -= 1  # cont_id is no longer terminal
    end
    if orig_child2 > 0
        topo.child2_id[cont_id] = orig_child2
        topo.parent_id[orig_child2] = cont_id
    end

    # New terminal segment: bifurcation → terminal point
    new_id = add_segment!(tree, (bx, by, bz), terminal_point, terminal_radius, Int32(parent_idx))

    return Int(cont_id), Int(new_id)
end

"""
    sample_terminal_candidate(domain, tree, params, distances_buf, rng; max_attempts=100)

Sample a candidate terminal point in the domain that is not too close to existing segments.
Returns the point as (x, y, z) or `nothing` if no valid point found.
"""
function sample_terminal_candidate(
    domain::AbstractDomain,
    tree::VascularTree,
    min_distance::Float64,
    distances_buf::Vector{Float64},
    rng::AbstractRNG;
    max_attempts::Int=100,
)
    n = tree.segments.n

    for _ in 1:max_attempts
        p = sample_point(domain, rng)

        if n == 0
            return p
        end

        # Check distance to all existing segments
        compute_all_distances!(distances_buf, tree.segments, p[1], p[2], p[3], n)
        min_dist = AK.minimum(@view distances_buf[1:n])

        if min_dist >= min_distance
            return p
        end
    end
    return nothing
end

"""
    grow_tree!(tree, domain, n_terminals, params; rng=Random.default_rng(), verbose=false, kassab=false)

Main CCO growth loop. Adds `n_terminals` terminal segments to the tree.

When `kassab=true`, uses Kassab asymmetry sampling for daughter radii instead
of symmetric bifurcations.

The tree must have a root segment already added before calling this function.
"""
function grow_tree!(
    tree::VascularTree,
    domain::AbstractDomain,
    n_terminals::Int,
    params::MorphometricParams;
    rng::AbstractRNG=Random.default_rng(),
    verbose::Bool=false,
    kassab::Bool=false,
)
    capacity = tree.segments.capacity
    gamma = params.gamma

    # Work buffers (pre-allocated)
    distances_buf = zeros(capacity)
    costs_buf = zeros(capacity)
    bifurc_t_buf = zeros(capacity)
    intersect_buf = Vector{Bool}(undef, capacity)

    # Initial terminal radius estimate
    terminal_radius = params.diameter_mean[1] / 2.0 / 1000.0  # um → mm (or whatever units)

    # Minimum distance between segments (scales with tree size)
    domain_size = 1.0
    if domain isa SphereDomain
        domain_size = domain.radius
    elseif domain isa BoxDomain
        domain_size = max(
            domain.max_corner[1] - domain.min_corner[1],
            domain.max_corner[2] - domain.min_corner[2],
            domain.max_corner[3] - domain.min_corner[3],
        )
    elseif domain isa EllipsoidDomain
        domain_size = max(domain.semi_axes[1], domain.semi_axes[2], domain.semi_axes[3])
    end

    added = 0
    attempts = 0
    max_total_attempts = n_terminals * 100

    while added < n_terminals && attempts < max_total_attempts
        attempts += 1
        n = tree.segments.n

        # Adaptive minimum distance
        d_thresh = domain_size / (10.0 * (n + 1)^(1.0 / 3.0))

        # Sample candidate point
        candidate = sample_terminal_candidate(domain, tree, d_thresh, distances_buf, rng)
        if candidate === nothing
            continue
        end

        tx, ty, tz = candidate

        # Find best connection via Kamiya optimization
        evaluate_all_connections!(costs_buf, bifurc_t_buf, tree.segments, n, tx, ty, tz, gamma)
        seg_idx, t_opt, cost = select_best_connection(costs_buf, bifurc_t_buf, n)

        # Bifurcation point
        bx = tree.segments.proximal_x[seg_idx] + t_opt * (tree.segments.distal_x[seg_idx] - tree.segments.proximal_x[seg_idx])
        by = tree.segments.proximal_y[seg_idx] + t_opt * (tree.segments.distal_y[seg_idx] - tree.segments.proximal_y[seg_idx])
        bz = tree.segments.proximal_z[seg_idx] + t_opt * (tree.segments.distal_z[seg_idx] - tree.segments.proximal_z[seg_idx])

        # Check intersection: new segment (bifurcation → terminal)
        check_intersections!(intersect_buf, tree.segments, bx, by, bz, tx, ty, tz, d_thresh * 0.1, n)
        # Exclude the parent segment from intersection check
        intersect_buf[seg_idx] = false
        if has_any_intersection(intersect_buf, n)
            continue
        end

        # Check domain crossing
        if check_domain_crossing(domain, (bx, by, bz), (tx, ty, tz))
            continue
        end

        # Add bifurcation
        cont_id, new_id = add_bifurcation!(tree, seg_idx, t_opt, (tx, ty, tz), terminal_radius)

        if kassab
            # Kassab-constrained: sample asymmetry and set daughter radii
            parent_radius = tree.segments.radius[seg_idx]
            asymmetry = sample_asymmetry(params, rng)
            r_large, r_small = compute_daughter_radii(parent_radius, asymmetry, gamma)
            # Continuation segment gets larger radius, new terminal gets smaller
            tree.segments.radius[cont_id] = r_large
            tree.segments.radius[new_id] = r_small
        end

        # Update radii via Murray's law
        update_radii!(tree, gamma)

        added += 1

        if verbose && added % 100 == 0
            println("  grow_tree!: $added / $n_terminals terminals added ($(tree.segments.n) segments)")
        end
    end

    if verbose
        println("  grow_tree! complete: $added terminals, $(tree.segments.n) segments, $attempts attempts")
    end

    return tree
end
