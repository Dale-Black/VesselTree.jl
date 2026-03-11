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

    # Save original children of parent (if it was already a bifurcation/trifurcation)
    orig_child1 = topo.child1_id[parent_idx]
    orig_child2 = topo.child2_id[parent_idx]
    orig_child3 = topo.child3_id[parent_idx]

    # Reset parent children (they will be set by add_segment!)
    topo.child1_id[parent_idx] = Int32(-1)
    topo.child2_id[parent_idx] = Int32(-1)
    topo.child3_id[parent_idx] = Int32(-1)
    topo.is_terminal[parent_idx] = true
    if topo.junction_type[parent_idx] == :bifurcation
        tree.n_bifurcations -= 1
        tree.n_terminals += 1
    elseif topo.junction_type[parent_idx] == :trifurcation
        tree.n_trifurcations -= 1
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
    if orig_child3 > 0
        topo.child3_id[cont_id] = orig_child3
        topo.parent_id[orig_child3] = cont_id
        topo.junction_type[cont_id] = :trifurcation
        tree.n_bifurcations -= 1
        tree.n_trifurcations += 1
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
    _domain_size(domain) -> Float64

Compute representative size of a domain for adaptive thresholding.
"""
function _domain_size(domain::AbstractDomain)
    if domain isa SphereDomain
        return domain.radius
    elseif domain isa BoxDomain
        return max(
            domain.max_corner[1] - domain.min_corner[1],
            domain.max_corner[2] - domain.min_corner[2],
            domain.max_corner[3] - domain.min_corner[3],
        )
    elseif domain isa EllipsoidDomain
        return max(domain.semi_axes[1], domain.semi_axes[2], domain.semi_axes[3])
    elseif domain isa EllipsoidShellDomain
        return max(domain.semi_axes[1], domain.semi_axes[2], domain.semi_axes[3])
    elseif domain isa CSVVolumeDomain
        return max(
            domain.max_corner[1] - domain.min_corner[1],
            domain.max_corner[2] - domain.min_corner[2],
            domain.max_corner[3] - domain.min_corner[3],
        )
    end
    return 1.0
end

"""
    _grid_cell_size(domain_size, n_segments) -> Float64

Adaptive cell size: shrinks as tree grows to maintain good locality.
"""
function _grid_cell_size(domain_size::Float64, n_segments::Int)
    # Target ~10-20 segments per cell on average
    # cell volume ~ domain^3 / n_cells, n_cells ~ n_segments / 10
    return domain_size / max(1.0, (n_segments / 10.0)^(1.0 / 3.0))
end

"""
    _insert_segment_to_grid!(grid, segments, idx)

Insert segment `idx` into spatial grid by its midpoint.
"""
function _insert_segment_to_grid!(grid::SpatialGrid, segments::SegmentData, idx::Int)
    mx = (segments.proximal_x[idx] + segments.distal_x[idx]) / 2.0
    my = (segments.proximal_y[idx] + segments.distal_y[idx]) / 2.0
    mz = (segments.proximal_z[idx] + segments.distal_z[idx]) / 2.0
    insert!(grid, idx, mx, my, mz)
end

# Threshold: only use grid when tree has enough segments for it to matter.
# Below this, brute force is fast and gives exact K-nearest.
const _GRID_ACTIVATION_THRESHOLD = 200

"""
    _find_nearest_via_grid(grid, segments, tx, ty, tz, search_radius, K, n) -> Vector{Int}

Find K nearest segments to point (tx,ty,tz) using spatial grid. Falls back to
brute-force scan if grid returns fewer than K candidates.
"""
function _find_nearest_via_grid(
    grid::SpatialGrid,
    segments::SegmentData,
    distances_buf::Vector{Float64},
    tx::Float64, ty::Float64, tz::Float64,
    search_radius::Float64,
    K::Int,
    n::Int,
)
    # For small trees, brute force is fast and exact
    if n < _GRID_ACTIVATION_THRESHOLD
        compute_all_distances!(distances_buf, segments, tx, ty, tz, n)
        K_actual = min(K, n)
        return collect(partialsortperm(@view(distances_buf[1:n]), 1:K_actual))
    end

    nearby = query_nearby(grid, tx, ty, tz, search_radius)

    if length(nearby) >= K
        # Compute distances only for nearby segments
        nearby_dists = Vector{Float64}(undef, length(nearby))
        for (j, si) in enumerate(nearby)  # sequential: small subset from grid query
            nearby_dists[j] = point_segment_distance(
                tx, ty, tz,
                segments.proximal_x[si], segments.proximal_y[si], segments.proximal_z[si],
                segments.distal_x[si], segments.distal_y[si], segments.distal_z[si],
            )
        end
        k_actual = min(K, length(nearby))
        perm = partialsortperm(nearby_dists, 1:k_actual)
        return [nearby[perm[j]] for j in 1:k_actual]
    else
        # Fall back to full scan
        compute_all_distances!(distances_buf, segments, tx, ty, tz, n)
        K_actual = min(K, n)
        return collect(partialsortperm(@view(distances_buf[1:n]), 1:K_actual))
    end
end

"""
    grow_tree!(tree, domain, n_terminals, params; rng, verbose, kassab, trifurcation, barabasi)

Main CCO growth loop. Adds `n_terminals` terminal segments to the tree.
Uses SpatialGrid for O(1) nearest-neighbor queries instead of O(n) full scan.

Keywords:
- `kassab=false`: Kassab asymmetry sampling for daughter radii
- `trifurcation=false`: Trifurcation merges when chi > 0.83
- `barabasi=false`: Apply Barabasi junction geometry (sprouting/branching angles)

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
    trifurcation::Bool=false,
    barabasi::Bool=false,
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

    domain_size = _domain_size(domain)

    # Build initial spatial grid
    n_initial = tree.segments.n
    cell_size = _grid_cell_size(domain_size, max(n_initial, 10))
    grid = build_grid(tree.segments, n_initial, domain, cell_size)
    last_rebuild_n = n_initial
    grid_rebuild_interval = 100

    added = 0
    attempts = 0
    max_total_attempts = n_terminals * 100

    while added < n_terminals && attempts < max_total_attempts
        attempts += 1
        n = tree.segments.n

        # Rebuild grid periodically
        if n - last_rebuild_n >= grid_rebuild_interval
            cell_size = _grid_cell_size(domain_size, n)
            grid = build_grid(tree.segments, n, domain, cell_size)
            last_rebuild_n = n
        end

        # OpenCCO Eq. 15: d_thresh = (V / (k * (i+1)))^(1/3)
        # Uses domain volume and target terminal count for principled spacing
        domain_vol = domain_size^3  # approximate; for shells use _effective_domain_volume
        d_thresh = (domain_vol / (n_terminals * max(1, added + 1)))^(1.0 / 3.0)

        # Progressive expansion: search radius grows from near-root to full domain
        # Must be >= d_thresh * 2 so the annulus [d_thresh, search_radius] is non-empty
        growth_frac = min(1.0, (added + 1.0) / n_terminals)
        search_radius = max(domain_size * (0.1 + 0.9 * cbrt(growth_frac)), d_thresh * 2.0)

        # Sample candidate point (checks min distance via d_thresh internally)
        candidate = sample_terminal_candidate(domain, tree, d_thresh, distances_buf, rng)
        if candidate === nothing
            continue
        end

        tx, ty, tz = candidate

        # Find best connection among K nearest segments using spatial grid
        K = min(n, 40)
        nearest = _find_nearest_via_grid(grid, tree.segments, distances_buf, tx, ty, tz, search_radius, K, n)

        seg_idx = 0
        t_opt = 0.5
        best_cost = Inf
        found = false

        for ki in 1:length(nearest)
            ci = nearest[ki]  # candidate segment index

            # Optimize bifurcation point on this segment
            ti, costi = optimize_bifurcation_point(
                tree.segments.proximal_x[ci], tree.segments.proximal_y[ci], tree.segments.proximal_z[ci],
                tree.segments.distal_x[ci], tree.segments.distal_y[ci], tree.segments.distal_z[ci],
                tree.segments.radius[ci], tx, ty, tz, gamma)

            costi >= best_cost && continue

            # Compute bifurcation point
            bx = tree.segments.proximal_x[ci] + ti * (tree.segments.distal_x[ci] - tree.segments.proximal_x[ci])
            by = tree.segments.proximal_y[ci] + ti * (tree.segments.distal_y[ci] - tree.segments.proximal_y[ci])
            bz = tree.segments.proximal_z[ci] + ti * (tree.segments.distal_z[ci] - tree.segments.proximal_z[ci])

            # Check intersection
            parent_r = tree.segments.radius[ci]
            intersect_dist = max(parent_r * 2.0, d_thresh * 0.01)
            check_intersections!(intersect_buf, tree.segments, bx, by, bz, tx, ty, tz, intersect_dist, n)

            # Exclude topological neighbors
            intersect_buf[ci] = false
            topo = tree.topology
            pp = topo.parent_id[ci]
            if pp > 0
                intersect_buf[pp] = false
                for c in (topo.child1_id[pp], topo.child2_id[pp], topo.child3_id[pp])
                    c > 0 && (intersect_buf[c] = false)
                end
            end
            for c in (topo.child1_id[ci], topo.child2_id[ci], topo.child3_id[ci])
                c > 0 && (intersect_buf[c] = false)
            end

            has_any_intersection(intersect_buf, n) && continue

            # Check domain crossing
            check_domain_crossing(domain, (bx, by, bz), (tx, ty, tz)) && continue

            best_cost = costi
            seg_idx = ci
            t_opt = ti
            found = true
        end

        !found && continue

        # Bifurcation point (recompute for selected segment)
        bx = tree.segments.proximal_x[seg_idx] + t_opt * (tree.segments.distal_x[seg_idx] - tree.segments.proximal_x[seg_idx])
        by = tree.segments.proximal_y[seg_idx] + t_opt * (tree.segments.distal_y[seg_idx] - tree.segments.proximal_y[seg_idx])
        bz = tree.segments.proximal_z[seg_idx] + t_opt * (tree.segments.distal_z[seg_idx] - tree.segments.proximal_z[seg_idx])

        # Check for trifurcation merge before creating bifurcation
        if trifurcation
            merge_target = check_trifurcation_merge(tree, (bx, by, bz), params)
            if merge_target !== nothing
                # Merge into existing bifurcation as trifurcation
                new_id = merge_to_trifurcation!(tree, merge_target, (tx, ty, tz), terminal_radius, params)

                if kassab
                    parent_radius = tree.segments.radius[Int(merge_target)]
                    asymmetry = sample_asymmetry(params, rng)
                    _, r_small = compute_daughter_radii(parent_radius, asymmetry, gamma)
                    tree.segments.radius[Int(new_id)] = r_small
                end

                update_radii!(tree, gamma)

                # Insert new segment into grid
                _insert_segment_to_grid!(grid, tree.segments, Int(new_id))

                added += 1

                if verbose && added % 100 == 0
                    println("  grow_tree!: $added / $n_terminals terminals added ($(tree.segments.n) segments)")
                end
                continue
            end
        end

        n_before = tree.segments.n

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

        # Apply Barabasi junction geometry
        if barabasi
            apply_junction_geometry!(tree, Int32(seg_idx), params)
        end

        # Insert new segments into grid (parent was modified, add new daughters)
        _insert_segment_to_grid!(grid, tree.segments, seg_idx)  # re-insert shortened parent
        for new_seg_idx in (n_before + 1):tree.segments.n
            _insert_segment_to_grid!(grid, tree.segments, new_seg_idx)
        end

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
