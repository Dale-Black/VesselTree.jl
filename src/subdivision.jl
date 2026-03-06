# Statistical subdivision engine — Kassab connectivity matrix-based
# Recursively subdivides CCO skeleton terminals into full microvasculature

"""
    estimate_total_segments(parent_order, params) -> Int

Estimate total segments produced by recursively subdividing a parent of given
Strahler order using the connectivity matrix. Uses expected values (no randomness).
"""
function estimate_total_segments(parent_order::Int, params::MorphometricParams)
    parent_order <= 0 && return 1
    CM = params.connectivity_matrix
    n_orders = params.n_orders
    total = 1.0  # the parent itself
    for m in 0:(parent_order - 1)
        if m + 1 <= n_orders && parent_order + 1 <= n_orders
            avg = CM[m + 1, parent_order + 1]
            avg > 0 && (total += avg * estimate_total_segments(m, params))
        end
    end
    return ceil(Int, total)
end

"""
    estimate_subdivision_capacity(tree, params) -> Int

Estimate total capacity needed after subdividing all current terminals.
"""
function estimate_subdivision_capacity(tree::VascularTree, params::MorphometricParams)
    assign_strahler_orders!(tree, params)
    total = tree.segments.n  # existing segments
    for i in 1:tree.segments.n
        if tree.topology.is_terminal[i]
            ord = Int(tree.topology.strahler_order[i])
            if ord > 0
                total += estimate_total_segments(ord, params) - 1  # subtract parent (already counted)
            end
        end
    end
    # Add 50% buffer for Poisson variance
    return ceil(Int, total * 1.5)
end

"""
    _random_daughter_direction(parent_dir, is_first, d_order, parent_order, rng)

Compute a direction for a daughter segment based on parent direction.
First daughter (continuation) gets small perturbation. Branches get larger offsets.
"""
function _random_daughter_direction(
    parent_dir::NTuple{3,Float64},
    is_first::Bool,
    d_order::Int,
    parent_order::Int,
    rng::AbstractRNG,
)
    dx, dy, dz = parent_dir

    if is_first
        # Continuation: small perturbation (< 15 degrees)
        max_angle = 0.26  # ~15 degrees
    else
        # Branch: angle depends on order ratio (Barabasi-like)
        rho = parent_order > 0 ? d_order / parent_order : 0.0
        if rho < 0.83
            # Sprouting: near-perpendicular
            max_angle = Float64(π) / 2.0 * (0.7 + 0.6 * rand(rng))
        else
            # Branching: moderate angle
            max_angle = Float64(π) / 4.0 * (0.7 + 0.6 * rand(rng))
        end
    end

    # Find perpendicular basis
    perp_x, perp_y, perp_z = _find_perpendicular(dx, dy, dz)
    # Second perpendicular via cross product
    perp2_x = dy * perp_z - dz * perp_y
    perp2_y = dz * perp_x - dx * perp_z
    perp2_z = dx * perp_y - dy * perp_x

    # Random azimuthal angle
    phi = 2.0 * Float64(π) * rand(rng)

    # Tilt direction by max_angle in random azimuthal plane
    ca = cos(max_angle)
    sa = sin(max_angle)
    cp = cos(phi)
    sp = sin(phi)

    new_dx = dx * ca + (perp_x * cp + perp2_x * sp) * sa
    new_dy = dy * ca + (perp_y * cp + perp2_y * sp) * sa
    new_dz = dz * ca + (perp_z * cp + perp2_z * sp) * sa

    # Normalize
    len = sqrt(new_dx^2 + new_dy^2 + new_dz^2)
    if len > 1e-15
        new_dx /= len; new_dy /= len; new_dz /= len
    end

    return (new_dx, new_dy, new_dz)
end

"""
    subdivide_terminals!(tree, params; rng)

Recursively subdivide each terminal of order > 0 down to order 0 (capillaries)
using Kassab's connectivity matrix. This is the statistical fast-path for
generating millions of microvasculature segments without intersection checking.

Algorithm per terminal of order k:
1. Look up CM column k for expected daughter counts per order
2. Sample actual counts from Poisson(CM[m+1, k+1])
3. Create daughter segments with sampled diameter/length
4. Recursively subdivide each daughter of order > 0
"""
function subdivide_terminals!(
    tree::VascularTree,
    params::MorphometricParams;
    rng::AbstractRNG=Random.default_rng(),
    max_order::Int=params.n_orders - 1,
)
    assign_strahler_orders!(tree, params)

    # Snapshot current terminals (don't iterate over segments we'll add)
    n = tree.segments.n
    term_indices = Int[]
    term_orders = Int[]
    for i in 1:n  # sequential: collect terminal list
        if tree.topology.is_terminal[i]
            ord = Int(tree.topology.strahler_order[i])
            if ord > 0 && ord <= max_order
                push!(term_indices, i)
                push!(term_orders, ord)
            end
        end
    end

    for k in 1:length(term_indices)
        _subdivide_recursive!(tree, term_indices[k], term_orders[k], params, rng)
    end

    return tree
end

"""
    _subdivide_recursive!(tree, parent_idx, parent_order, params, rng)

Recursively subdivide a single segment into daughters per Kassab CM.
"""
function _subdivide_recursive!(
    tree::VascularTree,
    parent_idx::Int,
    parent_order::Int,
    params::MorphometricParams,
    rng::AbstractRNG,
)
    parent_order <= 0 && return

    CM = params.connectivity_matrix
    gamma = params.gamma
    n_orders = params.n_orders
    seg = tree.segments

    # Parent direction
    dir_x = seg.distal_x[parent_idx] - seg.proximal_x[parent_idx]
    dir_y = seg.distal_y[parent_idx] - seg.proximal_y[parent_idx]
    dir_z = seg.distal_z[parent_idx] - seg.proximal_z[parent_idx]
    dir_len = sqrt(dir_x^2 + dir_y^2 + dir_z^2)
    if dir_len > 1e-15
        dir_x /= dir_len; dir_y /= dir_len; dir_z /= dir_len
    else
        dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
    end
    parent_dir = (dir_x, dir_y, dir_z)

    # Build daughter list from CM column
    daughters = Int[]  # list of daughter orders
    for m in 0:(parent_order - 1)
        if m + 1 <= n_orders && parent_order + 1 <= n_orders
            lambda = CM[m + 1, parent_order + 1]
            lambda <= 0.0 && continue
            n_d = rand(rng, Poisson(lambda))
            for _ in 1:n_d
                push!(daughters, m)
            end
        end
    end

    isempty(daughters) && return
    sort!(daughters, rev=true)  # highest order first (continuation)

    # Add daughters with cascading if > available child slots
    daughter_queue = copy(daughters)
    current_parent = parent_idx
    daughter_idx = 0  # tracks position for is_first logic

    while !isempty(daughter_queue)
        topo = tree.topology
        # Count existing children on current_parent
        n_existing = 0
        topo.child1_id[current_parent] > 0 && (n_existing += 1)
        topo.child2_id[current_parent] > 0 && (n_existing += 1)
        topo.child3_id[current_parent] > 0 && (n_existing += 1)
        available = 3 - n_existing

        available <= 0 && break  # shouldn't happen if cascading works

        if length(daughter_queue) <= available
            # All remaining fit directly
            for d_order in daughter_queue
                daughter_idx += 1
                _add_subdivision_daughter!(tree, current_parent, d_order, parent_dir, daughter_idx == 1, parent_order, params, rng)
            end
            break
        end

        # Need cascading: add (available - 1) daughters + 1 stub
        n_to_add = max(available - 1, 0)
        for _ in 1:n_to_add
            d_order = popfirst!(daughter_queue)
            daughter_idx += 1
            _add_subdivision_daughter!(tree, current_parent, d_order, parent_dir, daughter_idx == 1, parent_order, params, rng)
        end

        # Create cascade stub in last slot
        stub_len = params.vessel_cutoff_um / 1000.0 * 0.01
        cp_x = seg.distal_x[current_parent]
        cp_y = seg.distal_y[current_parent]
        cp_z = seg.distal_z[current_parent]
        stub_distal = (
            cp_x + parent_dir[1] * stub_len,
            cp_y + parent_dir[2] * stub_len,
            cp_z + parent_dir[3] * stub_len,
        )
        stub_radius = seg.radius[current_parent]
        stub_id = add_segment!(tree, (cp_x, cp_y, cp_z), stub_distal, stub_radius, Int32(current_parent))
        current_parent = Int(stub_id)
    end
end

"""
    _add_subdivision_daughter!(tree, parent_idx, d_order, parent_dir, is_first, parent_order, params, rng)

Add a single subdivision daughter segment and recursively subdivide it.
"""
function _add_subdivision_daughter!(
    tree::VascularTree,
    parent_idx::Int,
    d_order::Int,
    parent_dir::NTuple{3,Float64},
    is_first::Bool,
    parent_order::Int,
    params::MorphometricParams,
    rng::AbstractRNG,
)
    seg = tree.segments

    # Sample daughter diameter and length from per-order distributions
    d_diam_um = max(
        params.diameter_mean[d_order + 1] + randn(rng) * params.diameter_sd[d_order + 1],
        params.vessel_cutoff_um,
    )
    d_length_um = max(
        params.length_mean[d_order + 1] + randn(rng) * params.length_sd[d_order + 1],
        params.vessel_cutoff_um,
    )

    d_radius = d_diam_um / 2.0 / 1000.0  # um → mm
    d_length = d_length_um / 1000.0        # um → mm

    # Compute direction
    d_dir = _random_daughter_direction(parent_dir, is_first, d_order, parent_order, rng)

    # Compute endpoints
    cp_x = seg.distal_x[parent_idx]
    cp_y = seg.distal_y[parent_idx]
    cp_z = seg.distal_z[parent_idx]

    proximal = (cp_x, cp_y, cp_z)
    distal = (
        cp_x + d_dir[1] * d_length,
        cp_y + d_dir[2] * d_length,
        cp_z + d_dir[3] * d_length,
    )

    new_id = add_segment!(tree, proximal, distal, d_radius, Int32(parent_idx))

    # Recursively subdivide if order > 0
    if d_order > 0
        _subdivide_recursive!(tree, Int(new_id), d_order, params, rng)
    end
end
