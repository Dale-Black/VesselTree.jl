# Statistical subdivision engine — Kassab connectivity matrix-based
# Recursively subdivides CCO skeleton terminals into full microvasculature
#
# Uses bifurcation chains (not cascade stubs) to produce clean binary trees.
# Each element's daughters are peeled off one at a time in a chain of bifurcations,
# with continuation segments forming the parent element. This naturally produces
# the S/E ratio (segments per element) observed in Kassab 1993.

"""
    estimate_total_segments(parent_order, params) -> Int

Estimate total segments produced by recursively subdividing a parent of given
Strahler order using the connectivity matrix. Accounts for continuation segments
in bifurcation chains. Uses expected values (no randomness).
"""
function estimate_total_segments(parent_order::Int, params::MorphometricParams)
    parent_order <= 0 && return 1
    CM = params.connectivity_matrix
    n_orders = params.n_orders
    total = 1.0  # the parent itself
    for m in 0:(parent_order - 1)
        if m + 1 <= n_orders && parent_order + 1 <= n_orders
            avg = CM[m + 1, parent_order + 1]
            if avg > 0
                # Each daughter creates 1 daughter segment + 1 continuation segment
                # The daughter recursively expands; the continuation stays as-is
                total += avg * (1.0 + estimate_total_segments(m, params))
            end
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
    _random_daughter_direction(parent_dir, is_continuation, d_order, parent_order, rng)

Compute a direction for a daughter segment based on parent direction.
Continuation segments get small perturbation. Branches get larger offsets
based on Barabasi sprouting/branching regimes.
"""
function _random_daughter_direction(
    parent_dir::NTuple{3,Float64},
    is_continuation::Bool,
    d_order::Int,
    parent_order::Int,
    rng::AbstractRNG,
)
    dx, dy, dz = parent_dir

    if is_continuation
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
    _cm_implied_asymmetry(d_order, parent_order, params, rng) -> Float64

Compute asymmetry ratio (r_small/r_large) from Kassab element diameters.
Uses D_elem(daughter)/D_elem(parent) as mean, with noise derived from
the empirical coefficient of variation of both diameter distributions.

This replaces the generic Beta(2.5, 0.8) with order-specific ratios
grounded in Kassab 1993 Table 1 element data.
"""
function _cm_implied_asymmetry(
    d_order::Int, parent_order::Int,
    params::MorphometricParams, rng::AbstractRNG,
)
    d_idx = d_order + 1
    p_idx = parent_order + 1
    # Bounds check — fall back to generic 0.5 if out of range
    if d_idx < 1 || d_idx > length(params.diameter_mean_elem) ||
       p_idx < 1 || p_idx > length(params.diameter_mean_elem)
        return 0.5
    end

    d_mean = params.diameter_mean_elem[d_idx]
    p_mean = params.diameter_mean_elem[p_idx]
    d_sd = params.diameter_sd_elem[d_idx]
    p_sd = params.diameter_sd_elem[p_idx]

    # Mean asymmetry = daughter_diameter / parent_diameter
    mean_asym = d_mean / p_mean

    # Propagate uncertainty: CV_ratio = sqrt(CV_d^2 + CV_p^2)
    cv_d = d_sd / max(d_mean, 1e-10)
    cv_p = p_sd / max(p_mean, 1e-10)
    cv_ratio = sqrt(cv_d^2 + cv_p^2)

    # Sample with truncated Normal noise
    for _ in 1:10
        asym = mean_asym * (1.0 + cv_ratio * randn(rng))
        if asym > 0.02 && asym < 0.99
            return asym
        end
    end
    return clamp(mean_asym, 0.02, 0.99)
end

"""
    _sample_segment_length(order, params, rng) -> Float64

Sample segment length (mm) from per-order distribution. Uses truncated Normal
to avoid negative or unreasonably small lengths.
"""
function _sample_segment_length(order::Int, params::MorphometricParams, rng::AbstractRNG)
    idx = order + 1
    if idx < 1 || idx > length(params.length_mean)
        return params.vessel_cutoff_um / 1000.0
    end
    mean_um = params.length_mean[idx]
    sd_um = params.length_sd[idx]
    min_um = max(params.vessel_cutoff_um, mean_um * 0.1)  # at least 10% of mean
    # Truncated Normal: sample until we get a valid value (simple rejection)
    for _ in 1:10
        l = mean_um + randn(rng) * sd_um
        l >= min_um && return l / 1000.0  # um → mm
    end
    return mean_um / 1000.0  # fallback to mean
end

"""
    _sample_element_diameter(order, params, rng) -> Float64

Sample element diameter (um) from Kassab element-level distribution.
Returns a positive diameter, truncated to avoid unreasonably small values.
"""
function _sample_element_diameter(order::Int, params::MorphometricParams, rng::AbstractRNG)
    idx = order + 1
    if idx < 1 || idx > length(params.diameter_mean_elem)
        return params.vessel_cutoff_um
    end
    mean_um = params.diameter_mean_elem[idx]
    sd_um = params.diameter_sd_elem[idx]
    min_um = max(params.vessel_cutoff_um * 0.5, mean_um * 0.1)
    for _ in 1:10
        d = mean_um + randn(rng) * sd_um
        d >= min_um && return d
    end
    return mean_um
end

"""
    subdivide_terminals!(tree, params; rng, max_order)

Recursively subdivide each terminal of order > 0 down to order 0 (capillaries)
using Kassab's connectivity matrix. Uses bifurcation chains: each daughter
peeled off creates a proper bifurcation with a continuation segment, producing
trees with 100% bifurcations (no cascade stubs).

Asymmetry is baked in: at each bifurcation, CM-implied asymmetry from
Kassab element diameters determines daughter radii via Murray's law.
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
Uses a bifurcation chain: each daughter peels off via a proper bifurcation
with a continuation segment that maintains the parent's order.

For N daughters, creates N bifurcations (each with daughter + continuation).
The continuation segments form the parent element (S/E = N+1).
Only the daughter segments are recursively subdivided.
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
    current_dir = (dir_x, dir_y, dir_z)

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
    # Sort daughters by order ascending — peel off smallest first (most asymmetric)
    sort!(daughters)

    # Build bifurcation chain: each daughter gets its own bifurcation point
    current_parent = parent_idx

    for i in eachindex(daughters)
        d_order = daughters[i]
        r_parent = seg.radius[current_parent]

        # CM-implied asymmetry: D_elem(daughter)/D_elem(parent) with noise
        # Grounded in Kassab 1993 element diameters, not generic Beta distribution.
        # No floor clamp: Murray's law holds exactly at every junction.
        asymmetry = _cm_implied_asymmetry(d_order, parent_order, params, rng)
        r_large, r_small = compute_daughter_radii(r_parent, asymmetry, gamma)

        # Daughter (branch) gets the small radius; continuation gets the large radius
        d_radius = r_small
        c_radius = r_large

        # Bifurcation point: distal end of current parent
        bp_x = seg.distal_x[current_parent]
        bp_y = seg.distal_y[current_parent]
        bp_z = seg.distal_z[current_parent]
        bp = (bp_x, bp_y, bp_z)

        # Create daughter segment (branch)
        d_length = _sample_segment_length(d_order, params, rng)
        d_dir = _random_daughter_direction(current_dir, false, d_order, parent_order, rng)
        d_distal = (bp_x + d_dir[1] * d_length, bp_y + d_dir[2] * d_length, bp_z + d_dir[3] * d_length)
        d_id = add_segment!(tree, bp, d_distal, d_radius, Int32(current_parent))
        tree.topology.strahler_order[d_id] = Int32(d_order)

        # Create continuation segment (same order as parent element)
        c_length = _sample_segment_length(parent_order, params, rng)
        c_dir = _random_daughter_direction(current_dir, true, parent_order, parent_order, rng)
        c_distal = (bp_x + c_dir[1] * c_length, bp_y + c_dir[2] * c_length, bp_z + c_dir[3] * c_length)
        c_id = add_segment!(tree, bp, c_distal, c_radius, Int32(current_parent))
        tree.topology.strahler_order[c_id] = Int32(parent_order)

        # Recursively subdivide the daughter
        if d_order > 0
            _subdivide_recursive!(tree, Int(d_id), d_order, params, rng)
        end

        # Move to continuation for next bifurcation in the chain
        current_parent = Int(c_id)
        current_dir = c_dir
    end

    # current_parent is the terminal end of the element — NOT further subdivided
    # (all daughters for this element have been accounted for by the CM entry)
end
