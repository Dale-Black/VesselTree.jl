# Statistical subdivision engine — Kassab connectivity matrix-based
# Recursively subdivides CCO skeleton terminals into full microvasculature
#
# Uses bifurcation chains (not cascade stubs) to produce clean binary trees.
# Each element's daughters are peeled off one at a time in a chain of bifurcations,
# with continuation segments forming the parent element. This naturally produces
# the S/E ratio (segments per element) observed in Kassab 1993.

"""
    _shell_tangent_project(domain::EllipsoidShellDomain, point, dir) -> NTuple{3,Float64}

Project a direction vector onto the tangent plane of the ellipsoid shell at the
given point. Removes the radial (surface-normal) component so that segments
follow the shell surface rather than escaping through it.
"""
function _shell_tangent_project(domain::EllipsoidShellDomain, point, dir)
    a, b, c = domain.semi_axes
    cx, cy, cz = domain.center
    # Outward normal (gradient of x²/a² + y²/b² + z²/c²)
    nx = (point[1] - cx) / (a * a)
    ny = (point[2] - cy) / (b * b)
    nz = (point[3] - cz) / (c * c)
    nlen = sqrt(nx^2 + ny^2 + nz^2)
    nlen < 1e-15 && return dir
    nx /= nlen; ny /= nlen; nz /= nlen

    # Remove radial component: d_tangent = d - (d·n)n
    dot = dir[1] * nx + dir[2] * ny + dir[3] * nz
    tx = dir[1] - dot * nx
    ty = dir[2] - dot * ny
    tz = dir[3] - dot * nz
    tlen = sqrt(tx^2 + ty^2 + tz^2)
    tlen < 1e-15 && return dir  # pure radial — keep original
    return (tx / tlen, ty / tlen, tz / tlen)
end

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
Continuation segments get small perturbation. Branches get offsets based on
Poiseuille-optimal angles (Murray's law) with some noise.

Angle ranges (physiological):
- Continuation: 0–10° perturbation
- Sprouting (rho < 0.83): 40–75° (side branch)
- Branching (rho ≥ 0.83): 20–40° (Y-fork)
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
        # Continuation: small perturbation (< 10 degrees)
        max_angle = 0.17  # ~10 degrees
    else
        # Branch: angle depends on order ratio (Barabasi-like)
        rho = parent_order > 0 ? Float64(d_order) / Float64(parent_order) : 0.0
        if rho < 0.83
            # Sprouting: moderate-to-large angle (40°–75°)
            # Murray's law: small daughter deflects ~50-70° for typical asymmetry
            base_angle = 0.87  # ~50 degrees
            noise = 0.44 * rand(rng) - 0.17  # range: -10° to +15°
            max_angle = clamp(base_angle + noise, 0.70, 1.31)  # 40°–75°
        else
            # Branching: moderate angle (20°–40°)
            # Murray's law: ~25° for symmetric, up to ~35° for mild asymmetry
            base_angle = 0.44  # ~25 degrees
            noise = 0.26 * rand(rng) - 0.09  # range: -5° to +10°
            max_angle = clamp(base_angle + noise, 0.35, 0.70)  # 20°–40°
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

Sample segment length (mm) from per-order lognormal distribution.
Kassab length data is right-skewed; lognormal matches the empirical shape.
Parameters are converted from the stored Normal (mean, sd) to lognormal.
"""
function _sample_segment_length(order::Int, params::MorphometricParams, rng::AbstractRNG)
    idx = order + 1
    if idx < 1 || idx > length(params.length_mean)
        return params.vessel_cutoff_um / 1000.0
    end
    mean_um = params.length_mean[idx]
    sd_um = params.length_sd[idx]
    # Convert Normal(mean, sd) → LogNormal(μ_ln, σ_ln)
    cv2 = (sd_um / mean_um)^2
    σ_ln = sqrt(log1p(cv2))
    μ_ln = log(mean_um) - 0.5 * σ_ln^2
    min_um = max(params.vessel_cutoff_um, mean_um * 0.1)
    for _ in 1:10
        l = exp(μ_ln + σ_ln * randn(rng))
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

function _sample_cm_daughter_order(parent_order::Int, params::MorphometricParams, rng::AbstractRNG)
    CM = params.connectivity_matrix
    n_orders = params.n_orders
    daughter_orders = Int[]
    weights = Float64[]
    for m in 0:(parent_order - 1)
        if m + 1 <= n_orders && parent_order + 1 <= n_orders
            lambda = CM[m + 1, parent_order + 1]
            lambda > 0.0 || continue
            push!(daughter_orders, m)
            push!(weights, lambda)
        end
    end

    isempty(daughter_orders) && return max(parent_order - 1, 0)

    total = sum(weights)
    draw = rand(rng) * total
    accum = 0.0
    for (idx, w) in enumerate(weights)
        accum += w
        if draw <= accum
            return daughter_orders[idx]
        end
    end
    return daughter_orders[end]
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
    domain::Union{AbstractDomain, Nothing}=nothing,
    enforce_full_cutoff::Bool=false,
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
        _subdivide_recursive!(
            tree,
            term_indices[k],
            term_orders[k],
            params,
            rng,
            domain;
            enforce_full_cutoff=enforce_full_cutoff,
        )
    end

    return tree
end

"""
    _subdivide_recursive!(tree, parent_idx, parent_order, params, rng)

Recursively subdivide a single segment into daughters per Kassab CM.
Uses a bifurcation chain with a terminal bifurcation:
- Internal junctions: daughter peels off + continuation maintains parent order
- Terminal junction: both children are new daughter elements (no continuation)

For N daughters: N-2 internal bifurcations + 1 terminal → S/E = N-1 (for N>=2).
This matches Kassab's S/E ratios where CM column sum ≈ S/E + 1.
"""
function _subdivide_recursive!(
    tree::VascularTree,
    parent_idx::Int,
    parent_order::Int,
    params::MorphometricParams,
    rng::AbstractRNG,
    domain::Union{AbstractDomain, Nothing}=nothing,
    ;
    enforce_full_cutoff::Bool=false,
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

    if enforce_full_cutoff && parent_order > 0
        while length(daughters) < 2
            push!(daughters, _sample_cm_daughter_order(parent_order, params, rng))
        end
    elseif isempty(daughters)
        return
    end
    # Sort daughters by order ascending — peel off smallest first (most asymmetric)
    sort!(daughters)

    # Build bifurcation chain with terminal bifurcation.
    # In Kassab's framework, an element of K segments has:
    #   - K-1 internal junctions (each peels off 1 daughter + 1 continuation)
    #   - 1 terminal junction (both children are new elements, no continuation)
    # This gives S/E = K, and total daughters = K+1 ≈ CM column sum.
    #
    # For N daughters drawn from CM:
    #   N >= 2: first N-2 internal bifurcations + 1 terminal (2 daughters)
    #           S/E = (N-2) + 1 = N-1
    #   N = 1:  single internal bifurcation (daughter + continuation), S/E = 2
    #   N = 0:  terminal leaf (handled above)
    current_parent = parent_idx
    n_daughters = length(daughters)

    # Number of internal bifurcations (daughter + continuation)
    n_internal = n_daughters >= 2 ? n_daughters - 2 : n_daughters

    for i in 1:n_internal
        d_order = daughters[i]
        r_parent = seg.radius[current_parent]

        # CM-implied asymmetry: D_elem(daughter)/D_elem(parent) with noise
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

        if domain isa EllipsoidShellDomain
            d_dir = _shell_tangent_project(domain, bp, d_dir)
        end

        d_distal = (bp_x + d_dir[1] * d_length, bp_y + d_dir[2] * d_length, bp_z + d_dir[3] * d_length)

        if domain !== nothing && !in_domain(domain, d_distal)
            d_distal = project_to_domain(domain, d_distal)
        end

        d_id = add_segment!(tree, bp, d_distal, d_radius, Int32(current_parent))
        tree.topology.strahler_order[d_id] = Int32(d_order)

        # Create continuation segment (same order as parent element)
        c_length = _sample_segment_length(parent_order, params, rng)
        c_dir = _random_daughter_direction(current_dir, true, parent_order, parent_order, rng)

        if domain isa EllipsoidShellDomain
            c_dir = _shell_tangent_project(domain, bp, c_dir)
        end

        c_distal = (bp_x + c_dir[1] * c_length, bp_y + c_dir[2] * c_length, bp_z + c_dir[3] * c_length)

        if domain !== nothing && !in_domain(domain, c_distal)
            c_distal = project_to_domain(domain, c_distal)
        end

        c_id = add_segment!(tree, bp, c_distal, c_radius, Int32(current_parent))
        tree.topology.strahler_order[c_id] = Int32(parent_order)

        # Update direction from actual (possibly projected) continuation segment
        actual_cx = c_distal[1] - bp_x
        actual_cy = c_distal[2] - bp_y
        actual_cz = c_distal[3] - bp_z
        actual_len = sqrt(actual_cx^2 + actual_cy^2 + actual_cz^2)
        if actual_len > 1e-15
            c_dir = (actual_cx / actual_len, actual_cy / actual_len, actual_cz / actual_len)
        end

        # Recursively subdivide the daughter
        if d_order > 0
            _subdivide_recursive!(tree, Int(d_id), d_order, params, rng, domain; enforce_full_cutoff=enforce_full_cutoff)
        end

        # Move to continuation for next bifurcation in the chain
        current_parent = Int(c_id)
        current_dir = c_dir
    end

    # Terminal bifurcation: last 2 daughters share a junction (no continuation).
    # The element ends here — both children start new elements.
    if n_daughters >= 2
        d1_order = daughters[n_daughters - 1]
        d2_order = daughters[n_daughters]
        r_parent = seg.radius[current_parent]

        bp_x = seg.distal_x[current_parent]
        bp_y = seg.distal_y[current_parent]
        bp_z = seg.distal_z[current_parent]
        bp = (bp_x, bp_y, bp_z)

        # Split radius between 2 daughters using their relative element diameters
        if d1_order >= d2_order
            asym = _cm_implied_asymmetry(d2_order, max(d1_order, 1), params, rng)
            r_large, r_small = compute_daughter_radii(r_parent, asym, gamma)
            d1_radius, d2_radius = r_large, r_small
        else
            asym = _cm_implied_asymmetry(d1_order, max(d2_order, 1), params, rng)
            r_large, r_small = compute_daughter_radii(r_parent, asym, gamma)
            d1_radius, d2_radius = r_small, r_large
        end

        # Daughter 1
        d1_length = _sample_segment_length(d1_order, params, rng)
        d1_dir = _random_daughter_direction(current_dir, false, d1_order, parent_order, rng)
        if domain isa EllipsoidShellDomain
            d1_dir = _shell_tangent_project(domain, bp, d1_dir)
        end
        d1_distal = (bp_x + d1_dir[1] * d1_length, bp_y + d1_dir[2] * d1_length, bp_z + d1_dir[3] * d1_length)
        if domain !== nothing && !in_domain(domain, d1_distal)
            d1_distal = project_to_domain(domain, d1_distal)
        end
        d1_id = add_segment!(tree, bp, d1_distal, d1_radius, Int32(current_parent))
        tree.topology.strahler_order[d1_id] = Int32(d1_order)

        # Daughter 2
        d2_length = _sample_segment_length(d2_order, params, rng)
        d2_dir = _random_daughter_direction(current_dir, false, d2_order, parent_order, rng)
        if domain isa EllipsoidShellDomain
            d2_dir = _shell_tangent_project(domain, bp, d2_dir)
        end
        d2_distal = (bp_x + d2_dir[1] * d2_length, bp_y + d2_dir[2] * d2_length, bp_z + d2_dir[3] * d2_length)
        if domain !== nothing && !in_domain(domain, d2_distal)
            d2_distal = project_to_domain(domain, d2_distal)
        end
        d2_id = add_segment!(tree, bp, d2_distal, d2_radius, Int32(current_parent))
        tree.topology.strahler_order[d2_id] = Int32(d2_order)

        # Recursively subdivide both daughters
        d1_order > 0 && _subdivide_recursive!(tree, Int(d1_id), d1_order, params, rng, domain; enforce_full_cutoff=enforce_full_cutoff)
        d2_order > 0 && _subdivide_recursive!(tree, Int(d2_id), d2_order, params, rng, domain; enforce_full_cutoff=enforce_full_cutoff)
    end
    # For N=1: the internal loop above creates one daughter plus a continuation
    # segment that is still part of the same parent-order element. Continue
    # subdividing that continuation so the chain does not terminate prematurely
    # at a high-order stub.
    if n_daughters == 1 && parent_order > 0
        _subdivide_recursive!(tree, current_parent, parent_order, params, rng, domain; enforce_full_cutoff=enforce_full_cutoff)
    end
    # For N=0: handled by the isempty check above
end
