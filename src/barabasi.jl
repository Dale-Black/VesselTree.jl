# Barabasi 2026 surface optimization: junction geometry

"""
    compute_chi(parent_radius, node_distance) -> Float64

Compute chi parameter: chi = 2*pi*r / d.
chi > 0.83 indicates trifurcation is stable.
"""
function compute_chi(parent_radius::Float64, node_distance::Float64)
    node_distance <= 0.0 && return Inf
    return 2π * parent_radius / node_distance
end

"""
    compute_rho(r_small, r_large) -> Float64

Compute rho parameter: rho = r_small / r_large.
Determines sprouting (rho < 0.83) vs branching (rho >= 0.83).
"""
function compute_rho(r_small::Float64, r_large::Float64)
    r_large <= 0.0 && return 0.0
    return r_small / r_large
end

"""
    classify_junction(rho, params) -> Symbol

Classify junction type based on rho parameter.
Returns :sprouting or :branching.
"""
function classify_junction(rho::Float64, params::MorphometricParams)
    return rho < params.sprouting_rho_th ? :sprouting : :branching
end

"""
    compute_junction_angles(rho, params) -> (angle_large, angle_small)

Compute daughter branch angles based on rho parameter.
- Sprouting (rho < rho_th): large continues straight (0), small at ~90 deg
- Branching (rho >= rho_th): both deflect, angles depend on rho

Returns (angle_large, angle_small) in radians.
"""
function compute_junction_angles(rho::Float64, params::MorphometricParams)
    rho_th = params.sprouting_rho_th

    if rho < rho_th
        # Sprouting regime: large daughter continues straight, small branches perpendicular
        angle_large = 0.0
        angle_small = π / 2.0
        return (angle_large, angle_small)
    end

    # Branching regime: both daughters deflect
    # Murray's law optimal angles using cos formula:
    # cos(alpha_large) = (r_parent^4 + r_large^4 - r_small^4) / (2 * r_parent^2 * r_large^2)
    # Simplified using rho = r_small/r_large:
    # For the branching regime, use a linear interpolation from threshold to symmetric
    t = (rho - rho_th) / (1.0 - rho_th + 1e-10)  # 0 at threshold, 1 at symmetric

    # At rho_th: minimal deflection; at rho=1: equal deflection
    # Murray angle for symmetric bifurcation with gamma=7/3:
    # cos(alpha) = (1 + 1 - 1) / (2*1*1) = 0.5 → alpha ≈ 60 deg → but with gamma=7/3 it differs
    # Use energy-optimal formula: angle ≈ acos((1 + r_ratio^4 - (rho*r_ratio)^4)/(2*r_ratio^2))
    # Simplified: symmetric angle ≈ 37 degrees (0.65 rad) for gamma=7/3

    symmetric_angle = 0.65  # ~37 deg for gamma=7/3, both daughters
    min_angle = 0.05        # near-zero at threshold

    angle_total = min_angle + t * (2 * symmetric_angle - 2 * min_angle)

    # Distribute: small daughter gets proportionally more deflection
    # Weight by 1/(r^2) — smaller vessel deflects more
    w_large = rho        # proportional to r_small (which means large vessel gets less)
    w_small = 1.0

    total_w = w_large + w_small
    angle_large = angle_total * w_large / total_w / 2.0
    angle_small = angle_total * w_small / total_w / 2.0

    # Ensure small >= large
    if angle_small < angle_large
        angle_small, angle_large = angle_large, angle_small
    end

    return (angle_large, angle_small)
end

"""
    apply_junction_geometry!(tree::VascularTree, bifurc_idx::Int32, params::MorphometricParams)

Apply Barabasi junction geometry to a bifurcation.
Rotates daughter segment endpoints to match sprouting/branching angles
while preserving segment lengths.
"""
function apply_junction_geometry!(tree::VascularTree, bifurc_idx::Int32, params::MorphometricParams)
    topo = tree.topology
    seg = tree.segments
    idx = Int(bifurc_idx)

    topo.junction_type[idx] != :bifurcation && return nothing

    c1 = Int(topo.child1_id[idx])
    c2 = Int(topo.child2_id[idx])
    (c1 <= 0 || c2 <= 0) && return nothing

    # Determine which is large (continuation) and which is small (new branch)
    r1 = seg.radius[c1]
    r2 = seg.radius[c2]
    if r1 >= r2
        large_id, small_id = c1, c2
    else
        large_id, small_id = c2, c1
    end

    r_large = seg.radius[large_id]
    r_small = seg.radius[small_id]
    rho = compute_rho(r_small, r_large)

    angle_large, angle_small = compute_junction_angles(rho, params)

    # Parent direction vector
    px = seg.proximal_x[idx]
    py = seg.proximal_y[idx]
    pz = seg.proximal_z[idx]
    dx = seg.distal_x[idx]
    dy = seg.distal_y[idx]
    dz = seg.distal_z[idx]

    parent_dx = dx - px
    parent_dy = dy - py
    parent_dz = dz - pz
    parent_len = sqrt(parent_dx^2 + parent_dy^2 + parent_dz^2)
    parent_len <= 0.0 && return nothing

    # Normalize parent direction
    parent_dx /= parent_len
    parent_dy /= parent_len
    parent_dz /= parent_len

    # Find a perpendicular direction (for rotation plane)
    perp_x, perp_y, perp_z = _find_perpendicular(parent_dx, parent_dy, parent_dz)

    # Bifurcation point (distal of parent = proximal of children)
    bx, by, bz = dx, dy, dz

    # Rotate large daughter
    len_large = seg.seg_length[large_id]
    new_dx_l = parent_dx * cos(angle_large) + perp_x * sin(angle_large)
    new_dy_l = parent_dy * cos(angle_large) + perp_y * sin(angle_large)
    new_dz_l = parent_dz * cos(angle_large) + perp_z * sin(angle_large)

    seg.distal_x[large_id] = bx + new_dx_l * len_large
    seg.distal_y[large_id] = by + new_dy_l * len_large
    seg.distal_z[large_id] = bz + new_dz_l * len_large

    # Rotate small daughter (opposite side of perpendicular)
    len_small = seg.seg_length[small_id]
    new_dx_s = parent_dx * cos(angle_small) - perp_x * sin(angle_small)
    new_dy_s = parent_dy * cos(angle_small) - perp_y * sin(angle_small)
    new_dz_s = parent_dz * cos(angle_small) - perp_z * sin(angle_small)

    seg.distal_x[small_id] = bx + new_dx_s * len_small
    seg.distal_y[small_id] = by + new_dy_s * len_small
    seg.distal_z[small_id] = bz + new_dz_s * len_small

    # Recompute segment lengths (should be preserved, but recalculate for safety)
    _update_seg_length!(seg, large_id)
    _update_seg_length!(seg, small_id)

    return nothing
end

"""
    _find_perpendicular(dx, dy, dz) -> (px, py, pz)

Find a unit vector perpendicular to (dx, dy, dz).
"""
function _find_perpendicular(dx::Float64, dy::Float64, dz::Float64)
    # Choose the axis with smallest component for cross product
    if abs(dx) <= abs(dy) && abs(dx) <= abs(dz)
        ax, ay, az = 1.0, 0.0, 0.0
    elseif abs(dy) <= abs(dz)
        ax, ay, az = 0.0, 1.0, 0.0
    else
        ax, ay, az = 0.0, 0.0, 1.0
    end

    # Cross product: (d) × (a)
    px = dy * az - dz * ay
    py = dz * ax - dx * az
    pz = dx * ay - dy * ax

    # Normalize
    len = sqrt(px^2 + py^2 + pz^2)
    return (px / len, py / len, pz / len)
end

"""
    _update_seg_length!(seg::SegmentData, idx::Int)

Recompute seg_length for segment idx from its endpoints.
"""
function _update_seg_length!(seg::SegmentData, idx::Int)
    ex = seg.distal_x[idx] - seg.proximal_x[idx]
    ey = seg.distal_y[idx] - seg.proximal_y[idx]
    ez = seg.distal_z[idx] - seg.proximal_z[idx]
    seg.seg_length[idx] = sqrt(ex^2 + ey^2 + ez^2)
end

# --- Trifurcation detection and handling ---

"""
    compute_trifurcation_angles(r1, r2, r3, params) -> (a1, a2, a3)

Compute deflection angles for three daughters at a trifurcation.
Daughters are sorted by radius (largest first). Larger daughters get smaller deflection.
Returns angles in radians ordered by (r1, r2, r3) input order.
"""
function compute_trifurcation_angles(r1::Float64, r2::Float64, r3::Float64, params::MorphometricParams)
    radii = [r1, r2, r3]
    perm = sortperm(radii, rev=true)  # largest first

    # Total deflection budget: ~120 deg spread, distributed by inverse radius
    # For symmetric case (all equal), each gets ~40 deg (~0.7 rad)
    symmetric_angle = 0.7  # ~40 deg per daughter for symmetric trifurcation

    # Weight by inverse radius: smaller vessels deflect more
    weights = [1.0 / (radii[perm[k]] + 1e-20) for k in 1:3]
    total_w = sum(weights)

    # Compute angles — distribute the total angular budget
    total_angle = 3.0 * symmetric_angle  # total budget
    angles_sorted = [(total_angle * weights[k] / total_w) for k in 1:3]

    # Un-permute to match input order
    angles = zeros(3)
    for k in 1:3
        angles[perm[k]] = angles_sorted[k]
    end

    return (angles[1], angles[2], angles[3])
end

"""
    check_trifurcation_merge(tree, new_point, params) -> Union{Nothing, Int32}

Check if a new bifurcation point should merge with a nearby existing bifurcation
to form a trifurcation. Returns the bifurcation segment ID to merge with, or nothing.

Criteria: the new point must be close to an existing bifurcation's distal point,
and chi = 2*pi*r / d must exceed the trifurcation threshold (0.83).
"""
function check_trifurcation_merge(tree::VascularTree, new_point::NTuple{3,Float64}, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n

    best_idx = nothing
    best_chi = -Inf

    for i in 1:n
        # Only consider existing bifurcations (not terminals, not already trifurcations)
        topo.junction_type[i] != :bifurcation && continue

        # Distance from new point to bifurcation point (distal of parent segment)
        ddx = new_point[1] - seg.distal_x[i]
        ddy = new_point[2] - seg.distal_y[i]
        ddz = new_point[3] - seg.distal_z[i]
        dist = sqrt(ddx^2 + ddy^2 + ddz^2)

        # Compute chi
        chi = compute_chi(seg.radius[i], dist)

        if chi > params.trifurcation_chi_th && chi > best_chi
            best_chi = chi
            best_idx = Int32(i)
        end
    end

    return best_idx
end

"""
    merge_to_trifurcation!(tree, bifurc_idx, terminal_point, terminal_radius, params) -> Int32

Convert an existing bifurcation to a trifurcation by adding a third daughter segment.
Returns the new segment ID.
"""
function merge_to_trifurcation!(
    tree::VascularTree,
    bifurc_idx::Int32,
    terminal_point::NTuple{3,Float64},
    terminal_radius::Float64,
    params::MorphometricParams,
)
    seg = tree.segments
    idx = Int(bifurc_idx)

    # Proximal point of new daughter = distal point of parent
    prox = (seg.distal_x[idx], seg.distal_y[idx], seg.distal_z[idx])

    # add_segment! handles the topology upgrade from bifurcation to trifurcation
    new_id = add_segment!(tree, prox, terminal_point, terminal_radius, bifurc_idx)

    return Int32(new_id)
end

"""
    apply_trifurcation_geometry!(tree, trifurc_idx, params)

Apply angle geometry to a trifurcation node. Rotates the three daughter endpoints
to computed angles while preserving segment lengths.
"""
function apply_trifurcation_geometry!(tree::VascularTree, trifurc_idx::Int32, params::MorphometricParams)
    topo = tree.topology
    seg = tree.segments
    idx = Int(trifurc_idx)

    topo.junction_type[idx] != :trifurcation && return nothing

    c1 = Int(topo.child1_id[idx])
    c2 = Int(topo.child2_id[idx])
    c3 = Int(topo.child3_id[idx])
    (c1 <= 0 || c2 <= 0 || c3 <= 0) && return nothing

    r1 = seg.radius[c1]
    r2 = seg.radius[c2]
    r3 = seg.radius[c3]

    a1, a2, a3 = compute_trifurcation_angles(r1, r2, r3, params)

    # Parent direction
    px = seg.proximal_x[idx]
    py = seg.proximal_y[idx]
    pz = seg.proximal_z[idx]
    dx = seg.distal_x[idx]
    dy = seg.distal_y[idx]
    dz = seg.distal_z[idx]

    parent_dx = dx - px
    parent_dy = dy - py
    parent_dz = dz - pz
    parent_len = sqrt(parent_dx^2 + parent_dy^2 + parent_dz^2)
    parent_len <= 0.0 && return nothing

    parent_dx /= parent_len
    parent_dy /= parent_len
    parent_dz /= parent_len

    # Two perpendicular directions for 3D angle placement
    perp1_x, perp1_y, perp1_z = _find_perpendicular(parent_dx, parent_dy, parent_dz)
    # Second perpendicular: parent × perp1
    perp2_x = parent_dy * perp1_z - parent_dz * perp1_y
    perp2_y = parent_dz * perp1_x - parent_dx * perp1_z
    perp2_z = parent_dx * perp1_y - parent_dy * perp1_x

    bx, by, bz = dx, dy, dz

    # Place three daughters at azimuthal angles 0, 2π/3, 4π/3 in the perpendicular plane
    children = [(c1, a1, 0.0), (c2, a2, 2π / 3), (c3, a3, 4π / 3)]

    for (cid, polar_angle, azimuth) in children
        len_c = seg.seg_length[cid]
        sp = sin(polar_angle)
        cp = cos(polar_angle)
        sa = sin(azimuth)
        ca = cos(azimuth)

        # Direction: rotate parent_dir by polar_angle, with azimuthal rotation in perp plane
        dir_x = parent_dx * cp + (perp1_x * ca + perp2_x * sa) * sp
        dir_y = parent_dy * cp + (perp1_y * ca + perp2_y * sa) * sp
        dir_z = parent_dz * cp + (perp1_z * ca + perp2_z * sa) * sp

        seg.distal_x[cid] = bx + dir_x * len_c
        seg.distal_y[cid] = by + dir_y * len_c
        seg.distal_z[cid] = bz + dir_z * len_c

        _update_seg_length!(seg, cid)
    end

    return nothing
end
