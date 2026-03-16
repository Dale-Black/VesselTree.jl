# Multi-tree forest: territory partitioning and collision avoidance

"""
    TreeConfig

Configuration for a single tree in a multi-tree forest.
"""
struct TreeConfig
    name::String
    root_position::NTuple{3, Float64}
    root_radius::Float64
    root_direction::NTuple{3, Float64}
    target_terminals::Int
    territory_fraction::Float64
end

"""
    TerritoryMap

3D voxel grid assigning each region to a tree by name.
"""
mutable struct TerritoryMap
    cell_size::Float64
    origin::NTuple{3, Float64}
    dims::NTuple{3, Int}
    cells::Vector{Int}           # flat array: cell -> tree index (1-based)
    tree_names::Vector{String}   # index -> name
end

function _territory_targets(configs::Vector{TreeConfig})
    raw = [max(cfg.territory_fraction, 0.0) for cfg in configs]
    total = sum(raw)
    total > 0.0 || return fill(1.0 / length(configs), length(configs))
    return raw ./ total
end

function _territory_owner(
    refs::Vector{NTuple{3, Float64}},
    biases::Vector{Float64},
    point::NTuple{3, Float64},
)
    px, py, pz = point
    best_idx = 1
    best_score = Inf
    for (ti, ref) in enumerate(refs)
        dx = px - ref[1]
        dy = py - ref[2]
        dz = pz - ref[3]
        d = sqrt(dx^2 + dy^2 + dz^2)
        score = d - biases[ti]
        if score < best_score
            best_score = score
            best_idx = ti
        end
    end
    return best_idx
end

function _inside_voxel_centers(
    domain::AbstractDomain,
    origin::NTuple{3, Float64},
    cell_size::Float64,
    dims::NTuple{3, Int},
)
    inside = Vector{Tuple{Int, NTuple{3, Float64}}}()
    nx, ny, nz = dims
    for iz in 1:nz
        for iy in 1:ny
            for ix in 1:nx
                cx = origin[1] + (ix - 0.5) * cell_size
                cy = origin[2] + (iy - 0.5) * cell_size
                cz = origin[3] + (iz - 0.5) * cell_size
                point = (cx, cy, cz)
                in_domain(domain, point) || continue
                flat_idx = ix + (iy - 1) * nx + (iz - 1) * nx * ny
                push!(inside, (flat_idx, point))
            end
        end
    end
    return inside
end

function _territory_reference_points(
    domain::AbstractDomain,
    configs::Vector{TreeConfig},
    anchor_length::Float64,
)
    refs = NTuple{3, Float64}[]
    for cfg in configs
        dx, dy, dz = cfg.root_direction
        dlen = sqrt(dx^2 + dy^2 + dz^2)
        if dlen <= 1e-12
            push!(refs, cfg.root_position)
            continue
        end
        ref = (
            cfg.root_position[1] + anchor_length * dx / dlen,
            cfg.root_position[2] + anchor_length * dy / dlen,
            cfg.root_position[3] + anchor_length * dz / dlen,
        )
        if !in_domain(domain, ref)
            ref = project_to_domain(domain, ref)
        end
        push!(refs, ref)
    end
    return refs
end

function _clear_local_intersection_neighborhood!(
    flags::Vector{Bool},
    topo,
    seg_idx::Int;
    hop_limit::Int=8,
)
    queue = Tuple{Int32, Int}[(Int32(seg_idx), 0)]
    seen = Set{Int32}()

    while !isempty(queue)
        current, depth = popfirst!(queue)
        current > 0 || continue
        current in seen && continue
        push!(seen, current)
        flags[current] = false

        depth >= hop_limit && continue

        parent = topo.parent_id[current]
        parent > 0 && push!(queue, (parent, depth + 1))

        for child in (topo.child1_id[current], topo.child2_id[current], topo.child3_id[current])
            child > 0 && push!(queue, (child, depth + 1))
        end
    end

    return flags
end

"""
    initialize_territories(domain, configs) -> TerritoryMap

Create a territory map by assigning each voxel to a root using a weighted
distance partition. The weights are iteratively adjusted so the realized
in-domain voxel fractions approximate each tree's configured
`territory_fraction`.
"""
function initialize_territories(domain::AbstractDomain, configs::Vector{TreeConfig})
    lo, hi = domain_bounds(domain)

    # Cell size: ~20 cells per dimension
    extent = max(hi[1] - lo[1], hi[2] - lo[2], hi[3] - lo[3])
    cell_size = extent / 20.0

    nx = max(1, ceil(Int, (hi[1] - lo[1]) / cell_size))
    ny = max(1, ceil(Int, (hi[2] - lo[2]) / cell_size))
    nz = max(1, ceil(Int, (hi[3] - lo[3]) / cell_size))

    tree_names = [cfg.name for cfg in configs]
    n_cells = nx * ny * nz
    cells = zeros(Int, n_cells)
    inside = _inside_voxel_centers(domain, lo, cell_size, (nx, ny, nz))
    anchor_length = 0.35 * extent
    refs = _territory_reference_points(domain, configs, anchor_length)

    biases = zeros(Float64, length(configs))
    targets = _territory_targets(configs)
    if !isempty(inside) && length(configs) > 1
        counts = zeros(Int, length(configs))
        gain = 0.35 * cell_size * max(nx, ny, nz)
        for _ in 1:48
            fill!(counts, 0)
            for (_, point) in inside
                owner = _territory_owner(refs, biases, point)
                counts[owner] += 1
            end

            total_inside = sum(counts)
            total_inside == 0 && break
            actual = counts ./ total_inside
            maximum(abs.(actual .- targets)) <= 0.02 && break

            for i in eachindex(biases)
                biases[i] += gain * (targets[i] - actual[i])
            end

            biases .-= sum(biases) / length(biases)
        end
    end

    # Assign each cell with the calibrated weighted partition.
    for iz in 1:nz
        for iy in 1:ny
            for ix in 1:nx
                cx = lo[1] + (ix - 0.5) * cell_size
                cy = lo[2] + (iy - 0.5) * cell_size
                cz = lo[3] + (iz - 0.5) * cell_size
                flat_idx = ix + (iy - 1) * nx + (iz - 1) * nx * ny
                cells[flat_idx] = _territory_owner(refs, biases, (cx, cy, cz))
            end
        end
    end

    return TerritoryMap(cell_size, lo, (nx, ny, nz), cells, tree_names)
end

"""
    query_territory(tmap, x, y, z) -> String

Return the tree name that owns the voxel containing (x, y, z).
"""
function query_territory(tmap::TerritoryMap, x::Float64, y::Float64, z::Float64)
    ix = clamp(floor(Int, (x - tmap.origin[1]) / tmap.cell_size) + 1, 1, tmap.dims[1])
    iy = clamp(floor(Int, (y - tmap.origin[2]) / tmap.cell_size) + 1, 1, tmap.dims[2])
    iz = clamp(floor(Int, (z - tmap.origin[3]) / tmap.cell_size) + 1, 1, tmap.dims[3])
    flat_idx = ix + (iy - 1) * tmap.dims[1] + (iz - 1) * tmap.dims[1] * tmap.dims[2]
    tree_idx = tmap.cells[flat_idx]
    return tmap.tree_names[tree_idx]
end

"""
    sample_in_territory(domain, tmap, tree_name, rng; max_attempts=100) -> Union{Nothing, NTuple{3,Float64}}

Sample a random point within the territory assigned to `tree_name`.
"""
function sample_in_territory(
    domain::AbstractDomain,
    tmap::TerritoryMap,
    tree_name::String,
    rng::AbstractRNG;
    max_attempts::Int=100,
)
    for _ in 1:max_attempts
        pt = sample_point(domain, rng)
        if query_territory(tmap, pt[1], pt[2], pt[3]) == tree_name
            return pt
        end
    end
    return nothing
end

"""
    check_inter_tree_collision(point, other_tree, min_distance) -> Bool

Check if a point is too close to any segment in another tree.
Uses AK distance kernel.
"""
function check_inter_tree_collision(
    point::NTuple{3, Float64},
    other_tree::VascularTree,
    min_distance::Float64,
)
    n = other_tree.segments.n
    n == 0 && return false

    distances_buf = zeros(n)
    compute_all_distances!(distances_buf, other_tree.segments, point[1], point[2], point[3], n)
    min_dist = AK.minimum(@view distances_buf[1:n])
    return min_dist < min_distance
end

"""
    coronary_tree_configs() -> Vector{TreeConfig}

Default configuration for coronary artery forest: LAD (40%), LCX (25%), RCA (35%).
Root positions approximate coronary ostia on a heart-sized domain.
"""
function coronary_tree_configs()
    return [
        TreeConfig("LAD", (-2.0, 1.0, 0.0), 1.5, (0.0, -1.0, -0.5), 1000, 0.40),
        TreeConfig("LCX", (-2.0, -1.0, 0.0), 1.2, (0.0, -1.0, 0.5), 600, 0.25),
        TreeConfig("RCA", (2.0, 0.0, 0.0), 1.3, (0.0, -1.0, 0.0), 800, 0.35),
    ]
end

coronary_tree_configs(::AbstractDomain) = coronary_tree_configs()

function _rotate_about_axis(v::NTuple{3, Float64}, axis::NTuple{3, Float64}, theta::Float64)
    ax, ay, az = axis
    anorm = sqrt(ax * ax + ay * ay + az * az)
    anorm > 1e-12 || return v
    ax /= anorm; ay /= anorm; az /= anorm

    vx, vy, vz = v
    c = cos(theta)
    s = sin(theta)
    dot_av = ax * vx + ay * vy + az * vz

    return (
        vx * c + (ay * vz - az * vy) * s + ax * dot_av * (1.0 - c),
        vy * c + (az * vx - ax * vz) * s + ay * dot_av * (1.0 - c),
        vz * c + (ax * vy - ay * vx) * s + az * dot_av * (1.0 - c),
    )
end

function _make_perp_axis(v::NTuple{3, Float64})
    vx, vy, vz = v
    vnorm = sqrt(vx * vx + vy * vy + vz * vz)
    vnorm > 1e-12 || return (0.0, 1.0, 0.0)
    vx /= vnorm; vy /= vnorm; vz /= vnorm

    ax = vy
    ay = -vx
    az = 0.0
    anorm = sqrt(ax * ax + ay * ay + az * az)
    if anorm < 1e-6
        ax, ay, az = vz, 0.0, -vx
        anorm = sqrt(ax * ax + ay * ay + az * az)
    end
    return (ax / anorm, ay / anorm, az / anorm)
end

"""
    _surface_tangent_down(domain::EllipsoidShellDomain, pos) -> NTuple{3, Float64}

Compute a surface-tangent direction that points roughly "downward" (toward the
apex, -z direction) from a point on the ellipsoid surface. Used to set initial
growth directions for coronary arteries.
"""
function _surface_tangent_down(domain::EllipsoidShellDomain, pos)
    a, b, c = domain.semi_axes
    cx, cy, cz = domain.center
    # Surface normal (gradient of x²/a² + y²/b² + z²/c²)
    nx = (pos[1] - cx) / (a * a)
    ny = (pos[2] - cy) / (b * b)
    nz = (pos[3] - cz) / (c * c)
    nlen = sqrt(nx^2 + ny^2 + nz^2)
    nlen < 1e-15 && return (0.0, 0.0, -1.0)
    nx /= nlen; ny /= nlen; nz /= nlen

    # Project "downward" (0,0,-1) onto tangent plane: d - (d·n)n
    dot_dn = -nz  # (0,0,-1) · n
    dx = -dot_dn * nx
    dy = -dot_dn * ny
    dz = -1.0 - dot_dn * nz
    dlen = sqrt(dx^2 + dy^2 + dz^2)
    dlen > 1e-10 && return (dx / dlen, dy / dlen, dz / dlen)
    return (0.0, 0.0, -1.0)
end

"""
    coronary_tree_configs(domain::EllipsoidShellDomain) -> Vector{TreeConfig}

Shell-aware coronary tree configuration. Places root positions ON the shell
surface at anatomically-motivated locations, with growth directions tangent
to the surface and pointing toward the apex.
"""
function coronary_tree_configs(domain::EllipsoidShellDomain)
    a, b, c = domain.semi_axes
    cx, cy, cz = domain.center

    # Place roots at MID-SHELL radius (not outer surface!) so the root segment
    # and subsequent bifurcation points stay inside the shell.
    # The shell spans radial fraction [1-thickness, 1.0]; mid-shell = 1 - thickness/2.
    f_mid = 1.0 - domain.thickness / 2.0

    # Spherical parameterization: (a*f sinθ cosφ, b*f sinθ sinφ, c*f cosθ)
    # θ=0 is north pole (z=+c, superior), θ=π is south pole (z=-c, apex)

    # LAD: anterior, slightly left — courses down the anterior surface
    θ_lad, φ_lad = 0.5, π / 6
    pos_lad = (cx + a * f_mid * sin(θ_lad) * cos(φ_lad),
               cy + b * f_mid * sin(θ_lad) * sin(φ_lad),
               cz + c * f_mid * cos(θ_lad))
    dir_lad = _surface_tangent_down(domain, pos_lad)

    # LCX: left-posterior — courses around the left/back
    θ_lcx, φ_lcx = 0.4, 2π / 3
    pos_lcx = (cx + a * f_mid * sin(θ_lcx) * cos(φ_lcx),
               cy + b * f_mid * sin(θ_lcx) * sin(φ_lcx),
               cz + c * f_mid * cos(θ_lcx))
    dir_lcx = _surface_tangent_down(domain, pos_lcx)

    # RCA: right-anterior — courses down the right side
    θ_rca, φ_rca = 0.5, -π / 4
    pos_rca = (cx + a * f_mid * sin(θ_rca) * cos(φ_rca),
               cy + b * f_mid * sin(θ_rca) * sin(φ_rca),
               cz + c * f_mid * cos(θ_rca))
    dir_rca = _surface_tangent_down(domain, pos_rca)

    return [
        TreeConfig("LAD", pos_lad, 1.5, dir_lad, 3000, 0.40),
        TreeConfig("LCX", pos_lcx, 1.2, dir_lcx, 2000, 0.25),
        TreeConfig("RCA", pos_rca, 1.3, dir_rca, 2500, 0.35),
    ]
end

function coronary_tree_configs(domain::CSVVolumeDomain)
    seed_specs = [
        ("LAD", (1.37769, 0.393763, 3.1823), 0.18, 3000, 0.40, +0.45),
        ("LCX", (1.36868, 0.383006, 3.00845), 0.165, 2000, 0.25, -0.45),
        ("RCA", (6.30361, 3.29052, 2.88375), 0.2, 2500, 0.35, +1.10),
    ]

    configs = TreeConfig[]
    cx, cy, cz = domain.center

    for (name, raw_root, raw_radius, terminals, territory, angle) in seed_specs
        root_guess = (
            raw_root[1] * domain.length_scale,
            raw_root[2] * domain.length_scale,
            raw_root[3] * domain.length_scale,
        )
        radius = raw_radius * domain.length_scale
        root = in_domain(domain, root_guess) ? root_guess : project_to_domain(domain, root_guess)

        dx = cx - root[1]
        dy = cy - root[2]
        dz = cz - root[3]
        dlen = sqrt(dx * dx + dy * dy + dz * dz)
        if dlen > 1e-12
            dir = (dx / dlen, dy / dlen, dz / dlen)
        else
            dir = (0.0, 0.0, -1.0)
        end

        axis = _make_perp_axis(dir)
        dir = _rotate_about_axis(dir, axis, angle)

        root_len = radius * 5.0
        trial = (
            root[1] + dir[1] * root_len,
            root[2] + dir[2] * root_len,
            root[3] + dir[3] * root_len,
        )
        if !in_domain(domain, trial)
            dir = (-dir[1], -dir[2], -dir[3])
        end

        push!(configs, TreeConfig(name, root, radius, dir, terminals, territory))
    end

    return configs
end

function fixed_tree_configs(
    trees::Dict{String, VascularTree};
    target_terminals::Dict{String, Int},
    territory_fractions::Dict{String, Float64}=Dict(
        "LAD" => 0.40,
        "LCX" => 0.25,
        "RCA" => 0.35,
    ),
)
    configs = TreeConfig[]
    for name in sort(collect(keys(trees)))
        tree = trees[name]
        root_id = tree.root_segment_id
        root_id > 0 || error("Tree `$name` has no root segment.")
        seg = tree.segments
        topo = tree.topology
        prox = (
            seg.proximal_x[root_id],
            seg.proximal_y[root_id],
            seg.proximal_z[root_id],
        )

        terminal_points = NTuple{3, Float64}[]
        for id in 1:topo.n
            topo.is_terminal[id] || continue
            push!(
                terminal_points,
                (seg.distal_x[id], seg.distal_y[id], seg.distal_z[id]),
            )
        end

        anchor =
            if isempty(terminal_points)
                (seg.distal_x[root_id], seg.distal_y[root_id], seg.distal_z[root_id])
            else
                n_terminal = length(terminal_points)
                (
                    sum(pt[1] for pt in terminal_points) / n_terminal,
                    sum(pt[2] for pt in terminal_points) / n_terminal,
                    sum(pt[3] for pt in terminal_points) / n_terminal,
                )
            end

        dx = anchor[1] - prox[1]
        dy = anchor[2] - prox[2]
        dz = anchor[3] - prox[3]
        dlen = sqrt(dx^2 + dy^2 + dz^2)
        if dlen <= 0
            dx = seg.distal_x[root_id] - prox[1]
            dy = seg.distal_y[root_id] - prox[2]
            dz = seg.distal_z[root_id] - prox[3]
            dlen = sqrt(dx^2 + dy^2 + dz^2)
        end
        dir = dlen > 0 ? (dx / dlen, dy / dlen, dz / dlen) : (1.0, 0.0, 0.0)
        terminals = get(target_terminals, name, tree.n_terminals)
        territory = get(territory_fractions, name, 1.0 / max(length(trees), 1))
        push!(configs, TreeConfig(name, prox, seg.radius[root_id], dir, terminals, territory))
    end
    return configs
end

"""
    CoronaryForest

Holds multiple vascular trees with territory partitioning.
`tree_params` maps tree name → MorphometricParams for per-artery validation.
`params` is the default/primary params (backward compatible).
"""
struct CoronaryForest
    trees::Dict{String, VascularTree}
    territory_map::TerritoryMap
    params::MorphometricParams
    tree_params::Dict{String, MorphometricParams}
end

# Backward-compatible constructor (no tree_params)
function CoronaryForest(trees, tmap, params::MorphometricParams)
    tp = Dict{String, MorphometricParams}(name => params for name in keys(trees))
    return CoronaryForest(trees, tmap, params, tp)
end

function continue_coronary_forest(
    domain::AbstractDomain,
    trees::Dict{String, VascularTree},
    params::MorphometricParams;
    tree_configs::Vector{TreeConfig},
    rng::AbstractRNG=Random.default_rng(),
    verbose::Bool=false,
    kassab::Bool=true,
)
    gamma = params.gamma
    tmap = initialize_territories(domain, tree_configs)
    domain_size = _domain_size(domain)

    buffers = Dict{String, NamedTuple}()
    grids = Dict{String, SpatialGrid}()
    last_rebuild = Dict{String, Int}()
    grid_rebuild_interval = 100

    for cfg in tree_configs
        tree = trees[cfg.name]
        cap = tree.segments.capacity
        buffers[cfg.name] = (
            distances = zeros(cap),
            costs = zeros(cap),
            bifurc_t = zeros(cap),
            intersect = Vector{Bool}(undef, cap),
        )
        cs = _grid_cell_size(domain_size, max(tree.segments.n, 10))
        grids[cfg.name] = build_grid(tree.segments, tree.segments.n, domain, cs)
        last_rebuild[cfg.name] = tree.segments.n
    end

    terminal_radius = params.diameter_mean[1] / 2.0 / 1000.0
    remaining_targets = Dict(cfg.name => max(cfg.target_terminals - trees[cfg.name].n_terminals, 0) for cfg in tree_configs)
    max_rounds = max(sum(values(remaining_targets)) * 60, 1)
    added = Dict(cfg.name => 0 for cfg in tree_configs)
    initial_terminals = Dict(cfg.name => trees[cfg.name].n_terminals for cfg in tree_configs)
    initial_segments = Dict(cfg.name => trees[cfg.name].segments.n for cfg in tree_configs)

    round = 0
    while round < max_rounds
        round += 1
        all_done = true

        for cfg in tree_configs
            remaining_targets[cfg.name] <= 0 && continue
            all_done = false

            tree = trees[cfg.name]
            buf = buffers[cfg.name]
            n = tree.segments.n

            if n - last_rebuild[cfg.name] >= grid_rebuild_interval
                cs = _grid_cell_size(domain_size, n)
                grids[cfg.name] = build_grid(tree.segments, n, domain, cs)
                last_rebuild[cfg.name] = n
            end
            grid = grids[cfg.name]

            d_thresh = domain_size / (10.0 * (n + 1)^(1.0 / 3.0))
            pt = sample_in_territory(domain, tmap, cfg.name, rng)
            pt === nothing && continue
            tx, ty, tz = pt

            collision = false
            for (other_name, other_tree) in trees
                other_name == cfg.name && continue
                if check_inter_tree_collision((tx, ty, tz), other_tree, d_thresh)
                    collision = true
                    break
                end
            end
            collision && continue

            K = min(n, 20)
            search_radius = d_thresh * 5.0
            nearest = _find_nearest_via_grid(grid, tree.segments, buf.distances, tx, ty, tz, search_radius, K, n)

            seg_idx = 0
            t_opt = 0.5
            best_cost = Inf
            found = false

            for ki in eachindex(nearest)
                ci = nearest[ki]

                ti, costi = optimize_bifurcation_point(
                    tree.segments.proximal_x[ci], tree.segments.proximal_y[ci], tree.segments.proximal_z[ci],
                    tree.segments.distal_x[ci], tree.segments.distal_y[ci], tree.segments.distal_z[ci],
                    tree.segments.radius[ci], tx, ty, tz, gamma)

                costi >= best_cost && continue

                bx = tree.segments.proximal_x[ci] + ti * (tree.segments.distal_x[ci] - tree.segments.proximal_x[ci])
                by = tree.segments.proximal_y[ci] + ti * (tree.segments.distal_y[ci] - tree.segments.proximal_y[ci])
                bz = tree.segments.proximal_z[ci] + ti * (tree.segments.distal_z[ci] - tree.segments.proximal_z[ci])

                parent_r = tree.segments.radius[ci]
                intersect_dist = max(parent_r * 2.0, d_thresh * 0.01)
                check_intersections!(buf.intersect, tree.segments, bx, by, bz, tx, ty, tz, intersect_dist, n)
                topo = tree.topology
                _clear_local_intersection_neighborhood!(buf.intersect, topo, ci)

                has_any_intersection(buf.intersect, n) && continue
                check_domain_crossing(domain, (bx, by, bz), (tx, ty, tz)) && continue

                best_cost = costi
                seg_idx = ci
                t_opt = ti
                found = true
            end

            !found && continue

            n_before = tree.segments.n
            cont_id, new_id = add_bifurcation!(tree, seg_idx, t_opt, (tx, ty, tz), terminal_radius)

            if kassab
                parent_radius = tree.segments.radius[seg_idx]
                asymmetry = sample_asymmetry(params, rng)
                r_large, r_small = compute_daughter_radii(parent_radius, asymmetry, gamma)
                tree.segments.radius[cont_id] = r_large
                tree.segments.radius[new_id] = r_small
            end

            update_radii!(tree, gamma)

            _insert_segment_to_grid!(grid, tree.segments, seg_idx)
            for new_seg_idx in (n_before + 1):tree.segments.n
                _insert_segment_to_grid!(grid, tree.segments, new_seg_idx)
            end

            added[cfg.name] += 1
            remaining_targets[cfg.name] -= 1
        end

        all_done && break
    end

    if verbose
        for cfg in tree_configs
            tree = trees[cfg.name]
            terminal_delta = tree.n_terminals - initial_terminals[cfg.name]
            segment_delta = tree.segments.n - initial_segments[cfg.name]
            println(
                "  $(cfg.name): +$(terminal_delta) terminals, +$(segment_delta) segments ",
                "(successful insertions=$(added[cfg.name])) -> total $(tree.n_terminals) terminals, $(tree.segments.n) segments",
            )
        end
    end

    return CoronaryForest(trees, tmap, params)
end

"""
    generate_coronary_forest(domain, params; tree_configs, rng) -> CoronaryForest

Generate a multi-tree forest via round-robin growth. Each tree grows in its
assigned territory, with inter-tree collision avoidance.
"""
function generate_coronary_forest(
    domain::AbstractDomain,
    params::MorphometricParams;
    tree_configs::Vector{TreeConfig}=coronary_tree_configs(domain),
    rng::AbstractRNG=Random.default_rng(),
    verbose::Bool=false,
    kassab::Bool=true,
)
    gamma = params.gamma
    tmap = initialize_territories(domain, tree_configs)

    # Create trees and add root segments
    trees = Dict{String, VascularTree}()
    capacity_per_tree = maximum(cfg.target_terminals for cfg in tree_configs) * 3
    for cfg in tree_configs
        tree = VascularTree(cfg.name, capacity_per_tree)
        # Root segment: from root_position in root_direction
        dx, dy, dz = cfg.root_direction
        dlen = sqrt(dx^2 + dy^2 + dz^2)
        if dlen > 0
            dx /= dlen; dy /= dlen; dz /= dlen
        else
            dx, dy, dz = 1.0, 0.0, 0.0
        end
        root_len = cfg.root_radius * 5.0  # root segment length proportional to radius
        distal = (
            cfg.root_position[1] + dx * root_len,
            cfg.root_position[2] + dy * root_len,
            cfg.root_position[3] + dz * root_len,
        )
        add_segment!(tree, cfg.root_position, distal, cfg.root_radius, Int32(-1))
        trees[cfg.name] = tree
    end

    # Round-robin growth
    terminal_radius = params.diameter_mean[1] / 2.0 / 1000.0
    max_rounds = maximum(cfg.target_terminals for cfg in tree_configs) * 50
    added = Dict(cfg.name => 0 for cfg in tree_configs)
    domain_size = _domain_size(domain)

    # Pre-allocate buffers per tree
    buffers = Dict{String, NamedTuple}()
    for cfg in tree_configs
        cap = capacity_per_tree
        buffers[cfg.name] = (
            distances = zeros(cap),
            costs = zeros(cap),
            bifurc_t = zeros(cap),
            intersect = Vector{Bool}(undef, cap),
        )
    end

    # Spatial grids per tree
    grids = Dict{String, SpatialGrid}()
    last_rebuild = Dict{String, Int}()
    grid_rebuild_interval = 100
    for cfg in tree_configs
        tree = trees[cfg.name]
        cs = _grid_cell_size(domain_size, max(tree.segments.n, 10))
        grids[cfg.name] = build_grid(tree.segments, tree.segments.n, domain, cs)
        last_rebuild[cfg.name] = tree.segments.n
    end

    round = 0
    while round < max_rounds
        round += 1
        all_done = true

        for cfg in tree_configs
            added[cfg.name] >= cfg.target_terminals && continue
            all_done = false

            tree = trees[cfg.name]
            buf = buffers[cfg.name]
            n = tree.segments.n

            # Rebuild grid if needed
            if n - last_rebuild[cfg.name] >= grid_rebuild_interval
                cs = _grid_cell_size(domain_size, n)
                grids[cfg.name] = build_grid(tree.segments, n, domain, cs)
                last_rebuild[cfg.name] = n
            end
            grid = grids[cfg.name]

            d_thresh = domain_size / (10.0 * (n + 1)^(1.0 / 3.0))

            # Sample in territory
            pt = sample_in_territory(domain, tmap, cfg.name, rng)
            pt === nothing && continue

            tx, ty, tz = pt

            # Check inter-tree collision
            collision = false
            for (other_name, other_tree) in trees
                other_name == cfg.name && continue
                if check_inter_tree_collision((tx, ty, tz), other_tree, d_thresh)
                    collision = true
                    break
                end
            end
            collision && continue

            # Find best connection using spatial grid
            K = min(n, 20)
            search_radius = d_thresh * 5.0
            nearest = _find_nearest_via_grid(grid, tree.segments, buf.distances, tx, ty, tz, search_radius, K, n)

            seg_idx = 0
            t_opt = 0.5
            best_cost = Inf
            found = false

            for ki in 1:length(nearest)
                ci = nearest[ki]

                ti, costi = optimize_bifurcation_point(
                    tree.segments.proximal_x[ci], tree.segments.proximal_y[ci], tree.segments.proximal_z[ci],
                    tree.segments.distal_x[ci], tree.segments.distal_y[ci], tree.segments.distal_z[ci],
                    tree.segments.radius[ci], tx, ty, tz, gamma)

                costi >= best_cost && continue

                bx = tree.segments.proximal_x[ci] + ti * (tree.segments.distal_x[ci] - tree.segments.proximal_x[ci])
                by = tree.segments.proximal_y[ci] + ti * (tree.segments.distal_y[ci] - tree.segments.proximal_y[ci])
                bz = tree.segments.proximal_z[ci] + ti * (tree.segments.distal_z[ci] - tree.segments.proximal_z[ci])

                parent_r = tree.segments.radius[ci]
                intersect_dist = max(parent_r * 2.0, d_thresh * 0.01)
                check_intersections!(buf.intersect, tree.segments, bx, by, bz, tx, ty, tz, intersect_dist, n)

                topo = tree.topology
                _clear_local_intersection_neighborhood!(buf.intersect, topo, ci)

                has_any_intersection(buf.intersect, n) && continue
                check_domain_crossing(domain, (bx, by, bz), (tx, ty, tz)) && continue

                best_cost = costi
                seg_idx = ci
                t_opt = ti
                found = true
            end

            !found && continue

            n_before = tree.segments.n

            # Add bifurcation
            cont_id, new_id = add_bifurcation!(tree, seg_idx, t_opt, (tx, ty, tz), terminal_radius)

            if kassab
                parent_radius = tree.segments.radius[seg_idx]
                asymmetry = sample_asymmetry(params, rng)
                r_large, r_small = compute_daughter_radii(parent_radius, asymmetry, gamma)
                tree.segments.radius[cont_id] = r_large
                tree.segments.radius[new_id] = r_small
            end

            update_radii!(tree, gamma)

            # Insert new segments into grid
            _insert_segment_to_grid!(grid, tree.segments, seg_idx)
            for new_seg_idx in (n_before + 1):tree.segments.n
                _insert_segment_to_grid!(grid, tree.segments, new_seg_idx)
            end

            added[cfg.name] += 1
        end

        all_done && break
    end

    if verbose
        for cfg in tree_configs
            println("  $(cfg.name): $(added[cfg.name]) / $(cfg.target_terminals) terminals, $(trees[cfg.name].segments.n) segments")
        end
    end

    return CoronaryForest(trees, tmap, params)
end

"""
    _project_tree_to_domain!(tree, domain)

Project all segment endpoints back into the domain. Applied after Barabasi
geometry to ensure no segments escape the domain boundary.
Processes in topological order (parent before child) so that proximal
adjustments propagate correctly.
"""
function _project_tree_to_domain!(tree::VascularTree, domain::EllipsoidShellDomain)
    seg = tree.segments
    topo = tree.topology

    # Multiple passes in topological order. For out-of-domain distals, redirect the
    # segment to be surface-tangent (preserving length) rather than just projecting
    # the endpoint. This keeps branches following the shell surface.
    for _pass in 1:5
        n_fixed = 0
        for i in 1:seg.n
            # Project proximal if outside (root segments, etc.)
            pp = (seg.proximal_x[i], seg.proximal_y[i], seg.proximal_z[i])
            if !in_domain(domain, pp)
                proj = project_to_domain(domain, pp)
                seg.proximal_x[i] = proj[1]
                seg.proximal_y[i] = proj[2]
                seg.proximal_z[i] = proj[3]
                pp = proj
                n_fixed += 1
            end

            dp = (seg.distal_x[i], seg.distal_y[i], seg.distal_z[i])
            if !in_domain(domain, dp)
                # Redirect: project direction to surface tangent, re-extend
                dx = dp[1] - pp[1]
                dy = dp[2] - pp[2]
                dz = dp[3] - pp[3]
                seg_len = sqrt(dx^2 + dy^2 + dz^2)

                if seg_len > 1e-15
                    dir = (dx / seg_len, dy / seg_len, dz / seg_len)
                    tdir = _shell_tangent_project(domain, pp, dir)
                    new_dp = (pp[1] + tdir[1] * seg_len,
                              pp[2] + tdir[2] * seg_len,
                              pp[3] + tdir[3] * seg_len)

                    # If still outside after tangent redirect, project endpoint
                    if !in_domain(domain, new_dp)
                        new_dp = project_to_domain(domain, new_dp)
                    end
                else
                    new_dp = project_to_domain(domain, dp)
                end

                seg.distal_x[i] = new_dp[1]
                seg.distal_y[i] = new_dp[2]
                seg.distal_z[i] = new_dp[3]
                n_fixed += 1

                # Sync children's proximal points
                for cid in (topo.child1_id[i], topo.child2_id[i], topo.child3_id[i])
                    if cid > 0
                        seg.proximal_x[cid] = new_dp[1]
                        seg.proximal_y[cid] = new_dp[2]
                        seg.proximal_z[cid] = new_dp[3]
                    end
                end
            end

            # Update segment length
            ex = seg.distal_x[i] - seg.proximal_x[i]
            ey = seg.distal_y[i] - seg.proximal_y[i]
            ez = seg.distal_z[i] - seg.proximal_z[i]
            seg.seg_length[i] = sqrt(ex * ex + ey * ey + ez * ez)
        end
        n_fixed == 0 && break  # converged
    end
end

function _project_tree_to_domain!(tree::VascularTree, domain::AbstractDomain)
    seg = tree.segments
    topo = tree.topology

    for _pass in 1:3
        n_fixed = 0
        for i in 1:seg.n
            pp = (seg.proximal_x[i], seg.proximal_y[i], seg.proximal_z[i])
            if !in_domain(domain, pp)
                proj = project_to_domain(domain, pp)
                seg.proximal_x[i] = proj[1]
                seg.proximal_y[i] = proj[2]
                seg.proximal_z[i] = proj[3]
                pp = proj
                n_fixed += 1
            end

            dp = (seg.distal_x[i], seg.distal_y[i], seg.distal_z[i])
            if !in_domain(domain, dp)
                proj = project_to_domain(domain, dp)
                seg.distal_x[i] = proj[1]
                seg.distal_y[i] = proj[2]
                seg.distal_z[i] = proj[3]
                n_fixed += 1

                for cid in (topo.child1_id[i], topo.child2_id[i], topo.child3_id[i])
                    if cid > 0
                        seg.proximal_x[cid] = proj[1]
                        seg.proximal_y[cid] = proj[2]
                        seg.proximal_z[cid] = proj[3]
                    end
                end
            end

            ex = seg.distal_x[i] - seg.proximal_x[i]
            ey = seg.distal_y[i] - seg.proximal_y[i]
            ez = seg.distal_z[i] - seg.proximal_z[i]
            seg.seg_length[i] = sqrt(ex * ex + ey * ey + ez * ez)
        end
        n_fixed == 0 && break
    end
    return tree
end

"""
    _effective_domain_volume(domain) -> Float64

Compute the actual volume of the domain. For shell domains, this is the volume
of the shell (not the bounding box). Used by the OpenCCO d_thresh formula.
"""
function _effective_domain_volume(domain::AbstractDomain)
    if domain isa EllipsoidShellDomain
        a, b, c = domain.semi_axes
        inner = 1.0 - domain.thickness
        return (4.0 / 3.0) * Float64(π) * a * b * c * (1.0 - inner^3)
    elseif domain isa EllipsoidDomain
        a, b, c = domain.semi_axes
        return (4.0 / 3.0) * Float64(π) * a * b * c
    elseif domain isa SphereDomain
        return (4.0 / 3.0) * Float64(π) * domain.radius^3
    elseif domain isa BoxDomain
        dx = domain.max_corner[1] - domain.min_corner[1]
        dy = domain.max_corner[2] - domain.min_corner[2]
        dz = domain.max_corner[3] - domain.min_corner[3]
        return dx * dy * dz
    elseif domain isa CSVVolumeDomain
        return domain.volume
    elseif domain isa CSVShellDomain
        return domain.volume
    end
    s = _domain_size(domain)
    return s^3
end

"""
    _params_for_tree(name, default_params) -> MorphometricParams

Auto-select per-artery Kassab params based on tree name.
Falls back to `default_params` for unrecognized names.
"""
function _params_for_tree(name::String, default_params::MorphometricParams)
    upper = uppercase(name)
    if contains(upper, "RCA")
        return kassab_rca_params()
    elseif contains(upper, "LAD")
        return kassab_lad_params()
    elseif contains(upper, "LCX")
        return kassab_lcx_params()
    else
        return default_params
    end
end

"""
    generate_kassab_coronary(domain, params; rng, verbose, handoff_order, tree_configs) -> CoronaryForest

Full pipeline: CCO skeleton + statistical subdivision + Barabasi geometry.

1. **CCO phase**: Grid-accelerated multi-tree growth with territory partitioning.
2. **Subdivision phase**: subdivide_terminals! per tree using per-artery Kassab CM.
   Radii and asymmetry are baked into subdivision (CM-implied asymmetry).
3. **Refinement phase**: Barabasi junction geometry only (no post-hoc radius assignment).
4. **Summary**: per-tree segment counts.
"""
function generate_kassab_coronary(
    domain::AbstractDomain,
    params::MorphometricParams;
    rng::AbstractRNG=Random.default_rng(),
    verbose::Bool=false,
    handoff_order::Int=5,
    tree_configs::Vector{TreeConfig}=coronary_tree_configs(domain),
)
    t_start = time()

    # Build per-artery params map
    per_artery = Dict{String, MorphometricParams}()
    for cfg in tree_configs
        per_artery[cfg.name] = _params_for_tree(cfg.name, params)
    end

    # Estimate capacity: CCO skeleton + subdivision expansion
    max_terminals = maximum(cfg.target_terminals for cfg in tree_configs)
    skeleton_capacity = max_terminals * 3
    expansion = estimate_total_segments(handoff_order, params)
    total_capacity = min(skeleton_capacity * expansion * 5, 10_000_000)

    cco_configs = [
        TreeConfig(cfg.name, cfg.root_position, cfg.root_radius, cfg.root_direction,
                   cfg.target_terminals, cfg.territory_fraction)
        for cfg in tree_configs
    ]

    tmap = initialize_territories(domain, cco_configs)
    gamma = params.gamma

    trees = Dict{String, VascularTree}()
    for cfg in cco_configs
        tree = VascularTree(cfg.name, total_capacity)
        dx, dy, dz = cfg.root_direction
        dlen = sqrt(dx^2 + dy^2 + dz^2)
        if dlen > 0
            dx /= dlen; dy /= dlen; dz /= dlen
        else
            dx, dy, dz = 1.0, 0.0, 0.0
        end
        root_len = cfg.root_radius * 5.0
        distal = (
            cfg.root_position[1] + dx * root_len,
            cfg.root_position[2] + dy * root_len,
            cfg.root_position[3] + dz * root_len,
        )
        # Ensure distal point stays inside the domain (tangent lines exit
        # curved shells quickly). Project back if needed.
        if !in_domain(domain, distal)
            distal = project_to_domain(domain, distal)
        end
        add_segment!(tree, cfg.root_position, distal, cfg.root_radius, Int32(-1))
        trees[cfg.name] = tree
    end

    # Phase 1: CCO skeleton growth
    # Uses OpenCCO-inspired progressive domain expansion (IPOL Eq. 14-15):
    #   - d_thresh from domain volume and target terminal count (not ad-hoc)
    #   - search_radius GROWS as tree expands (progressive outward growth)
    #   - min-distance check prevents degenerate short branches
    t_phase1 = time()
    if verbose
        println("Phase 1: CCO skeleton growth...")
    end

    # Per-artery terminal radius at handoff order
    terminal_radii = Dict{String, Float64}()
    for cfg in cco_configs
        tp = per_artery[cfg.name]
        ho = min(handoff_order, tp.n_orders - 1)
        terminal_radii[cfg.name] = tp.diameter_mean[ho + 1] / 2.0 / 1000.0
    end
    max_rounds = max_terminals * 50
    added = Dict(cfg.name => 0 for cfg in cco_configs)
    domain_size = _domain_size(domain)

    # Effective domain volume for OpenCCO d_thresh formula
    domain_vol = _effective_domain_volume(domain)

    # Use CCO-sized buffers (not total_capacity) to save memory during CCO phase
    cco_buf_cap = skeleton_capacity
    buffers = Dict{String, NamedTuple}()
    grids = Dict{String, SpatialGrid}()
    last_rebuild = Dict{String, Int}()
    grid_rebuild_interval = 100

    # Pre-allocate inter-tree collision buffers (avoid per-call allocation)
    collision_bufs = Dict{String, Vector{Float64}}()

    for cfg in cco_configs
        buffers[cfg.name] = (
            distances = zeros(cco_buf_cap),
            costs = zeros(cco_buf_cap),
            bifurc_t = zeros(cco_buf_cap),
            intersect = Vector{Bool}(undef, cco_buf_cap),
        )
        collision_bufs[cfg.name] = zeros(cco_buf_cap)
        tree = trees[cfg.name]
        cs = _grid_cell_size(domain_size, max(tree.segments.n, 10))
        grids[cfg.name] = build_grid(tree.segments, tree.segments.n, domain, cs)
        last_rebuild[cfg.name] = tree.segments.n
    end

    iter = 0
    while iter < max_rounds
        iter += 1
        all_done = true

        for cfg in cco_configs
            added[cfg.name] >= cfg.target_terminals && continue
            all_done = false

            tree = trees[cfg.name]
            buf = buffers[cfg.name]
            n = tree.segments.n

            if n - last_rebuild[cfg.name] >= grid_rebuild_interval
                cs = _grid_cell_size(domain_size, n)
                grids[cfg.name] = build_grid(tree.segments, n, domain, cs)
                last_rebuild[cfg.name] = n
            end
            grid = grids[cfg.name]

            n_added = added[cfg.name]

            # OpenCCO Eq. 15: d_thresh = (V / (k * (i+1)))^(1/3)
            # Per-tree territory volume ≈ territory_fraction * domain_vol
            territory_vol = cfg.territory_fraction * domain_vol
            d_thresh = (territory_vol / (cfg.target_terminals * max(1, n_added + 1)))^(1.0 / 3.0)

            # Search radius for spatial grid query (find K nearest segments).
            # Use full domain size so CCO can connect to any reachable segment.
            search_radius = domain_size * 2.0

            pt = sample_in_territory(domain, tmap, cfg.name, rng)
            pt === nothing && continue

            tx, ty, tz = pt

            # Distance check against own tree: prevents degenerate short branches.
            # No upper-bound rejection — let CCO fill the entire domain naturally.
            # The K-nearest cost optimization handles connection quality.
            if n > 0
                compute_all_distances!(buf.distances, tree.segments, tx, ty, tz, n)
                own_min_dist = AK.minimum(@view buf.distances[1:n])
                own_min_dist < d_thresh && continue
            end

            # Inter-tree collision check (using pre-allocated buffers)
            collision = false
            for (other_name, other_tree) in trees
                other_name == cfg.name && continue
                on = other_tree.segments.n
                on == 0 && continue
                cbuf = collision_bufs[other_name]
                compute_all_distances!(cbuf, other_tree.segments, tx, ty, tz, on)
                other_min = AK.minimum(@view cbuf[1:on])
                if other_min < d_thresh
                    collision = true
                    break
                end
            end
            collision && continue

            K = min(n, 40)
            nearest = _find_nearest_via_grid(grid, tree.segments, buf.distances, tx, ty, tz, search_radius, K, n)

            seg_idx = 0
            t_opt = 0.5
            best_cost = Inf
            found = false

            for ki in 1:length(nearest)
                ci = nearest[ki]

                ti, costi = optimize_bifurcation_point(
                    tree.segments.proximal_x[ci], tree.segments.proximal_y[ci], tree.segments.proximal_z[ci],
                    tree.segments.distal_x[ci], tree.segments.distal_y[ci], tree.segments.distal_z[ci],
                    tree.segments.radius[ci], tx, ty, tz, gamma)

                costi >= best_cost && continue

                bx = tree.segments.proximal_x[ci] + ti * (tree.segments.distal_x[ci] - tree.segments.proximal_x[ci])
                by = tree.segments.proximal_y[ci] + ti * (tree.segments.distal_y[ci] - tree.segments.proximal_y[ci])
                bz = tree.segments.proximal_z[ci] + ti * (tree.segments.distal_z[ci] - tree.segments.proximal_z[ci])

                parent_r = tree.segments.radius[ci]
                intersect_dist = max(parent_r * 2.0, d_thresh * 0.01)
                check_intersections!(buf.intersect, tree.segments, bx, by, bz, tx, ty, tz, intersect_dist, n)
                topo = tree.topology
                _clear_local_intersection_neighborhood!(buf.intersect, topo, ci)

                has_any_intersection(buf.intersect, n) && continue
                check_domain_crossing(domain, (bx, by, bz), (tx, ty, tz)) && continue

                best_cost = costi
                seg_idx = ci
                t_opt = ti
                found = true
            end

            !found && continue

            n_before = tree.segments.n
            cont_id, new_id = add_bifurcation!(tree, seg_idx, t_opt, (tx, ty, tz), terminal_radii[cfg.name])
            # Skip update_radii! during CCO — deferred to end of phase (radii only
            # matter for cost computation which uses parent radius, unaffected by Murray)

            _insert_segment_to_grid!(grid, tree.segments, seg_idx)
            for new_seg_idx in (n_before + 1):tree.segments.n
                _insert_segment_to_grid!(grid, tree.segments, new_seg_idx)
            end

            added[cfg.name] += 1
        end

        all_done && break
    end

    # Single Murray's law pass after all CCO growth complete
    for (name, tree) in trees
        update_radii!(tree, gamma)
    end

    t_phase1_end = time()
    if verbose
        for cfg in cco_configs
            println("  $(cfg.name): $(added[cfg.name]) / $(cfg.target_terminals) CCO terminals, $(trees[cfg.name].segments.n) skeleton segments")
        end
        println("  Phase 1 time: $(round(t_phase1_end - t_phase1, digits=1))s")
    end

    # Free CCO buffers before subdivision (no longer needed)
    empty!(buffers)

    # Phase 2: Statistical subdivision
    t_phase2 = time()
    if verbose
        println("Phase 2: Statistical subdivision...")
    end

    for (name, tree) in trees
        tp = per_artery[name]
        ho = min(handoff_order, tp.n_orders - 1)
        subdivide_terminals!(tree, tp; rng=rng, max_order=ho, domain=domain)
        if verbose
            println("  $(name): $(tree.segments.n) segments after subdivision")
        end
    end

    # Fix low-order terminal radii: set order-0 and order-1 terminals to Kassab
    # diameters, then propagate Murray's law upward. This corrects Murray-chain
    # attenuation at low orders (which produces ~2.7um instead of 8um) while
    # keeping Murray's law exact at every junction.
    for (name, tree) in trees
        tp = per_artery[name]
        seg = tree.segments
        topo = tree.topology
        for i in 1:seg.n
            if topo.is_terminal[i]
                ord = Int(topo.strahler_order[i])
                if ord <= 1 && ord + 1 <= length(tp.diameter_mean_elem)
                    seg.radius[i] = tp.diameter_mean_elem[ord + 1] / 2.0 / 1000.0
                end
            end
        end
        update_radii!(tree, gamma)
    end

    t_phase2_end = time()
    if verbose
        println("  Phase 2 time: $(round(t_phase2_end - t_phase2, digits=1))s")
    end

    # Phase 3: Barabasi junction geometry (radii already baked into subdivision)
    t_phase3 = time()
    if verbose
        println("Phase 3: Kassab radius refinement + Barabasi geometry...")
    end

    for (name, tree) in trees
        tp = per_artery[name]

        # Apply Barabasi junction geometry in topological (top-down) order
        # sequential: tree traversal for correct parent-before-child ordering
        topo_ord = _topo_order(tree)
        for k in 1:length(topo_ord)
            i = topo_ord[k]
            if tree.topology.junction_type[i] == :bifurcation
                apply_junction_geometry!(tree, Int32(i), tp)
            elseif tree.topology.junction_type[i] == :trifurcation
                apply_trifurcation_geometry!(tree, Int32(i), tp)
            end
        end

        if verbose
            println("  $(name): radii + geometry applied")
        end
    end

    t_phase3_end = time()
    if verbose
        println("  Phase 3 time: $(round(t_phase3_end - t_phase3, digits=1))s")
    end

    # Phase 3b: Domain projection — project any out-of-domain endpoints back
    # Barabasi geometry may have rotated segments outside the domain
    if domain isa EllipsoidShellDomain
        for (name, tree) in trees
            _project_tree_to_domain!(tree, domain)
        end
        if verbose
            println("  Domain projection applied")
        end
    end

    # Phase 4: Summary
    forest = CoronaryForest(trees, tmap, params, per_artery)
    t_total = time() - t_start

    if verbose
        println("Phase 4: Summary")
        for (name, tree) in sort(collect(forest.trees), by=x->x[1])
            println("  $(name): $(tree.segments.n) segments")
        end
        total = sum(t.segments.n for (_, t) in trees)
        println("  Total: $total segments in $(round(t_total, digits=1))s")
    end

    return forest
end

"""
    validate_forest(forest) -> Dict{Symbol, Any}
    validate_forest(forest, params) -> Dict{Symbol, Any}

Validate a multi-tree forest using per-tree params from the forest.
Returns a dict with per-tree and inter-tree metrics.
"""
function validate_forest(forest::CoronaryForest)
    report = Dict{Symbol, Any}()

    total_segments = 0
    tree_reports = Dict{String, ValidationReport}()

    for (name, tree) in forest.trees
        total_segments += tree.segments.n
        tp = get(forest.tree_params, name, forest.params)
        tree_reports[name] = validate_tree(tree, tp)
    end

    report[:total_segments] = total_segments
    report[:tree_reports] = tree_reports
    report[:n_trees] = length(forest.trees)

    return report
end

# Backward-compatible overload
validate_forest(forest::CoronaryForest, params::MorphometricParams) = validate_forest(forest)
