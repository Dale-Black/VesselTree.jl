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

"""
    initialize_territories(domain, configs) -> TerritoryMap

Create a territory map by assigning each voxel to the nearest tree root.
"""
function initialize_territories(domain::AbstractDomain, configs::Vector{TreeConfig})
    # Compute bounding box
    if domain isa SphereDomain
        lo = (domain.center[1] - domain.radius, domain.center[2] - domain.radius, domain.center[3] - domain.radius)
        hi = (domain.center[1] + domain.radius, domain.center[2] + domain.radius, domain.center[3] + domain.radius)
    elseif domain isa BoxDomain
        lo = domain.min_corner
        hi = domain.max_corner
    elseif domain isa EllipsoidDomain
        lo = (domain.center[1] - domain.semi_axes[1], domain.center[2] - domain.semi_axes[2], domain.center[3] - domain.semi_axes[3])
        hi = (domain.center[1] + domain.semi_axes[1], domain.center[2] + domain.semi_axes[2], domain.center[3] + domain.semi_axes[3])
    elseif domain isa EllipsoidShellDomain
        lo = (domain.center[1] - domain.semi_axes[1], domain.center[2] - domain.semi_axes[2], domain.center[3] - domain.semi_axes[3])
        hi = (domain.center[1] + domain.semi_axes[1], domain.center[2] + domain.semi_axes[2], domain.center[3] + domain.semi_axes[3])
    else
        error("Unsupported domain type")
    end

    # Cell size: ~20 cells per dimension
    extent = max(hi[1] - lo[1], hi[2] - lo[2], hi[3] - lo[3])
    cell_size = extent / 20.0

    nx = max(1, ceil(Int, (hi[1] - lo[1]) / cell_size))
    ny = max(1, ceil(Int, (hi[2] - lo[2]) / cell_size))
    nz = max(1, ceil(Int, (hi[3] - lo[3]) / cell_size))

    tree_names = [cfg.name for cfg in configs]
    n_cells = nx * ny * nz
    cells = zeros(Int, n_cells)

    # Assign each cell to nearest root
    for iz in 1:nz
        for iy in 1:ny
            for ix in 1:nx
                cx = lo[1] + (ix - 0.5) * cell_size
                cy = lo[2] + (iy - 0.5) * cell_size
                cz = lo[3] + (iz - 0.5) * cell_size

                best_idx = 1
                best_dist = Inf
                for (ti, cfg) in enumerate(configs)
                    dx = cx - cfg.root_position[1]
                    dy = cy - cfg.root_position[2]
                    dz = cz - cfg.root_position[3]
                    d = sqrt(dx^2 + dy^2 + dz^2)
                    if d < best_dist
                        best_dist = d
                        best_idx = ti
                    end
                end

                flat_idx = ix + (iy - 1) * nx + (iz - 1) * nx * ny
                cells[flat_idx] = best_idx
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

    # Place roots on the upper portion of the outer ellipsoid surface
    # using spherical parameterization: (a sinθ cosφ, b sinθ sinφ, c cosθ)
    # θ=0 is north pole (z=+c, superior), θ=π is south pole (z=-c, apex)

    # LAD: anterior, slightly left — courses down the anterior surface
    # θ≈0.5 (near base), φ≈π/6 (anterior-left)
    θ_lad, φ_lad = 0.5, π / 6
    pos_lad = (cx + a * sin(θ_lad) * cos(φ_lad),
               cy + b * sin(θ_lad) * sin(φ_lad),
               cz + c * cos(θ_lad))
    dir_lad = _surface_tangent_down(domain, pos_lad)

    # LCX: left-posterior — courses around the left/back
    # θ≈0.4 (near base), φ≈2π/3 (left-posterior)
    θ_lcx, φ_lcx = 0.4, 2π / 3
    pos_lcx = (cx + a * sin(θ_lcx) * cos(φ_lcx),
               cy + b * sin(θ_lcx) * sin(φ_lcx),
               cz + c * cos(θ_lcx))
    dir_lcx = _surface_tangent_down(domain, pos_lcx)

    # RCA: right-anterior — courses down the right side
    # θ≈0.5 (near base), φ≈-π/4 (right-anterior)
    θ_rca, φ_rca = 0.5, -π / 4
    pos_rca = (cx + a * sin(θ_rca) * cos(φ_rca),
               cy + b * sin(θ_rca) * sin(φ_rca),
               cz + c * cos(θ_rca))
    dir_rca = _surface_tangent_down(domain, pos_rca)

    return [
        TreeConfig("LAD", pos_lad, 1.5, dir_lad, 500, 0.40),
        TreeConfig("LCX", pos_lcx, 1.2, dir_lcx, 300, 0.25),
        TreeConfig("RCA", pos_rca, 1.3, dir_rca, 400, 0.35),
    ]
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

"""
    generate_coronary_forest(domain, params; tree_configs, rng) -> CoronaryForest

Generate a multi-tree forest via round-robin growth. Each tree grows in its
assigned territory, with inter-tree collision avoidance.
"""
function generate_coronary_forest(
    domain::AbstractDomain,
    params::MorphometricParams;
    tree_configs::Vector{TreeConfig}=coronary_tree_configs(),
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

                # Exclude topological neighbors
                buf.intersect[ci] = false
                topo = tree.topology
                pp = topo.parent_id[ci]
                if pp > 0
                    buf.intersect[pp] = false
                    for c in (topo.child1_id[pp], topo.child2_id[pp], topo.child3_id[pp])
                        c > 0 && (buf.intersect[c] = false)
                    end
                end
                for c in (topo.child1_id[ci], topo.child2_id[ci], topo.child3_id[ci])
                    c > 0 && (buf.intersect[c] = false)
                end

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
    tree_configs::Vector{TreeConfig}=(domain isa EllipsoidShellDomain ? coronary_tree_configs(domain) : coronary_tree_configs()),
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
        add_segment!(tree, cfg.root_position, distal, cfg.root_radius, Int32(-1))
        trees[cfg.name] = tree
    end

    # Phase 1: CCO skeleton growth
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

    # Use CCO-sized buffers (not total_capacity) to save memory during CCO phase
    cco_buf_cap = skeleton_capacity
    buffers = Dict{String, NamedTuple}()
    grids = Dict{String, SpatialGrid}()
    last_rebuild = Dict{String, Int}()
    grid_rebuild_interval = 100

    for cfg in cco_configs
        buffers[cfg.name] = (
            distances = zeros(cco_buf_cap),
            costs = zeros(cco_buf_cap),
            bifurc_t = zeros(cco_buf_cap),
            intersect = Vector{Bool}(undef, cco_buf_cap),
        )
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

                buf.intersect[ci] = false
                topo = tree.topology
                pp = topo.parent_id[ci]
                if pp > 0
                    buf.intersect[pp] = false
                    for c in (topo.child1_id[pp], topo.child2_id[pp], topo.child3_id[pp])
                        c > 0 && (buf.intersect[c] = false)
                    end
                end
                for c in (topo.child1_id[ci], topo.child2_id[ci], topo.child3_id[ci])
                    c > 0 && (buf.intersect[c] = false)
                end

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
