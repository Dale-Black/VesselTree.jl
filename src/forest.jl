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
    CoronaryForest

Holds multiple vascular trees with territory partitioning.
"""
struct CoronaryForest
    trees::Dict{String, VascularTree}
    territory_map::TerritoryMap
    params::MorphometricParams
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

            # Domain size for adaptive thresholding
            domain_size = domain isa SphereDomain ? domain.radius : 10.0

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

            # Find best connection
            evaluate_all_connections!(buf.costs, buf.bifurc_t, tree.segments, n, tx, ty, tz, gamma)
            seg_idx, t_opt, _ = select_best_connection(buf.costs, buf.bifurc_t, n)

            # Bifurcation point
            bx = tree.segments.proximal_x[seg_idx] + t_opt * (tree.segments.distal_x[seg_idx] - tree.segments.proximal_x[seg_idx])
            by = tree.segments.proximal_y[seg_idx] + t_opt * (tree.segments.distal_y[seg_idx] - tree.segments.proximal_y[seg_idx])
            bz = tree.segments.proximal_z[seg_idx] + t_opt * (tree.segments.distal_z[seg_idx] - tree.segments.proximal_z[seg_idx])

            # Self-intersection check
            check_intersections!(buf.intersect, tree.segments, bx, by, bz, tx, ty, tz, d_thresh * 0.1, n)
            buf.intersect[seg_idx] = false
            has_any_intersection(buf.intersect, n) && continue

            # Domain crossing
            check_domain_crossing(domain, (bx, by, bz), (tx, ty, tz)) && continue

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
    validate_forest(forest, params) -> Dict{Symbol, Any}

Validate a multi-tree forest. Returns a dict with per-tree and inter-tree metrics.
"""
function validate_forest(forest::CoronaryForest, params::MorphometricParams)
    report = Dict{Symbol, Any}()

    total_segments = 0
    tree_reports = Dict{String, ValidationReport}()

    for (name, tree) in forest.trees
        total_segments += tree.segments.n
        tree_reports[name] = validate_tree(tree, params)
    end

    report[:total_segments] = total_segments
    report[:tree_reports] = tree_reports
    report[:n_trees] = length(forest.trees)

    return report
end
