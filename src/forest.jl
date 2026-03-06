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
