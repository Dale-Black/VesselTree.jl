# Spatial indexing for fast segment queries

"""
    SpatialGrid

3D uniform grid for bucketing segment midpoints. Enables O(1) lookup of
nearby segments instead of O(n) full scan.
"""
mutable struct SpatialGrid
    cell_size::Float64
    origin::NTuple{3, Float64}
    dims::NTuple{3, Int}           # grid dimensions (nx, ny, nz)
    cells::Dict{Int, Vector{Int}}  # cell_index -> segment indices
end

"""
    cell_index(grid, x, y, z) -> Int

Compute the flat cell index for a 3D point.
"""
@inline function cell_index(grid::SpatialGrid, x::Float64, y::Float64, z::Float64)
    ix = clamp(floor(Int, (x - grid.origin[1]) / grid.cell_size) + 1, 1, grid.dims[1])
    iy = clamp(floor(Int, (y - grid.origin[2]) / grid.cell_size) + 1, 1, grid.dims[2])
    iz = clamp(floor(Int, (z - grid.origin[3]) / grid.cell_size) + 1, 1, grid.dims[3])
    return ix + (iy - 1) * grid.dims[1] + (iz - 1) * grid.dims[1] * grid.dims[2]
end

"""
    build_grid(segments::SegmentData, n::Int, domain, cell_size::Float64) -> SpatialGrid

Build a spatial grid from the first `n` segments, bucketing by midpoint.
"""
function build_grid(segments::SegmentData, n::Int, domain::AbstractDomain, cell_size::Float64)
    # Compute bounding box from domain
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
        error("Unsupported domain type for spatial grid")
    end

    nx = max(1, ceil(Int, (hi[1] - lo[1]) / cell_size))
    ny = max(1, ceil(Int, (hi[2] - lo[2]) / cell_size))
    nz = max(1, ceil(Int, (hi[3] - lo[3]) / cell_size))

    grid = SpatialGrid(cell_size, lo, (nx, ny, nz), Dict{Int, Vector{Int}}())

    for i in 1:n
        mx = (segments.proximal_x[i] + segments.distal_x[i]) / 2.0
        my = (segments.proximal_y[i] + segments.distal_y[i]) / 2.0
        mz = (segments.proximal_z[i] + segments.distal_z[i]) / 2.0
        idx = cell_index(grid, mx, my, mz)
        if !haskey(grid.cells, idx)
            grid.cells[idx] = Int[]
        end
        push!(grid.cells[idx], i)
    end

    return grid
end

"""
    insert!(grid::SpatialGrid, seg_idx::Int, mx, my, mz)

Insert a segment into the grid by its midpoint coordinates.
"""
function insert!(grid::SpatialGrid, seg_idx::Int, mx::Float64, my::Float64, mz::Float64)
    idx = cell_index(grid, mx, my, mz)
    if !haskey(grid.cells, idx)
        grid.cells[idx] = Int[]
    end
    push!(grid.cells[idx], seg_idx)
    return nothing
end

"""
    query_nearby(grid::SpatialGrid, x, y, z, search_radius) -> Vector{Int}

Return indices of all segments whose midpoint cell is within `search_radius`
of the query point. Searches a neighborhood of cells.
"""
function query_nearby(grid::SpatialGrid, x::Float64, y::Float64, z::Float64, search_radius::Float64)
    result = Int[]
    n_cells = ceil(Int, search_radius / grid.cell_size)

    cx = floor(Int, (x - grid.origin[1]) / grid.cell_size) + 1
    cy = floor(Int, (y - grid.origin[2]) / grid.cell_size) + 1
    cz = floor(Int, (z - grid.origin[3]) / grid.cell_size) + 1

    for dz in -n_cells:n_cells
        iz = cz + dz
        (iz < 1 || iz > grid.dims[3]) && continue
        for dy in -n_cells:n_cells
            iy = cy + dy
            (iy < 1 || iy > grid.dims[2]) && continue
            for dx in -n_cells:n_cells
                ix = cx + dx
                (ix < 1 || ix > grid.dims[1]) && continue
                idx = ix + (iy - 1) * grid.dims[1] + (iz - 1) * grid.dims[1] * grid.dims[2]
                if haskey(grid.cells, idx)
                    append!(result, grid.cells[idx])
                end
            end
        end
    end

    return result
end
