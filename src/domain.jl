# Perfusion domain types for vascular tree generation

using DelimitedFiles

"""
    AbstractDomain

Interface for perfusion domains. Subtypes must implement:
- `in_domain(domain, point) -> Bool`
- `sample_point(domain, rng) -> NTuple{3, Float64}`
- `signed_distance(domain, point) -> Float64` (negative inside, positive outside)
"""
abstract type AbstractDomain end

"""
    domain_bounds(domain) -> (min_corner, max_corner)

Axis-aligned bounding box for a perfusion domain.
"""
function domain_bounds end

struct PointCloudGrid
    cell_size::Float64
    origin::NTuple{3, Float64}
    dims::NTuple{3, Int}
    cells::Dict{Int, Vector{Int}}
end

# --- SphereDomain ---

struct SphereDomain <: AbstractDomain
    center::NTuple{3, Float64}
    radius::Float64
end

domain_bounds(d::SphereDomain) = (
    (d.center[1] - d.radius, d.center[2] - d.radius, d.center[3] - d.radius),
    (d.center[1] + d.radius, d.center[2] + d.radius, d.center[3] + d.radius),
)

function in_domain(d::SphereDomain, point)
    dx = point[1] - d.center[1]
    dy = point[2] - d.center[2]
    dz = point[3] - d.center[3]
    return dx * dx + dy * dy + dz * dz <= d.radius * d.radius
end

function signed_distance(d::SphereDomain, point)
    dx = point[1] - d.center[1]
    dy = point[2] - d.center[2]
    dz = point[3] - d.center[3]
    return sqrt(dx * dx + dy * dy + dz * dz) - d.radius
end

function sample_point(d::SphereDomain, rng::AbstractRNG)
    while true
        x = d.center[1] + d.radius * (2.0 * rand(rng) - 1.0)
        y = d.center[2] + d.radius * (2.0 * rand(rng) - 1.0)
        z = d.center[3] + d.radius * (2.0 * rand(rng) - 1.0)
        p = (x, y, z)
        in_domain(d, p) && return p
    end
end

# --- BoxDomain ---

struct BoxDomain <: AbstractDomain
    min_corner::NTuple{3, Float64}
    max_corner::NTuple{3, Float64}
end

domain_bounds(d::BoxDomain) = (d.min_corner, d.max_corner)

function in_domain(d::BoxDomain, point)
    return (d.min_corner[1] <= point[1] <= d.max_corner[1] &&
            d.min_corner[2] <= point[2] <= d.max_corner[2] &&
            d.min_corner[3] <= point[3] <= d.max_corner[3])
end

function signed_distance(d::BoxDomain, point)
    # Signed distance to axis-aligned box
    dx = max(d.min_corner[1] - point[1], point[1] - d.max_corner[1], 0.0)
    dy = max(d.min_corner[2] - point[2], point[2] - d.max_corner[2], 0.0)
    dz = max(d.min_corner[3] - point[3], point[3] - d.max_corner[3], 0.0)

    outside_dist = sqrt(dx * dx + dy * dy + dz * dz)

    if outside_dist > 0.0
        return outside_dist
    end

    # Inside: negative distance to nearest face
    inner_x = min(point[1] - d.min_corner[1], d.max_corner[1] - point[1])
    inner_y = min(point[2] - d.min_corner[2], d.max_corner[2] - point[2])
    inner_z = min(point[3] - d.min_corner[3], d.max_corner[3] - point[3])
    return -min(inner_x, inner_y, inner_z)
end

function sample_point(d::BoxDomain, rng::AbstractRNG)
    x = d.min_corner[1] + rand(rng) * (d.max_corner[1] - d.min_corner[1])
    y = d.min_corner[2] + rand(rng) * (d.max_corner[2] - d.min_corner[2])
    z = d.min_corner[3] + rand(rng) * (d.max_corner[3] - d.min_corner[3])
    return (x, y, z)
end

# --- EllipsoidDomain ---

struct EllipsoidDomain <: AbstractDomain
    center::NTuple{3, Float64}
    semi_axes::NTuple{3, Float64}   # (a, b, c) semi-axis lengths
end

domain_bounds(d::EllipsoidDomain) = (
    (d.center[1] - d.semi_axes[1], d.center[2] - d.semi_axes[2], d.center[3] - d.semi_axes[3]),
    (d.center[1] + d.semi_axes[1], d.center[2] + d.semi_axes[2], d.center[3] + d.semi_axes[3]),
)

function in_domain(d::EllipsoidDomain, point)
    nx = (point[1] - d.center[1]) / d.semi_axes[1]
    ny = (point[2] - d.center[2]) / d.semi_axes[2]
    nz = (point[3] - d.center[3]) / d.semi_axes[3]
    return nx * nx + ny * ny + nz * nz <= 1.0
end

function signed_distance(d::EllipsoidDomain, point)
    # Approximate SDF for ellipsoid (exact SDF requires iterative solve)
    # Uses the normalized distance - 1, scaled by the minimum semi-axis
    nx = (point[1] - d.center[1]) / d.semi_axes[1]
    ny = (point[2] - d.center[2]) / d.semi_axes[2]
    nz = (point[3] - d.center[3]) / d.semi_axes[3]
    normalized_dist = sqrt(nx * nx + ny * ny + nz * nz)
    min_axis = min(d.semi_axes[1], d.semi_axes[2], d.semi_axes[3])
    return (normalized_dist - 1.0) * min_axis
end

function sample_point(d::EllipsoidDomain, rng::AbstractRNG)
    while true
        x = d.center[1] + d.semi_axes[1] * (2.0 * rand(rng) - 1.0)
        y = d.center[2] + d.semi_axes[2] * (2.0 * rand(rng) - 1.0)
        z = d.center[3] + d.semi_axes[3] * (2.0 * rand(rng) - 1.0)
        p = (x, y, z)
        in_domain(d, p) && return p
    end
end

# --- EllipsoidShellDomain ---

"""
    EllipsoidShellDomain

Thin shell between two concentric ellipsoids. Models the myocardial wall
where coronary arteries course. The outer boundary is the epicardial surface;
the inner boundary is at `(1 - thickness)` fraction of the semi-axes.

Coronary arteries grow within this shell, naturally wrapping around the
heart-shaped surface rather than filling the interior volume.
"""
struct EllipsoidShellDomain <: AbstractDomain
    center::NTuple{3, Float64}
    semi_axes::NTuple{3, Float64}   # outer ellipsoid semi-axis lengths
    thickness::Float64               # shell thickness as fraction (0, 1)
end

domain_bounds(d::EllipsoidShellDomain) = (
    (d.center[1] - d.semi_axes[1], d.center[2] - d.semi_axes[2], d.center[3] - d.semi_axes[3]),
    (d.center[1] + d.semi_axes[1], d.center[2] + d.semi_axes[2], d.center[3] + d.semi_axes[3]),
)

function in_domain(d::EllipsoidShellDomain, point)
    nx = (point[1] - d.center[1]) / d.semi_axes[1]
    ny = (point[2] - d.center[2]) / d.semi_axes[2]
    nz = (point[3] - d.center[3]) / d.semi_axes[3]
    r2 = nx * nx + ny * ny + nz * nz
    inner = 1.0 - d.thickness
    return r2 <= 1.0 && r2 >= inner * inner
end

function signed_distance(d::EllipsoidShellDomain, point)
    nx = (point[1] - d.center[1]) / d.semi_axes[1]
    ny = (point[2] - d.center[2]) / d.semi_axes[2]
    nz = (point[3] - d.center[3]) / d.semi_axes[3]
    r = sqrt(nx * nx + ny * ny + nz * nz)
    inner = 1.0 - d.thickness
    min_axis = min(d.semi_axes[1], d.semi_axes[2], d.semi_axes[3])

    if r > 1.0
        return (r - 1.0) * min_axis            # outside outer
    elseif r < inner
        return (inner - r) * min_axis           # inside inner hole
    else
        # Inside shell: negative distance to nearest boundary
        d_to_outer = (1.0 - r) * min_axis
        d_to_inner = (r - inner) * min_axis
        return -min(d_to_outer, d_to_inner)
    end
end

function sample_point(d::EllipsoidShellDomain, rng::AbstractRNG)
    # Parametric sampling: 100% acceptance, volumetrically uniform in the shell.
    # The shell consists of nested similar ellipsoids at radial fraction f ∈ [1-t, 1].
    # Volume element ∝ f² df sinθ dθ dφ, so sample f³ uniformly for uniform volume.
    a, b, c = d.semi_axes
    t = d.thickness
    cosθ = 2.0 * rand(rng) - 1.0
    sinθ = sqrt(max(0.0, 1.0 - cosθ * cosθ))
    φ = 2.0 * Float64(π) * rand(rng)
    f_min3 = (1.0 - t)^3
    f = cbrt(f_min3 + (1.0 - f_min3) * rand(rng))
    return (
        d.center[1] + a * f * sinθ * cos(φ),
        d.center[2] + b * f * sinθ * sin(φ),
        d.center[3] + c * f * cosθ,
    )
end

# --- CSVVolumeDomain ---

"""
    CSVVolumeDomain

Volume domain reconstructed from a surface point cloud and outward normals.
The surface comes from Wenbo's CSV export (`heart_points_unique.csv`,
`heart_normals_unique.csv`). Interior sampling is precomputed once so the
existing CCO code can keep calling `sample_point(domain, rng)` unchanged.
"""
struct CSVVolumeDomain <: AbstractDomain
    surface_points::Matrix{Float64}
    surface_normals::Matrix{Float64}
    interior_points::Matrix{Float64}
    center::NTuple{3, Float64}
    min_corner::NTuple{3, Float64}
    max_corner::NTuple{3, Float64}
    volume::Float64
    characteristic_length::Float64
    length_scale::Float64
    query_grid::PointCloudGrid
end

domain_bounds(d::CSVVolumeDomain) = (d.min_corner, d.max_corner)

function _read_xyz_csv(path::AbstractString)
    raw = readdlm(path, ',', Float64)
    mat = raw isa Matrix ? Matrix{Float64}(raw) : reshape(Float64.(raw), 1, :)
    size(mat, 2) == 3 || error("Expected Nx3 CSV at $path, got size $(size(mat)).")
    return mat
end

function _row_tuple(mat::Matrix{Float64}, i::Int)
    return (mat[i, 1], mat[i, 2], mat[i, 3])
end

function _normalize_rows!(mat::Matrix{Float64})
    for i in axes(mat, 1)
        nx = mat[i, 1]
        ny = mat[i, 2]
        nz = mat[i, 3]
        nrm = sqrt(nx * nx + ny * ny + nz * nz)
        nrm > 1e-12 || error("Normal row $i has near-zero norm.")
        mat[i, 1] = nx / nrm
        mat[i, 2] = ny / nrm
        mat[i, 3] = nz / nrm
    end
    return mat
end

function _point_grid_index(grid::PointCloudGrid, x::Float64, y::Float64, z::Float64)
    ix = clamp(floor(Int, (x - grid.origin[1]) / grid.cell_size) + 1, 1, grid.dims[1])
    iy = clamp(floor(Int, (y - grid.origin[2]) / grid.cell_size) + 1, 1, grid.dims[2])
    iz = clamp(floor(Int, (z - grid.origin[3]) / grid.cell_size) + 1, 1, grid.dims[3])
    return ix + (iy - 1) * grid.dims[1] + (iz - 1) * grid.dims[1] * grid.dims[2]
end

function _point_grid_dims(lo::NTuple{3, Float64}, hi::NTuple{3, Float64}, cell_size::Float64)
    return (
        max(1, ceil(Int, (hi[1] - lo[1]) / cell_size)),
        max(1, ceil(Int, (hi[2] - lo[2]) / cell_size)),
        max(1, ceil(Int, (hi[3] - lo[3]) / cell_size)),
    )
end

function _default_point_grid_size(points::Matrix{Float64}, lo::NTuple{3, Float64}, hi::NTuple{3, Float64})
    dx = hi[1] - lo[1]
    dy = hi[2] - lo[2]
    dz = hi[3] - lo[3]
    extent = max(dx, dy, dz)
    bbox_vol = max(dx * dy * dz, 1e-12)
    nominal_spacing = cbrt(bbox_vol / max(size(points, 1), 1))
    return max(nominal_spacing * 2.5, extent / 64.0, 1e-4)
end

function _build_point_grid(points::Matrix{Float64}, lo::NTuple{3, Float64}, hi::NTuple{3, Float64})
    cell_size = _default_point_grid_size(points, lo, hi)
    dims = _point_grid_dims(lo, hi, cell_size)
    cells = Dict{Int, Vector{Int}}()
    grid = PointCloudGrid(cell_size, lo, dims, cells)

    for i in axes(points, 1)
        idx = _point_grid_index(grid, points[i, 1], points[i, 2], points[i, 3])
        push!(get!(cells, idx, Int[]), i)
    end
    return grid
end

function _surface_candidates(grid::PointCloudGrid, point; max_rings::Int=2, min_candidates::Int=32)
    x, y, z = point
    cx = clamp(floor(Int, (x - grid.origin[1]) / grid.cell_size) + 1, 1, grid.dims[1])
    cy = clamp(floor(Int, (y - grid.origin[2]) / grid.cell_size) + 1, 1, grid.dims[2])
    cz = clamp(floor(Int, (z - grid.origin[3]) / grid.cell_size) + 1, 1, grid.dims[3])

    candidates = Int[]
    for ring in 0:max_rings
        for dz in -ring:ring
            iz = cz + dz
            (iz < 1 || iz > grid.dims[3]) && continue
            for dy in -ring:ring
                iy = cy + dy
                (iy < 1 || iy > grid.dims[2]) && continue
                for dx in -ring:ring
                    ix = cx + dx
                    (ix < 1 || ix > grid.dims[1]) && continue
                    idx = ix + (iy - 1) * grid.dims[1] + (iz - 1) * grid.dims[1] * grid.dims[2]
                    haskey(grid.cells, idx) || continue
                    append!(candidates, grid.cells[idx])
                end
            end
        end
        length(candidates) >= min_candidates && break
    end
    return candidates
end

function _nearest_surface_info(
    points::Matrix{Float64},
    normals::Matrix{Float64},
    grid::PointCloudGrid,
    point,
)
    candidates = _surface_candidates(grid, point)
    if isempty(candidates)
        candidates = collect(axes(points, 1))
    end

    px, py, pz = point
    best_idx = first(candidates)
    best_d2 = Inf
    for idx in candidates
        dx = px - points[idx, 1]
        dy = py - points[idx, 2]
        dz = pz - points[idx, 3]
        d2 = dx * dx + dy * dy + dz * dz
        if d2 < best_d2
            best_d2 = d2
            best_idx = idx
        end
    end

    dx = px - points[best_idx, 1]
    dy = py - points[best_idx, 2]
    dz = pz - points[best_idx, 3]
    nx = normals[best_idx, 1]
    ny = normals[best_idx, 2]
    nz = normals[best_idx, 3]
    signed_dist = dx * nx + dy * ny + dz * nz
    return best_idx, signed_dist
end

function _nearest_surface_info(domain::CSVVolumeDomain, point)
    return _nearest_surface_info(domain.surface_points, domain.surface_normals, domain.query_grid, point)
end

function csv_volume_domain(
    points_csv::AbstractString,
    normals_csv::AbstractString;
    rng::AbstractRNG=Random.default_rng(),
    target_interior::Int=120_000,
    max_candidates::Int=400_000,
    batch_size::Int=8192,
    inside_slack::Float64=0.0,
    coordinate_scale::Float64=10.0,
)
    surface_points = _read_xyz_csv(points_csv)
    surface_normals = _normalize_rows!(_read_xyz_csv(normals_csv))
    size(surface_points, 1) == size(surface_normals, 1) ||
        error("Point/normal CSV row count mismatch: $(size(surface_points, 1)) vs $(size(surface_normals, 1)).")

    # Wenbo's CSV geometry and seeds are in centimeter-like units, while
    # VesselTree internally uses millimeters. Scale geometry into mm so the
    # existing growth heuristics and Murray-law radii remain physically consistent.
    coordinate_scale > 0.0 || error("coordinate_scale must be positive.")
    surface_points .*= coordinate_scale

    n = size(surface_points, 1)
    cx = sum(@view surface_points[:, 1]) / n
    cy = sum(@view surface_points[:, 2]) / n
    cz = sum(@view surface_points[:, 3]) / n
    center = (cx, cy, cz)

    lo = (
        minimum(@view surface_points[:, 1]),
        minimum(@view surface_points[:, 2]),
        minimum(@view surface_points[:, 3]),
    )
    hi = (
        maximum(@view surface_points[:, 1]),
        maximum(@view surface_points[:, 2]),
        maximum(@view surface_points[:, 3]),
    )

    grid = _build_point_grid(surface_points, lo, hi)

    keep = Matrix{Float64}(undef, target_interior, 3)
    kept = 0
    draws = 0

    while kept < target_interior && draws < max_candidates
        n_batch = min(batch_size, max_candidates - draws)
        for _ in 1:n_batch
            px = lo[1] + rand(rng) * (hi[1] - lo[1])
            py = lo[2] + rand(rng) * (hi[2] - lo[2])
            pz = lo[3] + rand(rng) * (hi[3] - lo[3])
            draws += 1

            _, sd = _nearest_surface_info(surface_points, surface_normals, grid, (px, py, pz))
            if sd <= inside_slack
                kept += 1
                keep[kept, 1] = px
                keep[kept, 2] = py
                keep[kept, 3] = pz
                kept >= target_interior && break
            end
            draws >= max_candidates && break
        end
    end

    kept > 0 || error("Failed to sample interior points from $points_csv.")
    interior_points = copy(keep[1:kept, :])

    bbox_dx = hi[1] - lo[1]
    bbox_dy = hi[2] - lo[2]
    bbox_dz = hi[3] - lo[3]
    bbox_vol = bbox_dx * bbox_dy * bbox_dz
    volume = bbox_vol * kept / max(draws, 1)
    characteristic_length = cbrt(max(volume, 1e-12))

    return CSVVolumeDomain(
        surface_points,
        surface_normals,
        interior_points,
        center,
        lo,
        hi,
        volume,
        characteristic_length,
        coordinate_scale,
        grid,
    )
end

function default_coronary_volume_domain(; kwargs...)
    repo_root = dirname(dirname(@__DIR__))
    return csv_volume_domain(
        joinpath(repo_root, "heart_points_unique.csv"),
        joinpath(repo_root, "heart_normals_unique.csv");
        kwargs...,
    )
end

function in_domain(d::CSVVolumeDomain, point)
    _, signed_dist = _nearest_surface_info(d, point)
    return signed_dist <= 0.0
end

function signed_distance(d::CSVVolumeDomain, point)
    _, signed_dist = _nearest_surface_info(d, point)
    return signed_dist
end

function sample_point(d::CSVVolumeDomain, rng::AbstractRNG)
    idx = rand(rng, axes(d.interior_points, 1))
    return _row_tuple(d.interior_points, idx)
end

# --- project_to_domain ---

"""
    project_to_domain(domain, point) -> NTuple{3, Float64}

Project a point to the nearest location inside the domain.
Returns the point unchanged if already inside.
"""
function project_to_domain(d::SphereDomain, point)
    dx = point[1] - d.center[1]
    dy = point[2] - d.center[2]
    dz = point[3] - d.center[3]
    r = sqrt(dx * dx + dy * dy + dz * dz)
    r <= d.radius && return point
    scale = d.radius / r
    return (d.center[1] + dx * scale, d.center[2] + dy * scale, d.center[3] + dz * scale)
end

function project_to_domain(d::BoxDomain, point)
    x = clamp(point[1], d.min_corner[1], d.max_corner[1])
    y = clamp(point[2], d.min_corner[2], d.max_corner[2])
    z = clamp(point[3], d.min_corner[3], d.max_corner[3])
    return (x, y, z)
end

function project_to_domain(d::EllipsoidDomain, point)
    nx = (point[1] - d.center[1]) / d.semi_axes[1]
    ny = (point[2] - d.center[2]) / d.semi_axes[2]
    nz = (point[3] - d.center[3]) / d.semi_axes[3]
    r = sqrt(nx * nx + ny * ny + nz * nz)
    r <= 1.0 && return point
    scale = 1.0 / r
    return (
        d.center[1] + nx * scale * d.semi_axes[1],
        d.center[2] + ny * scale * d.semi_axes[2],
        d.center[3] + nz * scale * d.semi_axes[3],
    )
end

function project_to_domain(d::EllipsoidShellDomain, point)
    nx = (point[1] - d.center[1]) / d.semi_axes[1]
    ny = (point[2] - d.center[2]) / d.semi_axes[2]
    nz = (point[3] - d.center[3]) / d.semi_axes[3]
    r = sqrt(nx * nx + ny * ny + nz * nz)
    inner = 1.0 - d.thickness

    if r < 1e-15
        # At center — push to middle of shell along +x
        mid = (inner + 1.0) / 2.0
        return (d.center[1] + mid * d.semi_axes[1], d.center[2], d.center[3])
    end

    if r > 1.0
        scale = 1.0 / r           # project to outer surface
    elseif r < inner
        scale = inner / r          # project to inner surface
    else
        return point               # already inside shell
    end

    return (
        d.center[1] + nx * scale * d.semi_axes[1],
        d.center[2] + ny * scale * d.semi_axes[2],
        d.center[3] + nz * scale * d.semi_axes[3],
    )
end

function project_to_domain(d::CSVVolumeDomain, point)
    idx, signed_dist = _nearest_surface_info(d, point)
    signed_dist <= 0.0 && return point

    px, py, pz = point
    proj = (
        px - signed_dist * d.surface_normals[idx, 1],
        py - signed_dist * d.surface_normals[idx, 2],
        pz - signed_dist * d.surface_normals[idx, 3],
    )

    for _ in 1:3
        in_domain(d, proj) && return proj
        idx2, sd2 = _nearest_surface_info(d, proj)
        proj = (
            proj[1] - (sd2 + 1e-6) * d.surface_normals[idx2, 1],
            proj[2] - (sd2 + 1e-6) * d.surface_normals[idx2, 2],
            proj[3] - (sd2 + 1e-6) * d.surface_normals[idx2, 3],
        )
    end

    return proj
end

"""
    default_coronary_domain() -> EllipsoidShellDomain

Default heart-shaped shell domain for coronary tree generation.
Semi-axes approximate a human heart (~100mm x 70mm x 90mm outer).
Shell thickness 0.3 gives ~10-15mm myocardial wall.
"""
function default_coronary_domain()
    return EllipsoidShellDomain((0.0, 0.0, 0.0), (50.0, 35.0, 45.0), 0.3)
end
