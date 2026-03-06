# Perfusion domain types for vascular tree generation

"""
    AbstractDomain

Interface for perfusion domains. Subtypes must implement:
- `in_domain(domain, point) -> Bool`
- `sample_point(domain, rng) -> NTuple{3, Float64}`
- `signed_distance(domain, point) -> Float64` (negative inside, positive outside)
"""
abstract type AbstractDomain end

# --- SphereDomain ---

struct SphereDomain <: AbstractDomain
    center::NTuple{3, Float64}
    radius::Float64
end

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
