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

"""
    default_coronary_domain() -> EllipsoidShellDomain

Default heart-shaped shell domain for coronary tree generation.
Semi-axes approximate a human heart (~100mm x 70mm x 90mm outer).
Shell thickness 0.3 gives ~10-15mm myocardial wall.
"""
function default_coronary_domain()
    return EllipsoidShellDomain((0.0, 0.0, 0.0), (50.0, 35.0, 45.0), 0.3)
end
