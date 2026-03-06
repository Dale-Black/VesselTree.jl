# AK-accelerated segment-segment intersection testing

"""
    segments_min_distance_sq(a1x,a1y,a1z, a2x,a2y,a2z, b1x,b1y,b1z, b2x,b2y,b2z) -> Float64

Squared minimum distance between two 3D line segments (a1→a2) and (b1→b2).
Uses the standard closest-point-on-two-segments algorithm.
"""
@inline function segments_min_distance_sq(
    a1x::Float64, a1y::Float64, a1z::Float64,
    a2x::Float64, a2y::Float64, a2z::Float64,
    b1x::Float64, b1y::Float64, b1z::Float64,
    b2x::Float64, b2y::Float64, b2z::Float64,
)
    # Direction vectors
    dx = a2x - a1x; dy = a2y - a1y; dz = a2z - a1z  # d = a2 - a1
    ex = b2x - b1x; ey = b2y - b1y; ez = b2z - b1z  # e = b2 - b1
    rx = a1x - b1x; ry = a1y - b1y; rz = a1z - b1z  # r = a1 - b1

    a = dx*dx + dy*dy + dz*dz  # |d|^2
    e_val = ex*ex + ey*ey + ez*ez  # |e|^2
    f = ex*rx + ey*ry + ez*rz  # e·r

    if a < 1e-30 && e_val < 1e-30
        # Both segments are points
        return rx*rx + ry*ry + rz*rz
    end

    if a < 1e-30
        # Segment a is a point
        t = 0.0
        s = clamp(f / e_val, 0.0, 1.0)
    else
        c = dx*rx + dy*ry + dz*rz  # d·r
        if e_val < 1e-30
            # Segment b is a point
            s = 0.0
            t = clamp(-c / a, 0.0, 1.0)
        else
            b_val = dx*ex + dy*ey + dz*ez  # d·e
            denom = a * e_val - b_val * b_val
            if abs(denom) < 1e-30
                # Parallel segments
                t = 0.0
            else
                t = clamp((b_val * f - c * e_val) / denom, 0.0, 1.0)
            end
            s = (b_val * t + f) / e_val
            if s < 0.0
                s = 0.0
                t = clamp(-c / a, 0.0, 1.0)
            elseif s > 1.0
                s = 1.0
                t = clamp((b_val - c) / a, 0.0, 1.0)
            end
        end
    end

    # Closest points
    px = a1x + t * dx - (b1x + s * ex)
    py = a1y + t * dy - (b1y + s * ey)
    pz = a1z + t * dz - (b1z + s * ez)

    return px*px + py*py + pz*pz
end

"""
    segments_intersect(a1x,a1y,a1z, a2x,a2y,a2z, b1x,b1y,b1z, b2x,b2y,b2z, min_dist) -> Bool

Returns true if the minimum distance between segments (a1→a2) and (b1→b2) is less than `min_dist`.
"""
@inline function segments_intersect(
    a1x::Float64, a1y::Float64, a1z::Float64,
    a2x::Float64, a2y::Float64, a2z::Float64,
    b1x::Float64, b1y::Float64, b1z::Float64,
    b2x::Float64, b2y::Float64, b2z::Float64,
    min_dist::Float64,
)
    return segments_min_distance_sq(a1x,a1y,a1z, a2x,a2y,a2z, b1x,b1y,b1z, b2x,b2y,b2z) < min_dist * min_dist
end

"""
    check_intersections!(results, segments::SegmentData, npx,npy,npz, ndx,ndy,ndz, min_dist, n)

AK kernel: check if the new segment (np→nd) intersects any of the first `n` existing segments.
`results[i]` is set to `true` if segment `i` is within `min_dist` of the new segment.
"""
function check_intersections!(
    results::AbstractVector{Bool},
    segments::SegmentData,
    npx::Float64, npy::Float64, npz::Float64,
    ndx::Float64, ndy::Float64, ndz::Float64,
    min_dist::Float64,
    n::Int,
)
    px = segments.proximal_x
    py = segments.proximal_y
    pz = segments.proximal_z
    dx = segments.distal_x
    dy = segments.distal_y
    dz = segments.distal_z

    res_view = @view results[1:n]
    AK.foreachindex(res_view) do i
        @inbounds res_view[i] = segments_intersect(
            px[i], py[i], pz[i],
            dx[i], dy[i], dz[i],
            npx, npy, npz,
            ndx, ndy, ndz,
            min_dist,
        )
    end
    return nothing
end

"""
    has_any_intersection(results, n) -> Bool

Check if any segment intersects, using AK.any for parallel reduction.
"""
function has_any_intersection(results::AbstractVector{Bool}, n::Int)
    res_view = @view results[1:n]
    return AK.any(identity, res_view)
end

"""
    check_domain_crossing(domain::AbstractDomain, p1, p2; n_samples=10) -> Bool

Check if a segment from p1 to p2 leaves the domain by sampling points along the segment.
Returns `true` if the segment crosses the domain boundary (i.e., some sampled point is outside).
"""
function check_domain_crossing(domain::AbstractDomain, p1, p2; n_samples::Int=10)
    for i in 1:n_samples
        t = i / (n_samples + 1)
        x = p1[1] + t * (p2[1] - p1[1])
        y = p1[2] + t * (p2[2] - p1[2])
        z = p1[3] + t * (p2[3] - p1[3])
        if !in_domain(domain, (x, y, z))
            return true
        end
    end
    return false
end
