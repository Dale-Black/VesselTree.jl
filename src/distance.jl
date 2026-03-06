# AK-accelerated point-to-segment distance computation

"""
    point_segment_distance_sq(cx, cy, cz, px, py, pz, dx, dy, dz) -> Float64

Squared minimum distance from point (c) to line segment (p -> d).
Pure scalar function suitable for use inside AK kernels.
"""
@inline function point_segment_distance_sq(
    cx::Float64, cy::Float64, cz::Float64,
    px::Float64, py::Float64, pz::Float64,
    dx::Float64, dy::Float64, dz::Float64,
)
    # Vector from p to d
    ex = dx - px
    ey = dy - py
    ez = dz - pz

    # Vector from p to c
    fx = cx - px
    fy = cy - py
    fz = cz - pz

    seg_len_sq = ex * ex + ey * ey + ez * ez

    if seg_len_sq < 1e-30
        # Degenerate (zero-length) segment — distance to point p
        return fx * fx + fy * fy + fz * fz
    end

    # Project c onto line p->d, clamped to [0, 1]
    t = (fx * ex + fy * ey + fz * ez) / seg_len_sq
    t = max(0.0, min(1.0, t))

    # Closest point on segment
    qx = px + t * ex
    qy = py + t * ey
    qz = pz + t * ez

    # Squared distance
    rx = cx - qx
    ry = cy - qy
    rz = cz - qz
    return rx * rx + ry * ry + rz * rz
end

"""
    point_segment_distance(cx, cy, cz, px, py, pz, dx, dy, dz) -> Float64

Minimum distance from point (c) to line segment (p -> d).
"""
@inline function point_segment_distance(
    cx::Float64, cy::Float64, cz::Float64,
    px::Float64, py::Float64, pz::Float64,
    dx::Float64, dy::Float64, dz::Float64,
)
    return sqrt(point_segment_distance_sq(cx, cy, cz, px, py, pz, dx, dy, dz))
end

"""
    compute_all_distances!(distances, segments::SegmentData, cx, cy, cz, n)

Compute distance from candidate point (cx, cy, cz) to each of the first `n`
segments using AK.foreachindex. Results stored in `distances[1:n]`.
"""
function compute_all_distances!(
    distances::AbstractVector{Float64},
    segments::SegmentData,
    cx::Float64, cy::Float64, cz::Float64,
    n::Int,
)
    px = segments.proximal_x
    py = segments.proximal_y
    pz = segments.proximal_z
    dx = segments.distal_x
    dy = segments.distal_y
    dz = segments.distal_z

    dist_view = @view distances[1:n]
    AK.foreachindex(dist_view) do i
        @inbounds dist_view[i] = point_segment_distance(
            cx, cy, cz,
            px[i], py[i], pz[i],
            dx[i], dy[i], dz[i],
        )
    end
    return nothing
end

"""
    find_nearest_segments(distances, n, k) -> Vector{Int}

Return indices of the `k` nearest segments from `distances[1:n]`.
Uses AK.sortperm for parallel sorting.
"""
function find_nearest_segments(
    distances::AbstractVector{Float64},
    n::Int,
    k::Int,
)
    k = min(k, n)
    dist_view = @view distances[1:n]
    perm = AK.sortperm(dist_view)
    return perm[1:k]
end
