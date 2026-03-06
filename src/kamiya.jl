# Kamiya bifurcation point optimization
# Finds optimal bifurcation point and radii when connecting a new terminal to an existing segment.

"""
    compute_radii(r_parent::Float64, gamma::Float64) -> (r_left, r_right)

Given parent radius and Murray's law gamma, compute daughter radii for a symmetric bifurcation.
r_parent^gamma = r_left^gamma + r_right^gamma
For symmetric: r_daughter = r_parent / 2^(1/gamma)
"""
@inline function compute_radii_symmetric(r_parent::Float64, gamma::Float64)
    r_daughter = r_parent / 2.0^(1.0 / gamma)
    return r_daughter, r_daughter
end

"""
    compute_radii_asymmetric(r_parent::Float64, asymmetry::Float64, gamma::Float64) -> (r_large, r_small)

Given parent radius, asymmetry ratio (r_small/r_large), and gamma:
r_parent^gamma = r_large^gamma + r_small^gamma
r_small = asymmetry * r_large
=> r_large = r_parent / (1 + asymmetry^gamma)^(1/gamma)
"""
@inline function compute_radii_asymmetric(r_parent::Float64, asymmetry::Float64, gamma::Float64)
    r_large = r_parent / (1.0 + asymmetry^gamma)^(1.0 / gamma)
    r_small = asymmetry * r_large
    return r_large, r_small
end

"""
    surface_cost_at_t(t, seg_px, seg_py, seg_pz, seg_dx, seg_dy, seg_dz,
                      seg_radius, tx, ty, tz, gamma) -> Float64

Compute surface cost for bifurcation at parametric position t in [0,1] along
the existing segment (seg_p → seg_d), connecting to new terminal point (tx,ty,tz).

The bifurcation point splits the segment into:
- Upper segment: parent → bifurcation (length l0, radius r0 = seg_radius)
- Lower segment: bifurcation → original distal (length l1, radius r1)
- New segment: bifurcation → terminal (length l2, radius r2)

Radii r1, r2 are determined by Murray's law: r0^gamma = r1^gamma + r2^gamma

Surface cost = 2π(r0*l0 + r1*l1 + r2*l2)
"""
@inline function surface_cost_at_t(
    t::Float64,
    seg_px::Float64, seg_py::Float64, seg_pz::Float64,
    seg_dx::Float64, seg_dy::Float64, seg_dz::Float64,
    seg_radius::Float64,
    tx::Float64, ty::Float64, tz::Float64,
    gamma::Float64,
)
    # Bifurcation point
    bx = seg_px + t * (seg_dx - seg_px)
    by = seg_py + t * (seg_dy - seg_py)
    bz = seg_pz + t * (seg_dz - seg_pz)

    # Lengths
    l0_x = bx - seg_px; l0_y = by - seg_py; l0_z = bz - seg_pz
    l0 = sqrt(l0_x*l0_x + l0_y*l0_y + l0_z*l0_z)

    l1_x = seg_dx - bx; l1_y = seg_dy - by; l1_z = seg_dz - bz
    l1 = sqrt(l1_x*l1_x + l1_y*l1_y + l1_z*l1_z)

    l2_x = tx - bx; l2_y = ty - by; l2_z = tz - bz
    l2 = sqrt(l2_x*l2_x + l2_y*l2_y + l2_z*l2_z)

    # Parent radius (upstream from bifurcation to original proximal)
    r0 = seg_radius

    # Daughter radii via Murray's law (symmetric split for now)
    r1, r2 = compute_radii_symmetric(r0, gamma)

    # Surface cost (proportional to 2*pi*r*l, drop 2*pi constant)
    return r0 * l0 + r1 * l1 + r2 * l2
end

"""
    optimize_bifurcation_point(seg_px, seg_py, seg_pz, seg_dx, seg_dy, seg_dz,
                               seg_radius, tx, ty, tz, gamma;
                               n_iter=20, t_min=0.01, t_max=0.99) -> (t_opt, cost)

Golden section search along segment to minimize surface cost.
Returns optimal parametric position and cost.
"""
function optimize_bifurcation_point(
    seg_px::Float64, seg_py::Float64, seg_pz::Float64,
    seg_dx::Float64, seg_dy::Float64, seg_dz::Float64,
    seg_radius::Float64,
    tx::Float64, ty::Float64, tz::Float64,
    gamma::Float64;
    n_iter::Int=20, t_min::Float64=0.01, t_max::Float64=0.99,
)
    phi = (sqrt(5.0) - 1.0) / 2.0  # golden ratio conjugate ≈ 0.618

    a = t_min
    b = t_max

    for _ in 1:n_iter
        c = b - phi * (b - a)
        d = a + phi * (b - a)

        fc = surface_cost_at_t(c, seg_px, seg_py, seg_pz, seg_dx, seg_dy, seg_dz,
                               seg_radius, tx, ty, tz, gamma)
        fd = surface_cost_at_t(d, seg_px, seg_py, seg_pz, seg_dx, seg_dy, seg_dz,
                               seg_radius, tx, ty, tz, gamma)

        if fc < fd
            b = d
        else
            a = c
        end
    end

    t_opt = (a + b) / 2.0
    cost = surface_cost_at_t(t_opt, seg_px, seg_py, seg_pz, seg_dx, seg_dy, seg_dz,
                             seg_radius, tx, ty, tz, gamma)
    return t_opt, cost
end

"""
    evaluate_all_connections!(costs, bifurc_t, segments, n, tx, ty, tz, gamma)

AK kernel: for each of the first `n` segments, run golden section search to find
optimal bifurcation point and store cost + position.
"""
function evaluate_all_connections!(
    costs::AbstractVector{Float64},
    bifurc_t::AbstractVector{Float64},
    segments::SegmentData,
    n::Int,
    tx::Float64, ty::Float64, tz::Float64,
    gamma::Float64,
)
    px = segments.proximal_x
    py = segments.proximal_y
    pz = segments.proximal_z
    dx = segments.distal_x
    dy = segments.distal_y
    dz = segments.distal_z
    radii = segments.radius

    costs_view = @view costs[1:n]
    t_view = @view bifurc_t[1:n]

    AK.foreachindex(costs_view) do i
        @inbounds t_opt, cost = optimize_bifurcation_point(
            px[i], py[i], pz[i],
            dx[i], dy[i], dz[i],
            radii[i],
            tx, ty, tz,
            gamma,
        )
        @inbounds costs_view[i] = cost
        @inbounds t_view[i] = t_opt
    end
    return nothing
end

"""
    select_best_connection(costs, bifurc_t, n) -> (seg_idx, t_opt, cost)

Find the connection with minimum cost. Uses AK.minimum for the parallel min,
then a linear scan for the index (argmin not available in AK).
"""
function select_best_connection(
    costs::AbstractVector{Float64},
    bifurc_t::AbstractVector{Float64},
    n::Int,
)
    costs_view = @view costs[1:n]
    min_cost = AK.minimum(costs_view)

    # Linear scan for index of minimum (AK has no argmin)
    min_idx = 1
    for i in 1:n
        if costs[i] == min_cost
            min_idx = i
            break
        end
    end

    return min_idx, bifurc_t[min_idx], costs[min_idx]
end
