using LinearAlgebra
using StaticArrays
using DelimitedFiles
using Random
using Statistics

struct XCATNurbsSurface
    name::String
    m::Int
    n::Int
    u_knots::Vector{Float64}
    v_knots::Vector{Float64}
    control_points::Array{SVector{3, Float64}, 2}
end

function Base.show(io::IO, s::XCATNurbsSurface)
    print(io, "XCATNurbsSurface($(repr(s.name)), size=$(size(s.control_points)))")
end

function xcat_bounds(surface::XCATNurbsSurface)
    xs = map(p -> p[1], surface.control_points)
    ys = map(p -> p[2], surface.control_points)
    zs = map(p -> p[3], surface.control_points)
    return (minimum(xs), minimum(ys), minimum(zs)), (maximum(xs), maximum(ys), maximum(zs))
end

function xcat_center(surface::XCATNurbsSurface)
    acc = SVector(0.0, 0.0, 0.0)
    for p in surface.control_points
        acc += p
    end
    return acc / length(surface.control_points)
end

function xcat_object_dict(surfaces::AbstractVector{XCATNurbsSurface})
    return Dict(surface.name => surface for surface in surfaces)
end

function xcat_select_objects(
    surfaces::AbstractVector{XCATNurbsSurface};
    include_regex::Union{Nothing, Regex}=nothing,
    names::Union{Nothing, AbstractVector{<:AbstractString}}=nothing,
)
    selected = surfaces
    if include_regex !== nothing
        selected = [surface for surface in selected if occursin(include_regex, surface.name)]
    end
    if names !== nothing
        wanted = Set(String(name) for name in names)
        selected = [surface for surface in selected if surface.name in wanted]
    end
    return selected
end

function xcat_summary_rows(surfaces::AbstractVector{XCATNurbsSurface})
    rows = NamedTuple[]
    for surface in surfaces
        lo, hi = xcat_bounds(surface)
        nu, nv = xcat_uv_counts(surface)
        pu, pv = xcat_degrees(surface)
        push!(rows, (
            name=surface.name,
            nu=nu,
            nv=nv,
            pu=pu,
            pv=pv,
            min_x=lo[1], min_y=lo[2], min_z=lo[3],
            max_x=hi[1], max_y=hi[2], max_z=hi[3],
        ))
    end
    return rows
end

function xcat_uv_counts(surface::XCATNurbsSurface)
    # XCAT stores control points in M x N order, where:
    # - V-direction count = M = number of rows
    # - U-direction count = N = number of columns
    n_v, n_u = size(surface.control_points)
    return n_u, n_v
end

function xcat_degrees(surface::XCATNurbsSurface)
    n_u, n_v = xcat_uv_counts(surface)
    p_u = length(surface.u_knots) - n_u - 1
    p_v = length(surface.v_knots) - n_v - 1
    p_u >= 0 || error("Invalid inferred U degree for $(surface.name)")
    p_v >= 0 || error("Invalid inferred V degree for $(surface.name)")
    return p_u, p_v
end

function _xcat_basis(i::Int, p::Int, u::Float64, knots::Vector{Float64}, n_basis::Int)
    if p == 0
        left = knots[i]
        right = knots[i + 1]
        return ((left <= u < right) || (u == knots[end] && i == n_basis && left <= u <= right)) ? 1.0 : 0.0
    end

    left_den = knots[i + p] - knots[i]
    right_den = knots[i + p + 1] - knots[i + 1]

    left_term = 0.0
    right_term = 0.0

    if left_den > 0
        left_term = ((u - knots[i]) / left_den) * _xcat_basis(i, p - 1, u, knots, n_basis)
    end
    if right_den > 0
        right_term = ((knots[i + p + 1] - u) / right_den) * _xcat_basis(i + 1, p - 1, u, knots, n_basis)
    end

    return left_term + right_term
end

function xcat_surface_point(surface::XCATNurbsSurface, u::Real, v::Real)
    u_f = clamp(Float64(u), 0.0, 1.0)
    v_f = clamp(Float64(v), 0.0, 1.0)
    n_u, n_v = xcat_uv_counts(surface)
    p_u, p_v = xcat_degrees(surface)

    basis_u = [_xcat_basis(i, p_u, u_f, surface.u_knots, n_u) for i in 1:n_u]
    basis_v = [_xcat_basis(j, p_v, v_f, surface.v_knots, n_v) for j in 1:n_v]

    acc = SVector(0.0, 0.0, 0.0)
    @inbounds for j in 1:n_v
        b_v = basis_v[j]
        b_v == 0.0 && continue
        for i in 1:n_u
            b = b_v * basis_u[i]
            b == 0.0 && continue
            acc += b * surface.control_points[j, i]
        end
    end
    return acc
end

function xcat_surface_normal(surface::XCATNurbsSurface, u::Real, v::Real; δ::Float64=1e-4)
    u0 = max(0.0, Float64(u) - δ)
    u1 = min(1.0, Float64(u) + δ)
    v0 = max(0.0, Float64(v) - δ)
    v1 = min(1.0, Float64(v) + δ)

    pu0 = xcat_surface_point(surface, u0, v)
    pu1 = xcat_surface_point(surface, u1, v)
    pv0 = xcat_surface_point(surface, u, v0)
    pv1 = xcat_surface_point(surface, u, v1)

    tu = pu1 - pu0
    tv = pv1 - pv0
    n = cross(tu, tv)
    n_norm = norm(n)
    n_norm > 0 || return SVector(0.0, 0.0, 1.0)
    return n / n_norm
end

function xcat_sample_surface(
    surface::XCATNurbsSurface;
    n_u::Int=80,
    n_v::Int=80,
    orient_outward::Bool=false,
)
    u_values = collect(range(0.0, 1.0; length=n_u))
    v_values = collect(range(0.0, 1.0; length=n_v))
    points = Array{SVector{3, Float64}, 2}(undef, n_v, n_u)
    normals = Array{SVector{3, Float64}, 2}(undef, n_v, n_u)

    for j in 1:n_v
        for i in 1:n_u
            u = u_values[i]
            v = v_values[j]
            points[j, i] = xcat_surface_point(surface, u, v)
            normals[j, i] = xcat_surface_normal(surface, u, v)
        end
    end

    if orient_outward
        center = let acc = SVector(0.0, 0.0, 0.0)
            for p in points
                acc += p
            end
            acc / length(points)
        end
        for j in 1:n_v, i in 1:n_u
            p = points[j, i]
            n = normals[j, i]
            if dot(n, p - center) < 0
                normals[j, i] = -n
            end
        end
    end

    return points, normals, u_values, v_values
end

function xcat_export_sampled_surface_csv(
    surface::XCATNurbsSurface,
    points_path::AbstractString,
    normals_path::AbstractString;
    n_u::Int=120,
    n_v::Int=120,
    orient_outward::Bool=false,
)
    points, normals, _, _ = xcat_sample_surface(surface; n_u=n_u, n_v=n_v, orient_outward=orient_outward)
    point_rows = Matrix{Float64}(undef, length(points), 3)
    normal_rows = Matrix{Float64}(undef, length(normals), 3)

    k = 1
    for j in 1:size(points, 1), i in 1:size(points, 2)
        p = points[j, i]
        n = normals[j, i]
        point_rows[k, 1] = p[1]
        point_rows[k, 2] = p[2]
        point_rows[k, 3] = p[3]
        normal_rows[k, 1] = n[1]
        normal_rows[k, 2] = n[2]
        normal_rows[k, 3] = n[3]
        k += 1
    end

    writedlm(points_path, point_rows, ',')
    writedlm(normals_path, normal_rows, ',')
    return points_path, normals_path
end

struct XCATCenterline
    name::String
    centers::Vector{SVector{3, Float64}}
    radii::Vector{Float64}
    axial_param::Symbol
end

function Base.show(io::IO, c::XCATCenterline)
    print(io, "XCATCenterline($(repr(c.name)), n=$(length(c.centers)), axial=$(c.axial_param))")
end

struct XCATTreeConnection
    parent_segment::String
    child_segment::String
    parent_index::Int
    child_index::Int
    gap_mm::Float64
end

struct XCATCenterlineTree
    name::String
    segments::Dict{String, XCATCenterline}
    root_segment::String
    connections::Vector{XCATTreeConnection}
end

function Base.show(io::IO, tree::XCATCenterlineTree)
    print(io, "XCATCenterlineTree($(repr(tree.name)), segments=$(length(tree.segments)), connections=$(length(tree.connections)))")
end

function xcat_reverse_centerline(centerline::XCATCenterline)
    return XCATCenterline(
        centerline.name,
        reverse(centerline.centers),
        reverse(centerline.radii),
        centerline.axial_param,
    )
end

function _xcat_slice_centerline_from(centerline::XCATCenterline, start_idx::Int)
    start_idx = clamp(start_idx, 1, length(centerline.centers))
    return XCATCenterline(
        centerline.name,
        centerline.centers[start_idx:end],
        centerline.radii[start_idx:end],
        centerline.axial_param,
    )
end

function xcat_sampled_surface_rows(
    surface::XCATNurbsSurface;
    n_u::Int=120,
    n_v::Int=120,
    orient_outward::Bool=false,
)
    points, normals, _, _ = xcat_sample_surface(surface; n_u=n_u, n_v=n_v, orient_outward=orient_outward)
    point_rows = Matrix{Float64}(undef, length(points), 3)
    normal_rows = Matrix{Float64}(undef, length(normals), 3)
    k = 1
    for j in 1:size(points, 1), i in 1:size(points, 2)
        p = points[j, i]
        n = normals[j, i]
        point_rows[k, 1] = p[1]
        point_rows[k, 2] = p[2]
        point_rows[k, 3] = p[3]
        normal_rows[k, 1] = n[1]
        normal_rows[k, 2] = n[2]
        normal_rows[k, 3] = n[3]
        k += 1
    end
    return point_rows, normal_rows
end

function xcat_surface_axis(surface::XCATNurbsSurface)
    pts = surface.control_points
    u_closure = mean(norm.(pts[:, 1] .- pts[:, end]))
    v_closure = mean(norm.(pts[1, :] .- pts[end, :]))

    # The closed/circumferential direction has the smaller first-vs-last gap.
    # The axial direction is the other one.
    if u_closure <= v_closure
        return :v
    else
        return :u
    end
end

function _trim_centerline_caps(
    centers::Vector{SVector{3, Float64}},
    radii::Vector{Float64};
    min_radius_fraction::Float64=0.2,
)
    positive = filter(>(0.0), radii)
    isempty(positive) && return centers, radii
    threshold = min_radius_fraction * median(positive)

    first_keep = findfirst(r -> r >= threshold, radii)
    last_keep = findlast(r -> r >= threshold, radii)
    (first_keep === nothing || last_keep === nothing || first_keep > last_keep) && return centers, radii

    return centers[first_keep:last_keep], radii[first_keep:last_keep]
end

function xcat_centerline_from_surface(
    surface::XCATNurbsSurface;
    circumferential_samples::Int=48,
    axial_samples::Union{Nothing, Int}=nothing,
)
    n_u, n_v = xcat_uv_counts(surface)
    axial = xcat_surface_axis(surface)

    if axial === :v
        n_axial = something(axial_samples, max(16, n_v * 2))
        n_circ = max(12, circumferential_samples)
        points, _, _, _ = xcat_sample_surface(surface; n_u=n_circ, n_v=n_axial, orient_outward=false)
        centers = Vector{SVector{3, Float64}}(undef, n_axial)
        radii = Vector{Float64}(undef, n_axial)
        for j in 1:n_axial
            row = points[j, :]
            center = let acc = SVector(0.0, 0.0, 0.0)
                for p in row
                    acc += p
                end
                acc / length(row)
            end
            centers[j] = center
            radii[j] = mean(norm.(row .- Ref(center)))
        end
        centers, radii = _trim_centerline_caps(centers, radii)
        return XCATCenterline(surface.name, centers, radii, :v)
    else
        n_axial = something(axial_samples, max(16, n_u * 2))
        n_circ = max(12, circumferential_samples)
        points, _, _, _ = xcat_sample_surface(surface; n_u=n_axial, n_v=n_circ, orient_outward=false)
        centers = Vector{SVector{3, Float64}}(undef, n_axial)
        radii = Vector{Float64}(undef, n_axial)
        for i in 1:n_axial
            col = points[:, i]
            center = let acc = SVector(0.0, 0.0, 0.0)
                for p in col
                    acc += p
                end
                acc / length(col)
            end
            centers[i] = center
            radii[i] = mean(norm.(col .- Ref(center)))
        end
        centers, radii = _trim_centerline_caps(centers, radii)
        return XCATCenterline(surface.name, centers, radii, :u)
    end
end

function xcat_centerline_summary_rows(centerlines::AbstractVector{XCATCenterline})
    rows = NamedTuple[]
    for line in centerlines
        radii = line.radii
        centers = line.centers
        push!(rows, (
            name=line.name,
            n=length(centers),
            axial=String(line.axial_param),
            radius_min=minimum(radii),
            radius_max=maximum(radii),
            radius_mean=mean(radii),
            start_x=first(centers)[1],
            start_y=first(centers)[2],
            start_z=first(centers)[3],
            end_x=last(centers)[1],
            end_y=last(centers)[2],
            end_z=last(centers)[3],
        ))
    end
    return rows
end

function xcat_export_centerline_csv(
    centerline::XCATCenterline,
    path::AbstractString,
)
    rows = Matrix{Float64}(undef, length(centerline.centers), 4)
    for i in eachindex(centerline.centers)
        p = centerline.centers[i]
        rows[i, 1] = p[1]
        rows[i, 2] = p[2]
        rows[i, 3] = p[3]
        rows[i, 4] = centerline.radii[i]
    end
    writedlm(path, rows, ',')
    return path
end

function _point_to_centerline_distance(point::SVector{3, Float64}, centerline::XCATCenterline)
    return minimum(norm(point - c) for c in centerline.centers)
end

function _connection_distance(a::XCATCenterline, b::XCATCenterline)
    return norm(last(a.centers) - first(b.centers))
end

function _orient_root_to_aorta(centerline::XCATCenterline, aorta::XCATCenterline)
    d_start = _point_to_centerline_distance(first(centerline.centers), aorta)
    d_end = _point_to_centerline_distance(last(centerline.centers), aorta)
    return d_start <= d_end ? centerline : xcat_reverse_centerline(centerline)
end

function _orient_to_previous(previous::XCATCenterline, current::XCATCenterline)
    forward = _connection_distance(previous, current)
    reversed = _connection_distance(previous, xcat_reverse_centerline(current))
    return forward <= reversed ? current : xcat_reverse_centerline(current)
end

function _nearest_centerline_pair(parent::XCATCenterline, child::XCATCenterline)
    best_dist = Inf
    best_parent_idx = 1
    best_child_idx = 1
    for (i, p) in enumerate(parent.centers)
        for (j, q) in enumerate(child.centers)
            d = norm(p - q)
            if d < best_dist
                best_dist = d
                best_parent_idx = i
                best_child_idx = j
            end
        end
    end
    return best_parent_idx, best_child_idx, best_dist
end

function _orient_child_to_parent(parent::XCATCenterline, child::XCATCenterline)
    _, _, forward_dist = _nearest_centerline_pair(parent, child)
    reversed = xcat_reverse_centerline(child)
    _, _, reversed_dist = _nearest_centerline_pair(parent, reversed)
    return forward_dist <= reversed_dist ? child : reversed
end

function xcat_merge_centerline_chain(name::AbstractString, chain::Vector{XCATCenterline}; merge_tol_mm::Float64=2.0)
    isempty(chain) && error("Cannot merge empty centerline chain for $name")
    centers = copy(first(chain).centers)
    radii = copy(first(chain).radii)

    for segment in chain[2:end]
        seg_centers = segment.centers
        seg_radii = segment.radii
        if norm(last(centers) - first(seg_centers)) <= merge_tol_mm
            append!(centers, seg_centers[2:end])
            append!(radii, seg_radii[2:end])
        else
            append!(centers, seg_centers)
            append!(radii, seg_radii)
        end
    end

    return XCATCenterline(String(name), centers, radii, :v)
end

function xcat_build_coronary_trunks(centerlines::AbstractVector{XCATCenterline})
    cmap = Dict(line.name => line for line in centerlines)
    haskey(cmap, "dias_aorta") || error("`dias_aorta` centerline is required to orient coronary trunks.")
    aorta = cmap["dias_aorta"]

    lad_chain = [
        _orient_root_to_aorta(cmap["dias_lad1"], aorta),
    ]
    push!(lad_chain, _orient_to_previous(lad_chain[end], cmap["dias_lad2"]))
    push!(lad_chain, _orient_to_previous(lad_chain[end], cmap["dias_lad3"]))

    rca_chain = [
        _orient_root_to_aorta(cmap["dias_rca1"], aorta),
    ]
    push!(rca_chain, _orient_to_previous(rca_chain[end], cmap["dias_rca2"]))

    lcx_chain = [
        _orient_root_to_aorta(cmap["dias_lcx"], aorta),
    ]

    return Dict(
        "AORTA" => aorta,
        "LAD" => xcat_merge_centerline_chain("LAD", lad_chain),
        "LCX" => xcat_merge_centerline_chain("LCX", lcx_chain),
        "RCA" => xcat_merge_centerline_chain("RCA", rca_chain),
    )
end

function _make_tree_connection(parent::XCATCenterline, child::XCATCenterline)
    oriented_child = _orient_child_to_parent(parent, child)
    parent_idx, child_idx, gap = _nearest_centerline_pair(parent, oriented_child)
    trimmed_child = _xcat_slice_centerline_from(oriented_child, child_idx)
    connection = XCATTreeConnection(parent.name, trimmed_child.name, parent_idx, 1, gap)
    return trimmed_child, connection
end

function xcat_build_coronary_trees(centerlines::AbstractVector{XCATCenterline})
    cmap = Dict(line.name => line for line in centerlines)
    haskey(cmap, "dias_aorta") || error("`dias_aorta` centerline is required to orient coronary trees.")
    aorta = cmap["dias_aorta"]

    lad_root = _orient_root_to_aorta(cmap["dias_lad1"], aorta)
    lad_mid = _orient_to_previous(lad_root, cmap["dias_lad2"])
    lad_distal = _orient_to_previous(lad_mid, cmap["dias_lad3"])
    lad_segments = Dict(lad_root.name => lad_root)
    lad_mid_oriented, lad_conn1 = _make_tree_connection(lad_root, lad_mid)
    lad_segments[lad_mid_oriented.name] = lad_mid_oriented
    lad_distal_oriented, lad_conn2 = _make_tree_connection(lad_mid_oriented, lad_distal)
    lad_segments[lad_distal_oriented.name] = lad_distal_oriented
    lad_connections = [lad_conn1, lad_conn2]

    lcx_root = _orient_root_to_aorta(cmap["dias_lcx"], aorta)
    lcx_segments = Dict(lcx_root.name => lcx_root)
    lcx_connections = XCATTreeConnection[]

    rca_root = _orient_root_to_aorta(cmap["dias_rca1"], aorta)
    rca_child = _orient_to_previous(rca_root, cmap["dias_rca2"])
    rca_segments = Dict(rca_root.name => rca_root)
    rca_child_oriented, rca_conn = _make_tree_connection(rca_root, rca_child)
    rca_segments[rca_child_oriented.name] = rca_child_oriented
    rca_connections = [rca_conn]

    return Dict(
        "AORTA" => XCATCenterlineTree("AORTA", Dict(aorta.name => aorta), aorta.name, XCATTreeConnection[]),
        "LAD" => XCATCenterlineTree("LAD", lad_segments, lad_root.name, lad_connections),
        "LCX" => XCATCenterlineTree("LCX", lcx_segments, lcx_root.name, lcx_connections),
        "RCA" => XCATCenterlineTree("RCA", rca_segments, rca_root.name, rca_connections),
    )
end

function xcat_tree_summary_rows(trees::AbstractVector{XCATCenterlineTree})
    rows = NamedTuple[]
    for tree in trees
        n_points = sum(length(seg.centers) for seg in values(tree.segments))
        mean_radius = mean(vcat([seg.radii for seg in values(tree.segments)]...))
        max_gap = isempty(tree.connections) ? 0.0 : maximum(conn.gap_mm for conn in tree.connections)
        push!(rows, (
            name=tree.name,
            n_segments=length(tree.segments),
            n_points=n_points,
            mean_radius=mean_radius,
            max_gap=max_gap,
        ))
    end
    return rows
end

function xcat_default_cavity_names(; phase::AbstractString="dias")
    prefix = String(phase)
    return [
        "$(prefix)_lv_0",
        "$(prefix)_lv_1",
        "$(prefix)_lv_2",
        "$(prefix)_lv_3",
        "$(prefix)_lv_4",
        "$(prefix)_la",
        "$(prefix)_ra",
        "$(prefix)_rv",
    ]
end

function xcat_myocardial_shell_domain(
    nrb_path::AbstractString;
    phase::AbstractString="dias",
    outer_name::AbstractString="$(phase)_pericardium",
    cavity_names::Vector{String}=xcat_default_cavity_names(; phase=phase),
    outer_n_u::Int=160,
    outer_n_v::Int=96,
    cavity_n_u::Int=120,
    cavity_n_v::Int=72,
    rng::AbstractRNG=Random.default_rng(),
    target_interior::Int=120_000,
    max_candidates::Int=400_000,
    batch_size::Int=8192,
)
    surfaces = parse_xcat_nrb(nrb_path)
    object_map = xcat_object_dict(surfaces)

    haskey(object_map, outer_name) || error("Outer surface `$outer_name` not found in $nrb_path")
    outer_surface = object_map[outer_name]
    outer_points, outer_normals = xcat_sampled_surface_rows(
        outer_surface;
        n_u=outer_n_u,
        n_v=outer_n_v,
        orient_outward=true,
    )

    cavity_points = Matrix{Float64}[]
    cavity_normals = Matrix{Float64}[]
    actual_names = String[]
    for name in cavity_names
        haskey(object_map, name) || error("Cavity surface `$name` not found in $nrb_path")
        pts, nrms = xcat_sampled_surface_rows(
            object_map[name];
            n_u=cavity_n_u,
            n_v=cavity_n_v,
            orient_outward=true,
        )
        push!(cavity_points, pts)
        push!(cavity_normals, nrms)
        push!(actual_names, name)
    end

    domain = _shell_domain_from_matrices(
        outer_points,
        outer_normals,
        cavity_points,
        cavity_normals;
        cavity_names=actual_names,
        rng=rng,
        target_interior=target_interior,
        max_candidates=max_candidates,
        batch_size=batch_size,
        coordinate_scale=1.0,
    )

    return domain, surfaces
end

function _parse_xcat_control_point(line::AbstractString)
    parts = split(strip(line), ',')
    length(parts) == 3 || error("Expected x,y,z control point line, got: $line")
    return SVector(parse(Float64, parts[1]), parse(Float64, parts[2]), parse(Float64, parts[3]))
end

function _parse_xcat_knot_vector(lines::Vector{String}, idx::Int)
    knots = Float64[]
    while idx <= length(lines)
        token = strip(lines[idx])
        (token == "Control Points" || token == "U Knot Vector" || token == "V Knot Vector") && break
        isempty(token) && break
        push!(knots, parse(Float64, token))
        idx += 1
    end
    return knots, idx
end

function _reshape_xcat_control_points(points::Vector{SVector{3, Float64}}, m::Int, n::Int)
    nu = m
    nv = n
    expected = nu * nv
    length(points) == expected || error("Expected $expected control points, found $(length(points))")

    grid = Array{SVector{3, Float64}, 2}(undef, nu, nv)
    k = 1
    for u in 1:nu
        for v in 1:nv
            grid[u, v] = points[k]
            k += 1
        end
    end
    return grid
end

"""
    parse_xcat_nrb(path) -> Vector{XCATNurbsSurface}

Parse an XCAT `.nrb` ASCII file into named NURBS-like surface objects.

This parser preserves the control net exactly as stored in the file. Downstream
helpers can infer spline degrees from the knot vectors and sample the surfaces
for inspection or conversion.
"""
function parse_xcat_nrb(path::AbstractString)
    lines = readlines(path)
    surfaces = XCATNurbsSurface[]
    i = 1
    while i <= length(lines)
        name = strip(lines[i])
        if isempty(name)
            i += 1
            continue
        end
        if i + 2 > length(lines) ||
           !occursin(":M", strip(lines[i + 1])) ||
           !occursin(":N", strip(lines[i + 2]))
            i += 1
            continue
        end

        m = parse(Int, strip(split(strip(lines[i + 1]), ':')[1]))
        n = parse(Int, strip(split(strip(lines[i + 2]), ':')[1]))
        i += 3

        strip(lines[i]) == "U Knot Vector" || error("Expected `U Knot Vector` after $name")
        i += 1
        u_knots, i = _parse_xcat_knot_vector(lines, i)

        strip(lines[i]) == "Control Points" || begin
            strip(lines[i]) == "V Knot Vector" || error("Expected `V Knot Vector` after U knots for $name")
        end

        strip(lines[i]) == "V Knot Vector" || error("Expected `V Knot Vector` after U knots for $name")
        i += 1
        v_knots, i = _parse_xcat_knot_vector(lines, i)

        strip(lines[i]) == "Control Points" || error("Expected `Control Points` for $name")
        i += 1

        n_points = m * n
        points = Vector{SVector{3, Float64}}(undef, n_points)
        for k in 1:n_points
            i <= length(lines) || error("Unexpected EOF while reading control points for $name")
            points[k] = _parse_xcat_control_point(lines[i])
            i += 1
        end

        control_points = _reshape_xcat_control_points(points, m, n)
        push!(surfaces, XCATNurbsSurface(name, m, n, u_knots, v_knots, control_points))
    end
    return surfaces
end
