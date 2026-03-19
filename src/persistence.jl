using Dates

struct SparseContrastVolume{T<:AbstractFloat}
    origin_mm::NTuple{3, Float64}
    spacing_mm::NTuple{3, Float64}
    dims::NTuple{3, Int}
    times_s::Vector{Float64}
    voxel_indices_ijk::Matrix{Int32}      # n_active x 3
    occupancy_fraction::Vector{T}
    concentration_mg_mL::Matrix{T}        # n_active x nt
end

function Base.show(io::IO, volume::SparseContrastVolume)
    print(
        io,
        "SparseContrastVolume(",
        "dims=", volume.dims,
        ", active_voxels=", size(volume.voxel_indices_ijk, 1),
        ", nt=", length(volume.times_s),
        ")",
    )
end

struct _SegmentVoxelContribution{T<:AbstractFloat}
    linear_index::Int32
    tree_index::Int32
    segment_id::Int32
    weight::T
end

function format_run_timestamp(ts::DateTime=Dates.now())
    return Dates.format(ts, dateformat"yyyymmddTHHMMSS")
end

function _normalize_timestamp(ts::Union{Nothing, DateTime, AbstractString})
    if ts === nothing
        return format_run_timestamp()
    elseif ts isa DateTime
        return format_run_timestamp(ts)
    else
        return String(ts)
    end
end

function _normalize_spacing_mm(spacing_mm::Real)
    s = Float64(spacing_mm)
    s > 0.0 || error("Expected positive spacing, got $spacing_mm")
    return (s, s, s)
end

function _normalize_spacing_mm(spacing_mm::NTuple{3, <:Real})
    s = (Float64(spacing_mm[1]), Float64(spacing_mm[2]), Float64(spacing_mm[3]))
    all(x -> x > 0.0, s) || error("Expected positive spacing, got $spacing_mm")
    return s
end

function _normalize_padding_mm(padding_mm::Real)
    p = Float64(padding_mm)
    p >= 0.0 || error("Expected non-negative padding, got $padding_mm")
    return (p, p, p)
end

function _normalize_padding_mm(padding_mm::NTuple{3, <:Real})
    p = (Float64(padding_mm[1]), Float64(padding_mm[2]), Float64(padding_mm[3]))
    all(x -> x >= 0.0, p) || error("Expected non-negative padding, got $padding_mm")
    return p
end

function _tree_segment_centers(tree::VascularTree)
    seg = tree.segments
    n = seg.n
    centers = Matrix{Float64}(undef, n, 3)
    for i in 1:n
        centers[i, 1] = 0.5 * (seg.proximal_x[i] + seg.distal_x[i])
        centers[i, 2] = 0.5 * (seg.proximal_y[i] + seg.distal_y[i])
        centers[i, 3] = 0.5 * (seg.proximal_z[i] + seg.distal_z[i])
    end
    return centers
end

function _forest_bounds_mm(forest::CoronaryForest)
    lo = (Inf, Inf, Inf)
    hi = (-Inf, -Inf, -Inf)
    for tree in values(forest.trees)
        seg = tree.segments
        for i in 1:seg.n
            r = seg.radius[i]
            lo = (
                min(lo[1], seg.proximal_x[i] - r, seg.distal_x[i] - r),
                min(lo[2], seg.proximal_y[i] - r, seg.distal_y[i] - r),
                min(lo[3], seg.proximal_z[i] - r, seg.distal_z[i] - r),
            )
            hi = (
                max(hi[1], seg.proximal_x[i] + r, seg.distal_x[i] + r),
                max(hi[2], seg.proximal_y[i] + r, seg.distal_y[i] + r),
                max(hi[3], seg.proximal_z[i] + r, seg.distal_z[i] + r),
            )
        end
    end
    return lo, hi
end

function _segment_voxel_hits(
    px::Float64, py::Float64, pz::Float64,
    dx::Float64, dy::Float64, dz::Float64,
    radius_mm::Float64,
    origin_mm::NTuple{3, Float64},
    spacing_mm::NTuple{3, Float64},
    dims::NTuple{3, Int};
    supersample::Int=2,
)
    supersample >= 1 || error("supersample must be >= 1")
    radius_sq = radius_mm^2
    min_x = min(px, dx) - radius_mm
    max_x = max(px, dx) + radius_mm
    min_y = min(py, dy) - radius_mm
    max_y = max(py, dy) + radius_mm
    min_z = min(pz, dz) - radius_mm
    max_z = max(pz, dz) + radius_mm

    ix_lo = clamp(floor(Int, (min_x - origin_mm[1]) / spacing_mm[1]) + 1, 1, dims[1])
    ix_hi = clamp(ceil(Int, (max_x - origin_mm[1]) / spacing_mm[1]), 1, dims[1])
    iy_lo = clamp(floor(Int, (min_y - origin_mm[2]) / spacing_mm[2]) + 1, 1, dims[2])
    iy_hi = clamp(ceil(Int, (max_y - origin_mm[2]) / spacing_mm[2]), 1, dims[2])
    iz_lo = clamp(floor(Int, (min_z - origin_mm[3]) / spacing_mm[3]) + 1, 1, dims[3])
    iz_hi = clamp(ceil(Int, (max_z - origin_mm[3]) / spacing_mm[3]), 1, dims[3])

    lin = LinearIndices(dims)
    hits = Tuple{Int32, Float32}[]
    total_subsamples = supersample^3

    for iz in iz_lo:iz_hi, iy in iy_lo:iy_hi, ix in ix_lo:ix_hi
        inside = 0
        for sz in 1:supersample, sy in 1:supersample, sx in 1:supersample
            cx = origin_mm[1] + (ix - 1 + (sx - 0.5) / supersample) * spacing_mm[1]
            cy = origin_mm[2] + (iy - 1 + (sy - 0.5) / supersample) * spacing_mm[2]
            cz = origin_mm[3] + (iz - 1 + (sz - 0.5) / supersample) * spacing_mm[3]
            if point_segment_distance_sq(cx, cy, cz, px, py, pz, dx, dy, dz) <= radius_sq
                inside += 1
            end
        end

        if inside > 0
            weight = Float32(inside / total_subsamples)
            push!(hits, (Int32(lin[ix, iy, iz]), weight))
        end
    end

    return hits
end

function sparse_contrast_to_dense(
    volume::SparseContrastVolume;
    storage_type::Type{T}=eltype(volume.concentration_mg_mL),
    fill_value::Real=0.0,
) where {T<:AbstractFloat}
    nx, ny, nz = volume.dims
    nt = length(volume.times_s)
    dense = fill(T(fill_value), nx, ny, nz, nt)
    for row in 1:size(volume.voxel_indices_ijk, 1)
        i = Int(volume.voxel_indices_ijk[row, 1])
        j = Int(volume.voxel_indices_ijk[row, 2])
        k = Int(volume.voxel_indices_ijk[row, 3])
        dense[i, j, k, :] .= T.(volume.concentration_mg_mL[row, :])
    end
    return dense
end

function rasterize_forest_contrast_sparse(
    forest::CoronaryForest,
    results::Dict{String, <:ContrastTransportResult};
    spacing_mm::Union{Real, NTuple{3, <:Real}}=(1.0, 1.0, 1.0),
    padding_mm::Union{Real, NTuple{3, <:Real}}=(2.0, 2.0, 2.0),
    supersample::Int=2,
    min_radius_um::Float64=0.0,
    storage_type::Type{T}=Float32,
) where {T<:AbstractFloat}
    names = sort!(collect(keys(forest.trees)))
    isempty(names) && error("Forest has no trees to rasterize.")

    spacing = _normalize_spacing_mm(spacing_mm)
    padding = _normalize_padding_mm(padding_mm)
    times = Float64.(results[first(names)].times)
    nt = length(times)
    for name in names[2:end]
        results[name].times == times || error("All contrast results must share the same time axis.")
    end

    lo, hi = _forest_bounds_mm(forest)
    origin = (lo[1] - padding[1], lo[2] - padding[2], lo[3] - padding[3])
    extent = (
        (hi[1] - lo[1]) + 2.0 * padding[1],
        (hi[2] - lo[2]) + 2.0 * padding[2],
        (hi[3] - lo[3]) + 2.0 * padding[3],
    )
    dims = (
        max(1, ceil(Int, extent[1] / spacing[1])),
        max(1, ceil(Int, extent[2] / spacing[2])),
        max(1, ceil(Int, extent[3] / spacing[3])),
    )

    occupancy_by_linear = Dict{Int32, Float64}()
    contributions = _SegmentVoxelContribution{Float32}[]

    for (tree_idx, name) in enumerate(names)
        tree = forest.trees[name]
        seg = tree.segments
        for seg_id in 1:seg.n
            radius_um = seg.radius[seg_id] * 2000.0
            radius_um < min_radius_um && continue
            hits = _segment_voxel_hits(
                seg.proximal_x[seg_id], seg.proximal_y[seg_id], seg.proximal_z[seg_id],
                seg.distal_x[seg_id], seg.distal_y[seg_id], seg.distal_z[seg_id],
                seg.radius[seg_id],
                origin,
                spacing,
                dims;
                supersample=supersample,
            )
            for (linear_index, weight) in hits
                occupancy_by_linear[linear_index] = get(occupancy_by_linear, linear_index, 0.0) + weight
                push!(contributions, _SegmentVoxelContribution{Float32}(linear_index, Int32(tree_idx), Int32(seg_id), weight))
            end
        end
    end

    active_linear = sort!(collect(keys(occupancy_by_linear)))
    n_active = length(active_linear)
    row_of_linear = Dict{Int32, Int32}(active_linear[i] => Int32(i) for i in eachindex(active_linear))
    occupancy = Vector{T}(undef, n_active)
    concentration = zeros(T, n_active, nt)
    ijk = Matrix{Int32}(undef, n_active, 3)
    cart = CartesianIndices(dims)

    for (row, linear_index) in enumerate(active_linear)
        occupancy[row] = T(occupancy_by_linear[linear_index])
        idx = Tuple(cart[Int(linear_index)])
        ijk[row, 1] = Int32(idx[1])
        ijk[row, 2] = Int32(idx[2])
        ijk[row, 3] = Int32(idx[3])
    end

    for contrib in contributions
        row = Int(row_of_linear[contrib.linear_index])
        tree_name = names[Int(contrib.tree_index)]
        seg_id = Int(contrib.segment_id)
        concentration[row, :] .+= T(contrib.weight) .* T.(results[tree_name].concentration[seg_id, :])
    end

    for row in 1:n_active
        occupancy[row] > 0 || continue
        concentration[row, :] ./= occupancy[row]
    end

    return SparseContrastVolume{T}(origin, spacing, dims, times, ijk, occupancy, concentration)
end

function save_tree_snapshot(
    directory::AbstractString,
    tree::VascularTree;
    timestamp::Union{Nothing, DateTime, AbstractString}=nothing,
    formats::AbstractVector{Symbol}=[:jld2, :csv, :graph_json, :wenbo_txt],
)
    ts = _normalize_timestamp(timestamp)
    isdir(directory) || mkpath(directory)
    stem = joinpath(directory, "$(tree.name)-$(ts)")
    paths = Dict{String, String}()

    for fmt in formats
        if fmt == :jld2
            path = stem * ".jld2"
            save_tree(path, tree)
            paths["jld2"] = path
        elseif fmt == :csv
            path = stem * ".csv"
            export_csv(tree, path)
            paths["csv"] = path
        elseif fmt == :graph_json
            path = stem * ".json"
            export_graph_json(tree, path)
            paths["graph_json"] = path
        elseif fmt == :wenbo_txt
            path = stem * ".txt"
            export_wenbo_txt(tree, path)
            paths["wenbo_txt"] = path
        elseif fmt == :vtp
            path = export_centerlines_vtp(tree, stem * "-centerlines")
            paths["vtp"] = path
        elseif fmt == :stl
            path = stem * ".stl"
            export_stl(tree, path)
            paths["stl"] = path
        else
            error("Unsupported tree snapshot format: $fmt")
        end
    end

    return paths
end

function save_forest_snapshot(
    directory::AbstractString,
    forest::CoronaryForest;
    timestamp::Union{Nothing, DateTime, AbstractString}=nothing,
    formats::AbstractVector{Symbol}=[:jld2, :csv, :graph_json, :wenbo_txt],
)
    ts = _normalize_timestamp(timestamp)
    isdir(directory) || mkpath(directory)
    tree_paths = Dict{String, Dict{String, String}}()
    for name in sort!(collect(keys(forest.trees)))
        tree_paths[name] = save_tree_snapshot(directory, forest.trees[name]; timestamp=ts, formats=formats)
    end

    manifest_path = joinpath(directory, "forest-$(ts).json")
    manifest = Dict(
        "timestamp" => ts,
        "trees" => [
            Dict(
                "name" => name,
                "n_segments" => forest.trees[name].segments.n,
                "files" => tree_paths[name],
            ) for name in sort!(collect(keys(tree_paths)))
        ],
    )
    open(manifest_path, "w") do io
        _write_json(io, manifest)
    end

    return (timestamp=ts, tree_paths=tree_paths, manifest_path=manifest_path)
end

function _segment_snapshot_payload(tree::VascularTree, result::ContrastTransportResult)
    seg = tree.segments
    topo = tree.topology
    return Dict(
        "times_s" => result.times,
        "concentration_mg_mL" => result.concentration,
        "transit_time_s" => result.transit_time_s,
        "segment_volume_mL" => result.segment_volume_mL,
        "segment_centers_mm" => _tree_segment_centers(tree),
        "segment_radius_mm" => copy(seg.radius[1:seg.n]),
        "segment_length_mm" => copy(seg.seg_length[1:seg.n]),
        "segment_flow_m3_s" => copy(seg.flow[1:seg.n]),
        "strahler_order" => copy(topo.strahler_order[1:seg.n]),
    )
end

function save_contrast_snapshot(
    directory::AbstractString,
    forest::CoronaryForest,
    results::Dict{String, <:ContrastTransportResult};
    timestamp::Union{Nothing, DateTime, AbstractString}=nothing,
    voxel_spacing_mm::Union{Real, NTuple{3, <:Real}}=(1.0, 1.0, 1.0),
    voxel_padding_mm::Union{Real, NTuple{3, <:Real}}=(2.0, 2.0, 2.0),
    voxel_supersample::Int=2,
    min_radius_um::Float64=0.0,
    storage_type::Type{T}=Float32,
    dense_mode::Symbol=:auto,
    dense_max_bytes::Int=1_500_000_000,
) where {T<:AbstractFloat}
    dense_mode in (:never, :auto, :always) || error("dense_mode must be :never, :auto, or :always")
    ts = _normalize_timestamp(timestamp)
    isdir(directory) || mkpath(directory)

    segment_path = joinpath(directory, "contrast-segments-$(ts).jld2")
    tree_payload = Dict{String, Any}()
    for name in sort!(collect(keys(results)))
        tree_payload[name] = _segment_snapshot_payload(forest.trees[name], results[name])
    end
    JLD2.jldsave(segment_path; tree_results=tree_payload)

    sparse_volume = rasterize_forest_contrast_sparse(
        forest,
        results;
        spacing_mm=voxel_spacing_mm,
        padding_mm=voxel_padding_mm,
        supersample=voxel_supersample,
        min_radius_um=min_radius_um,
        storage_type=storage_type,
    )

    dense_bytes = prod(sparse_volume.dims) * length(sparse_volume.times_s) * sizeof(storage_type)
    save_dense = dense_mode == :always || (dense_mode == :auto && dense_bytes <= dense_max_bytes)
    volume_path = joinpath(directory, "contrast-volume-$(ts).jld2")

    if save_dense
        concentration_4d = sparse_contrast_to_dense(sparse_volume; storage_type=storage_type)
        JLD2.jldsave(
            volume_path;
            origin_mm=collect(sparse_volume.origin_mm),
            spacing_mm=collect(sparse_volume.spacing_mm),
            dims=collect(sparse_volume.dims),
            times_s=sparse_volume.times_s,
            voxel_indices_ijk=sparse_volume.voxel_indices_ijk,
            occupancy_fraction=sparse_volume.occupancy_fraction,
            concentration_sparse_mg_mL=sparse_volume.concentration_mg_mL,
            concentration_4d_mg_mL=concentration_4d,
        )
    else
        JLD2.jldsave(
            volume_path;
            origin_mm=collect(sparse_volume.origin_mm),
            spacing_mm=collect(sparse_volume.spacing_mm),
            dims=collect(sparse_volume.dims),
            times_s=sparse_volume.times_s,
            voxel_indices_ijk=sparse_volume.voxel_indices_ijk,
            occupancy_fraction=sparse_volume.occupancy_fraction,
            concentration_sparse_mg_mL=sparse_volume.concentration_mg_mL,
        )
    end

    manifest_path = joinpath(directory, "contrast-$(ts).json")
    manifest = Dict(
        "timestamp" => ts,
        "segment_snapshot" => segment_path,
        "volume_snapshot" => volume_path,
        "dense_4d_saved" => save_dense,
        "spacing_mm" => collect(sparse_volume.spacing_mm),
        "dims" => collect(sparse_volume.dims),
        "times" => Dict(
            "n" => length(sparse_volume.times_s),
            "t0_s" => first(sparse_volume.times_s),
            "tend_s" => last(sparse_volume.times_s),
        ),
        "active_voxels" => size(sparse_volume.voxel_indices_ijk, 1),
    )
    open(manifest_path, "w") do io
        _write_json(io, manifest)
    end

    return (
        timestamp=ts,
        segment_path=segment_path,
        volume_path=volume_path,
        manifest_path=manifest_path,
        sparse_volume=sparse_volume,
        dense_4d_saved=save_dense,
    )
end

function _open_uniform_knots(n_ctrl::Int, degree::Int)
    n_ctrl > degree || error("Need n_ctrl > degree, got n_ctrl=$n_ctrl degree=$degree")
    knots = Float64[]
    append!(knots, zeros(Float64, degree + 1))
    n_internal = n_ctrl - degree - 1
    if n_internal > 0
        denom = n_ctrl - degree
        for i in 1:n_internal
            push!(knots, i / denom)
        end
    end
    append!(knots, ones(Float64, degree + 1))
    return knots
end

function _format_xcat_scalar(x::Real; digits::Int=9)
    y = round(Float64(x), digits=digits)
    if abs(y) < 10.0^(-digits)
        y = 0.0
    end
    s = string(y)
    return endswith(s, ".0") ? s[1:end-2] : s
end

function _tree_segment_surface(
    tree::VascularTree,
    seg_id::Int;
    name_prefix::AbstractString="grown",
    circumferential_points::Int=8,
)
    seg = tree.segments
    p0 = SVector(seg.proximal_x[seg_id], seg.proximal_y[seg_id], seg.proximal_z[seg_id])
    p1 = SVector(seg.distal_x[seg_id], seg.distal_y[seg_id], seg.distal_z[seg_id])
    r0 = seg.radius[seg_id]
    r1 = seg.radius[seg_id]

    axis = p1 - p0
    axis_norm = norm(axis)
    if axis_norm < 1e-12
        axis = SVector(1.0, 0.0, 0.0)
        axis_norm = 1.0
    end
    axis = axis / axis_norm
    u1, u2, u3 = _find_perpendicular(axis[1], axis[2], axis[3])
    u = SVector(u1, u2, u3)
    v = cross(axis, u)

    n_u = circumferential_points + 1
    control_points = Array{SVector{3, Float64}, 2}(undef, 2, n_u)
    for j in 1:n_u
        θ = 2.0 * π * (j - 1) / circumferential_points
        direction = cos(θ) * u + sin(θ) * v
        control_points[1, j] = p0 + r0 * direction
        control_points[2, j] = p1 + r1 * direction
    end

    return XCATNurbsSurface(
        string(name_prefix, "_", tree.name, "_", lpad(seg_id, 7, '0')),
        2,
        n_u,
        _open_uniform_knots(n_u, 1),
        [0.0, 0.0, 1.0, 1.0],
        control_points,
    )
end

function write_xcat_surface(io::IO, surface::XCATNurbsSurface)
    println(io, surface.name)
    println(io, "$(surface.m) :M")
    println(io, "$(surface.n) :N")
    println(io, "U Knot Vector")
    for knot in surface.u_knots
        println(io, _format_xcat_scalar(knot; digits=6))
    end
    println(io, "V Knot Vector")
    for knot in surface.v_knots
        println(io, _format_xcat_scalar(knot; digits=6))
    end
    println(io, "Control Points")
    for i in 1:surface.m, j in 1:surface.n
        p = surface.control_points[i, j]
        println(
            io,
            _format_xcat_scalar(p[1]), ",",
            _format_xcat_scalar(p[2]), ",",
            _format_xcat_scalar(p[3]),
        )
    end
    println(io)
    return io
end

function write_xcat_nrb(path::AbstractString, surfaces::AbstractVector{XCATNurbsSurface})
    open(path, "w") do io
        for surface in surfaces
            write_xcat_surface(io, surface)
        end
    end
    return path
end

function export_fused_nrb_model(
    original_nrb_path::AbstractString,
    output_path::AbstractString,
    forest::CoronaryForest;
    name_prefix::AbstractString="grown",
    circumferential_points::Int=8,
)
    cp(original_nrb_path, output_path; force=true)
    open(output_path, "a") do io
        println(io)
        for name in sort!(collect(keys(forest.trees)))
            tree = forest.trees[name]
            for seg_id in 1:tree.segments.n
                surface = _tree_segment_surface(
                    tree,
                    seg_id;
                    name_prefix=name_prefix,
                    circumferential_points=circumferential_points,
                )
                write_xcat_surface(io, surface)
            end
        end
    end
    return output_path
end

function save_xcat_run_artifacts(
    output_root::AbstractString,
    original_nrb_path::AbstractString,
    forest::CoronaryForest,
    contrast_results::Dict{String, <:ContrastTransportResult};
    timestamp::Union{Nothing, DateTime, AbstractString}=nothing,
    tree_formats::AbstractVector{Symbol}=[:jld2, :csv, :graph_json, :wenbo_txt],
    voxel_spacing_mm::Union{Real, NTuple{3, <:Real}}=(1.0, 1.0, 1.0),
    voxel_padding_mm::Union{Real, NTuple{3, <:Real}}=(2.0, 2.0, 2.0),
    voxel_supersample::Int=2,
    min_radius_um::Float64=0.0,
    dense_mode::Symbol=:auto,
    dense_max_bytes::Int=1_500_000_000,
    nrb_name_prefix::AbstractString="grown",
    nrb_circumferential_points::Int=8,
)
    ts = _normalize_timestamp(timestamp)
    run_dir = joinpath(output_root, ts)
    tree_dir = joinpath(run_dir, "trees")
    contrast_dir = joinpath(run_dir, "contrast")
    model_dir = joinpath(run_dir, "model")
    mkpath(tree_dir)
    mkpath(contrast_dir)
    mkpath(model_dir)

    tree_snapshot = save_forest_snapshot(tree_dir, forest; timestamp=ts, formats=tree_formats)
    contrast_snapshot = save_contrast_snapshot(
        contrast_dir,
        forest,
        contrast_results;
        timestamp=ts,
        voxel_spacing_mm=voxel_spacing_mm,
        voxel_padding_mm=voxel_padding_mm,
        voxel_supersample=voxel_supersample,
        min_radius_um=min_radius_um,
        dense_mode=dense_mode,
        dense_max_bytes=dense_max_bytes,
    )

    fused_nrb_path = joinpath(model_dir, "xcat-grown-$(ts).nrb")
    export_fused_nrb_model(
        original_nrb_path,
        fused_nrb_path,
        forest;
        name_prefix=nrb_name_prefix,
        circumferential_points=nrb_circumferential_points,
    )

    manifest_path = joinpath(run_dir, "run-$(ts).json")
    manifest = Dict(
        "timestamp" => ts,
        "original_nrb" => original_nrb_path,
        "trees_manifest" => tree_snapshot.manifest_path,
        "contrast_manifest" => contrast_snapshot.manifest_path,
        "fused_nrb" => fused_nrb_path,
    )
    open(manifest_path, "w") do io
        _write_json(io, manifest)
    end

    return (
        timestamp=ts,
        run_dir=run_dir,
        tree_snapshot=tree_snapshot,
        contrast_snapshot=contrast_snapshot,
        fused_nrb_path=fused_nrb_path,
        manifest_path=manifest_path,
    )
end
