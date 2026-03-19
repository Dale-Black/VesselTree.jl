using Mmap
using Unitful
using XLSX
import XrayAttenuation as XA

function _create_material_mixture(
    materials::Vector{<:XA.Material},
    fractions::Vector{Float64};
    by_volume::Bool=false,
    name::String="Custom Mixture",
)
    length(materials) == length(fractions) || error("materials and fractions must have same length")

    mass_fractions =
        if by_volume
            densities = [ustrip(u"g/cm^3", m.density) for m in materials]
            mass = fractions .* densities
            mass ./ sum(mass)
        else
            fractions ./ sum(fractions)
        end

    composition = Dict{Int, Float64}()
    for (material, fraction) in zip(materials, mass_fractions)
        for (Z, weight) in material.composition
            composition[Z] = get(composition, Z, 0.0) + weight * fraction
        end
    end

    density = sum(ustrip(u"g/cm^3", material.density) * fraction for (material, fraction) in zip(materials, mass_fractions))
    ZA = sum(material.ZA_ratio * fraction for (material, fraction) in zip(materials, mass_fractions))
    I = sum(ustrip(u"eV", material.I) * fraction for (material, fraction) in zip(materials, mass_fractions))
    return XA.Material(name, ZA, I * u"eV", density * u"g/cm^3", composition)
end

function _iodine_contrast_material(
    base_material::XA.Material,
    concentration_mg_per_mL::Float64;
    name::String="Iodine Contrast",
)
    iodine_g_per_cm3 = concentration_mg_per_mL / 1000.0
    base_density = ustrip(u"g/cm^3", base_material.density)
    total_density = base_density + iodine_g_per_cm3
    f_iodine = iodine_g_per_cm3 / total_density
    f_base = 1.0 - f_iodine

    composition = Dict{Int, Float64}()
    for (Z, weight) in base_material.composition
        composition[Z] = weight * f_base
    end
    composition[53] = get(composition, 53, 0.0) + f_iodine

    iodine = XA.Materials.iodine
    ZA = base_material.ZA_ratio * f_base + iodine.ZA_ratio * f_iodine
    I = ustrip(u"eV", base_material.I) * f_base + ustrip(u"eV", iodine.I) * f_iodine
    return XA.Material(name, ZA, I * u"eV", total_density * u"g/cm^3", composition)
end

struct XCATRawLabels{T<:Unsigned, A<:AbstractArray{T,3}}
    mask::A
    voxel_size_mm::NTuple{3, Float64}
    source_path::String
    label_log_path::Union{Nothing, String}
end

struct XCATRawAlignment
    origin_mm::NTuple{3, Float64}
    voxel_size_mm::NTuple{3, Float64}
    coronary_center_mm::NTuple{3, Float64}
    pericardium_center_mm::NTuple{3, Float64}
    raw_coronary_center_idx::NTuple{3, Float64}
    raw_pericardium_center_idx::NTuple{3, Float64}
end

function Base.show(io::IO, x::XCATRawAlignment)
    print(io, "XCATRawAlignment(origin_mm=$(x.origin_mm), voxel_size_mm=$(x.voxel_size_mm))")
end

function load_xcat_raw_labels(
    raw_path::AbstractString;
    dims::NTuple{3, Int}=(1600, 1400, 500),
    dtype::Type{T}=UInt8,
    reverse_dims::Tuple{Vararg{Int}}=(2, 3),
    voxel_size_mm::NTuple{3, <:Real}=(0.2, 0.2, 0.2),
    log_path::Union{Nothing, AbstractString}=nothing,
) where {T<:Unsigned}
    expected = prod(dims) * sizeof(dtype)
    actual = filesize(raw_path)
    actual == expected || error("File size mismatch for $raw_path: expected $expected bytes, found $actual")

    io = open(raw_path, "r")
    data = Mmap.mmap(io, Vector{dtype}, prod(dims))
    close(io)
    arr = reshape(data, dims)
    mask = copy(arr)
    for d in reverse_dims
        mask = reverse(mask; dims=d)
    end

    return XCATRawLabels(mask, Float64.(voxel_size_mm), String(raw_path), isnothing(log_path) ? nothing : String(log_path))
end

function _compute_ZA_ratio(composition::Dict{Int, Float64})
    atomic_masses = Dict(
        1=>1.008, 6=>12.011, 7=>14.007, 8=>15.999, 11=>22.990, 12=>24.305,
        15=>30.974, 16=>32.06, 17=>35.45, 19=>39.098, 20=>40.078, 26=>55.845, 53=>126.904,
    )
    zsum = 0.0
    asum = 0.0
    for (Z, mass_frac) in composition
        A = get(atomic_masses, Z, Float64(Z) * 2.0)
        zsum += mass_frac * Z / A
        asum += mass_frac / A
    end
    return zsum / asum
end

function _compute_mean_excitation_energy(composition::Dict{Int, Float64})
    I_values = Dict(
        1=>19.2, 6=>81.0, 7=>82.0, 8=>95.0, 11=>149.0, 12=>156.0,
        15=>173.0, 16=>180.0, 17=>174.0, 19=>190.0, 20=>191.0, 26=>286.0, 53=>491.0,
    )
    atomic_masses = Dict(
        1=>1.008, 6=>12.011, 7=>14.007, 8=>15.999, 11=>22.990, 12=>24.305,
        15=>30.974, 16=>32.06, 17=>35.45, 19=>39.098, 20=>40.078, 26=>55.845, 53=>126.904,
    )
    log_I_sum = 0.0
    ZA_sum = 0.0
    for (Z, mass_frac) in composition
        A = get(atomic_masses, Z, Float64(Z) * 2.0)
        I = get(I_values, Z, 10.0 * Z)
        ZA = mass_frac * Z / A
        log_I_sum += ZA * log(I)
        ZA_sum += ZA
    end
    return exp(log_I_sum / ZA_sum) * u"eV"
end

function load_xcat_materials_from_xlsx(xlsx_path::AbstractString)
    xf = XLSX.readxlsx(xlsx_path)
    sheet = xf["Sheet1"]
    data = sheet["A2:P34"]
    materials = Dict{Int, XA.Material}()
    for i in 1:size(data, 1)
        name = String(data[i, 1])
        organ_id = Int(data[i, 16])
        density = Float64(data[i, 15]) * u"g/cm^3"
        comp = Dict{Int, Float64}(
            1 => Float64(data[i, 2]),
            6 => Float64(data[i, 3]),
            7 => Float64(data[i, 4]),
            8 => Float64(data[i, 5]),
            11 => Float64(data[i, 6]),
            12 => Float64(data[i, 7]),
            15 => Float64(data[i, 8]),
            16 => Float64(data[i, 9]),
            17 => Float64(data[i, 10]),
            19 => Float64(data[i, 11]),
            20 => Float64(data[i, 12]),
            26 => Float64(data[i, 13]),
            53 => Float64(data[i, 14]),
        )
        filter!(p -> p.second > 0.0, comp)
        ZA = _compute_ZA_ratio(comp)
        I = _compute_mean_excitation_energy(comp)
        materials[organ_id] = XA.Material(name, ZA, I, density, comp)
    end
    return materials
end

function xcat_label_bbox(mask::AbstractArray{<:Unsigned,3}, labels::AbstractVector{<:Integer})
    wanted = Set(UInt16.(labels))
    minx = typemax(Int); miny = typemax(Int); minz = typemax(Int)
    maxx = 0; maxy = 0; maxz = 0
    found = false
    for k in axes(mask, 3), j in axes(mask, 2), i in axes(mask, 1)
        if UInt16(mask[i, j, k]) in wanted
            found = true
            minx = min(minx, i); miny = min(miny, j); minz = min(minz, k)
            maxx = max(maxx, i); maxy = max(maxy, j); maxz = max(maxz, k)
        end
    end
    found || error("Labels $labels not found in raw mask")
    return (min=(minx, miny, minz), max=(maxx, maxy, maxz))
end

function _bbox_center_mm(bounds)
    lo, hi = bounds
    return (
        0.5 * (lo[1] + hi[1]),
        0.5 * (lo[2] + hi[2]),
        0.5 * (lo[3] + hi[3]),
    )
end

function _union_bounds_mm(surfaces::AbstractVector{XCATNurbsSurface})
    lo = (Inf, Inf, Inf)
    hi = (-Inf, -Inf, -Inf)
    for surface in surfaces
        slo, shi = xcat_bounds(surface)
        lo = (min(lo[1], slo[1]), min(lo[2], slo[2]), min(lo[3], slo[3]))
        hi = (max(hi[1], shi[1]), max(hi[2], shi[2]), max(hi[3], shi[3]))
    end
    return lo, hi
end

function estimate_xcat_raw_alignment(
    raw::XCATRawLabels,
    surfaces::AbstractVector{XCATNurbsSurface};
    coronary_surface_names::Vector{String}=String[
        "dias_lad1", "dias_lad2", "dias_lad3", "dias_lcx", "dias_rca1", "dias_rca2",
    ],
    coronary_labels::Vector{Int}=[26],
    pericardium_surface_name::String="dias_pericardium",
    pericardium_labels::Vector{Int}=[29],
)
    object_map = xcat_object_dict(surfaces)
    coronary_surfaces = [object_map[name] for name in coronary_surface_names if haskey(object_map, name)]
    isempty(coronary_surfaces) && error("No coronary surfaces found for alignment.")
    haskey(object_map, pericardium_surface_name) || error("Missing pericardium surface `$pericardium_surface_name`.")

    coronary_bounds = _union_bounds_mm(coronary_surfaces)
    pericardium_bounds = xcat_bounds(object_map[pericardium_surface_name])
    coronary_center_mm = _bbox_center_mm(coronary_bounds)
    pericardium_center_mm = _bbox_center_mm(pericardium_bounds)

    raw_cor = xcat_label_bbox(raw.mask, coronary_labels)
    raw_peri = xcat_label_bbox(raw.mask, pericardium_labels)
    raw_cor_idx = (
        0.5 * (raw_cor.min[1] + raw_cor.max[1]),
        0.5 * (raw_cor.min[2] + raw_cor.max[2]),
        0.5 * (raw_cor.min[3] + raw_cor.max[3]),
    )
    raw_peri_idx = (
        0.5 * (raw_peri.min[1] + raw_peri.max[1]),
        0.5 * (raw_peri.min[2] + raw_peri.max[2]),
        0.5 * (raw_peri.min[3] + raw_peri.max[3]),
    )

    spacing = raw.voxel_size_mm
    cor_origin = (
        coronary_center_mm[1] - (raw_cor_idx[1] - 1.0) * spacing[1],
        coronary_center_mm[2] - (raw_cor_idx[2] - 1.0) * spacing[2],
        coronary_center_mm[3] - (raw_cor_idx[3] - 1.0) * spacing[3],
    )
    peri_origin = (
        pericardium_center_mm[1] - (raw_peri_idx[1] - 1.0) * spacing[1],
        pericardium_center_mm[2] - (raw_peri_idx[2] - 1.0) * spacing[2],
        pericardium_center_mm[3] - (raw_peri_idx[3] - 1.0) * spacing[3],
    )
    origin = (
        0.5 * (cor_origin[1] + peri_origin[1]),
        0.5 * (cor_origin[2] + peri_origin[2]),
        0.5 * (cor_origin[3] + peri_origin[3]),
    )

    return XCATRawAlignment(origin, raw.voxel_size_mm, coronary_center_mm, pericardium_center_mm, raw_cor_idx, raw_peri_idx)
end

function _segment_voxel_hits_raw(
    px::Float64, py::Float64, pz::Float64,
    dx::Float64, dy::Float64, dz::Float64,
    radius_mm::Float64,
    origin_mm::NTuple{3, Float64},
    spacing_mm::NTuple{3, Float64},
    dims::NTuple{3, Int};
    supersample::Int=1,
)
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
    total_sub = supersample^3
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
            push!(hits, (Int32(lin[ix, iy, iz]), Float32(inside / total_sub)))
        end
    end
    return hits
end

function _quantize_index(value::Float64, bins::AbstractVector{<:Real})
    vals = Float64.(bins)
    idx = argmin(abs.(vals .- value))
    return idx, vals[idx]
end

function _material_summary_row(label::Int, material::XA.Material)
    return (label=label, name=material.name, density_g_cm3=ustrip(u"g/cm^3", material.density))
end

function compile_basis_material_frame(
    raw_mask::AbstractArray{<:Unsigned,3},
    base_materials::AbstractDict{<:Integer, <:XA.Material},
    forest::CoronaryForest,
    results::Dict{String, <:ContrastTransportResult},
    frame_index::Int;
    fixed_prefix_segments::Dict{String, Int}=Dict{String, Int}(),
    alignment::XCATRawAlignment,
    blood_fraction_bins::AbstractVector{<:Real}=[0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0],
    iodine_concentration_bins_mg_mL::AbstractVector{<:Real}=[0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 12.5, 15.0],
    min_fraction::Float64=1e-4,
    min_radius_um::Float64=0.0,
    supersample::Int=1,
    blood_material::XA.Material=XA.Materials.blood,
    iodine_material_builder::Function=_iodine_contrast_material,
    mixture_builder::Function=_create_material_mixture,
)
    dims = size(raw_mask)
    spacing = alignment.voxel_size_mm
    origin = alignment.origin_mm
    union_prod = Dict{Int32, Float64}()
    conc_num = Dict{Int32, Float64}()
    conc_den = Dict{Int32, Float64}()

    for name in sort!(collect(keys(forest.trees)))
        tree = forest.trees[name]
        result = results[name]
        seg = tree.segments
        start_seg = get(fixed_prefix_segments, name, 0) + 1
        for seg_id in start_seg:seg.n
            radius_um = seg.radius[seg_id] * 2000.0
            radius_um < min_radius_um && continue
            c_seg = Float64(result.concentration[seg_id, frame_index])
            hits = _segment_voxel_hits_raw(
                seg.proximal_x[seg_id], seg.proximal_y[seg_id], seg.proximal_z[seg_id],
                seg.distal_x[seg_id], seg.distal_y[seg_id], seg.distal_z[seg_id],
                seg.radius[seg_id],
                origin,
                spacing,
                dims;
                supersample=supersample,
            )
            for (linear_idx, frac) in hits
                f = Float64(frac)
                union_prod[linear_idx] = get(union_prod, linear_idx, 1.0) * (1.0 - f)
                conc_num[linear_idx] = get(conc_num, linear_idx, 0.0) + f * c_seg
                conc_den[linear_idx] = get(conc_den, linear_idx, 0.0) + f
            end
        end
    end

    overlay_rows = Tuple{Int32, Int, Float64, Float64}[]
    for linear_idx in keys(union_prod)
        phi = 1.0 - union_prod[linear_idx]
        phi < min_fraction && continue
        c = conc_den[linear_idx] > 0.0 ? conc_num[linear_idx] / conc_den[linear_idx] : 0.0
        host_label = Int(raw_mask[Int(linear_idx)])
        push!(overlay_rows, (linear_idx, host_label, phi, c))
    end

    material_dict = Dict{Int, XA.Material}(base_materials)
    dynamic_labels = Dict{Tuple{Int, Int, Int}, Int}()
    next_label = maximum(keys(material_dict)) + 1

    for (_, host_label, phi, c) in overlay_rows
        haskey(base_materials, host_label) || error("Missing host material for label $host_label")
        phi_idx, phi_bin = _quantize_index(phi, blood_fraction_bins)
        c_idx, c_bin = _quantize_index(c, iodine_concentration_bins_mg_mL)
        key = (host_label, phi_idx, c_idx)
        if !haskey(dynamic_labels, key)
            iodinated_blood = iodine_material_builder(
                blood_material,
                c_bin;
                name="blood_iodine_$(round(c_bin, digits=3))mgmL",
            )
            host_material = base_materials[host_label]
            mixed = mixture_builder(
                XA.Material[host_material, iodinated_blood],
                Float64[1.0 - phi_bin, phi_bin];
                by_volume=true,
                name="$(host_material.name)_blood$(round(phi_bin, digits=4))_iodine$(round(c_bin, digits=3))",
            )
            dynamic_labels[key] = next_label
            material_dict[next_label] = mixed
            next_label += 1
        end
    end

    max_label = max(maximum(raw_mask), maximum(keys(material_dict)))
    frame_mask = max_label <= typemax(UInt8) ? copy(UInt8.(raw_mask)) : UInt16.(raw_mask)
    for (linear_idx, host_label, phi, c) in overlay_rows
        phi_idx, _ = _quantize_index(phi, blood_fraction_bins)
        c_idx, _ = _quantize_index(c, iodine_concentration_bins_mg_mL)
        frame_mask[Int(linear_idx)] = eltype(frame_mask)(dynamic_labels[(host_label, phi_idx, c_idx)])
    end

    rows = [_material_summary_row(label, material_dict[label]) for label in sort!(collect(keys(material_dict)))]
    return (
        mask=frame_mask,
        materials_dict=material_dict,
        material_rows=rows,
        overlay_voxels=length(overlay_rows),
        dynamic_label_count=length(dynamic_labels),
    )
end

function save_basis_frame(
    output_dir::AbstractString,
    frame_mask::AbstractArray{<:Unsigned,3},
    materials_dict::AbstractDict{<:Integer, <:XA.Material};
    time_s::Real,
    voxel_size_cm::NTuple{3, <:Real},
    origin_mm::NTuple{3, <:Real},
    prefix::AbstractString="xcat_basis",
)
    isdir(output_dir) || mkpath(output_dir)
    stem = joinpath(output_dir, "$(prefix)_t$(lpad(round(Int, time_s), 3, '0'))")
    raw_ext = eltype(frame_mask) == UInt8 ? ".raw" : ".raw16"
    raw_path = stem * raw_ext
    open(raw_path, "w") do io
        write(io, frame_mask)
    end

    materials_path = stem * "_materials.jld2"
    JLD2.jldsave(materials_path; materials_dict=materials_dict)

    csv_path = stem * "_materials.csv"
    open(csv_path, "w") do io
        println(io, "label,name,density_g_cm3")
        for label in sort!(collect(keys(materials_dict)))
            mat = materials_dict[label]
            println(io, label, ",", replace(mat.name, ","=>"_"), ",", ustrip(u"g/cm^3", mat.density))
        end
    end

    manifest_path = stem * ".json"
    manifest = Dict(
        "time_s" => Float64(time_s),
        "raw_path" => raw_path,
        "materials_jld2" => materials_path,
        "materials_csv" => csv_path,
        "dims" => collect(size(frame_mask)),
        "dtype" => string(eltype(frame_mask)),
        "voxel_size_cm" => collect(Float64.(voxel_size_cm)),
        "origin_mm" => collect(Float64.(origin_mm)),
        "n_materials" => length(materials_dict),
    )
    open(manifest_path, "w") do io
        _write_json(io, manifest)
    end

    return (raw_path=raw_path, materials_path=materials_path, csv_path=csv_path, manifest_path=manifest_path)
end

function generate_basis_frames(
    output_dir::AbstractString,
    raw::XCATRawLabels,
    base_materials::AbstractDict{<:Integer, <:XA.Material},
    forest::CoronaryForest,
    results::Dict{String, <:ContrastTransportResult};
    frame_times_s::AbstractVector{<:Real},
    fixed_prefix_segments::Dict{String, Int}=Dict{String, Int}(),
    alignment::XCATRawAlignment,
    blood_fraction_bins::AbstractVector{<:Real}=[0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0],
    iodine_concentration_bins_mg_mL::AbstractVector{<:Real}=[0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 12.5, 15.0],
    min_fraction::Float64=1e-4,
    min_radius_um::Float64=0.0,
    supersample::Int=1,
    prefix::AbstractString="xcat_basis",
)
    times = results[first(sort!(collect(keys(results))))].times
    voxel_size_cm = ntuple(i -> raw.voxel_size_mm[i] / 10.0, 3)
    manifests = Any[]
    for time_s in frame_times_s
        idx = argmin(abs.(times .- Float64(time_s)))
        compiled = compile_basis_material_frame(
            raw.mask,
            base_materials,
            forest,
            results,
            idx;
            fixed_prefix_segments=fixed_prefix_segments,
            alignment=alignment,
            blood_fraction_bins=blood_fraction_bins,
            iodine_concentration_bins_mg_mL=iodine_concentration_bins_mg_mL,
            min_fraction=min_fraction,
            min_radius_um=min_radius_um,
            supersample=supersample,
        )
        saved = save_basis_frame(
            output_dir,
            compiled.mask,
            compiled.materials_dict;
            time_s=Float64(times[idx]),
            voxel_size_cm=voxel_size_cm,
            origin_mm=alignment.origin_mm,
            prefix=prefix,
        )
        push!(manifests, Dict(
            "time_s" => Float64(times[idx]),
            "manifest_path" => saved.manifest_path,
            "overlay_voxels" => compiled.overlay_voxels,
            "dynamic_label_count" => compiled.dynamic_label_count,
        ))
    end
    run_manifest = joinpath(output_dir, "$(prefix)_frames.json")
    open(run_manifest, "w") do io
        _write_json(io, Dict(
            "source_raw" => raw.source_path,
            "frame_manifests" => manifests,
            "origin_mm" => collect(alignment.origin_mm),
            "voxel_size_cm" => collect(voxel_size_cm),
        ))
    end
    return run_manifest
end
