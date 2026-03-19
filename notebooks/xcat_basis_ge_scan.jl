### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00001
begin
    VERSION >= v"1.11" || error("xcat_basis_ge_scan.jl requires Julia 1.11+ because BasisSimulator has compat julia = 1.11")

    import Pkg
    env_dir = joinpath(dirname(@__DIR__), ".basis_scan_env")
    v4_dir = dirname(@__DIR__)
    basis_dir = normpath(joinpath(v4_dir, "..", "basis_simulator"))
    skip_bootstrap = lowercase(get(ENV, "XCAT_GE_SKIP_BOOTSTRAP", "false")) in ("1", "true", "yes")

    Pkg.activate(env_dir)
    if !skip_bootstrap
        Pkg.develop([Pkg.PackageSpec(path=v4_dir), Pkg.PackageSpec(path=basis_dir)])
        Pkg.add([Pkg.PackageSpec(name="JLD2"), Pkg.PackageSpec(name="JSON"), Pkg.PackageSpec(name="CUDA")])
        Pkg.instantiate()
    end
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00002
using VesselTree

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00003
import BasisSimulator as BS

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00004
using JLD2

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00005
import JSON

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00006
using Dates

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00007
using Statistics

const _CUDA = try
    @eval using CUDA
    CUDA
catch
    nothing
end

function to_accel_array(x::AbstractArray)
    if _CUDA !== nothing && CUDA.functional()
        return CUDA.CuArray(x)
    end
    return x
end

function clear_accel!()
    GC.gc(true)
    if _CUDA !== nothing && CUDA.functional()
        CUDA.reclaim()
    end
    return nothing
end


# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00009
md"""
# XCAT Basis Frames -> GE Scanner

This notebook consumes the time-resolved phantom frames compiled from `VesselTree` and runs them through the existing **BasisSimulator GE workflow** without modifying the simulator.

The default path is:
1. load one or more precompiled XCAT basis frames,
2. crop around the dynamic coronary overlay,
3. optionally downsample for faster simulation,
4. build a `BS.Phantom`,
5. run a GE Revolution Apex style single-kVp simulation,
6. reconstruct FDK output and save attenuation / HU volumes.

Use **Julia 1.12** (or any `>= 1.11`) for this notebook.
"""

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00010
begin
    const DEFAULT_FRAME_RUN_MANIFEST = joinpath(dirname(@__DIR__), "output", "basis_frames", "20260318T125831", "xcat_basis_frames.json")
    const DEFAULT_OUTPUT_DIR = joinpath(dirname(@__DIR__), "output", "ge_scans")
    const ENV_TRUE = ("1", "true", "yes")
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00011
begin
    run_scan = lowercase(get(ENV, "XCAT_RUN_GE_SCAN", "false")) in ENV_TRUE
    scan_all_frames = lowercase(get(ENV, "XCAT_SCAN_ALL_FRAMES", "false")) in ENV_TRUE
    calibrate_water = lowercase(get(ENV, "XCAT_GE_CALIBRATE_WATER", "true")) in ENV_TRUE
    frame_run_manifest_path = get(ENV, "XCAT_FRAME_RUN_MANIFEST", DEFAULT_FRAME_RUN_MANIFEST)
    selected_frame_index = parse(Int, get(ENV, "XCAT_FRAME_INDEX", "1"))
    scan_all_frames = lowercase(get(ENV, "XCAT_SCAN_ALL_FRAMES", "false")) in ("1", "true", "yes")
    crop_padding_cm = parse(Float64, get(ENV, "XCAT_GE_CROP_PADDING_CM", "3.0"))
    downsample_factor = parse(Int, get(ENV, "XCAT_GE_DOWNSAMPLE_FACTOR", "4"))
    recon_xy_cap = parse(Int, get(ENV, "XCAT_GE_RECON_XY_CAP", "512"))
    recon_slices_cap = parse(Int, get(ENV, "XCAT_GE_RECON_SLICES_CAP", "128"))
    protocol_collimation_mm = parse(Float64, get(ENV, "XCAT_GE_PROTOCOL_COLLIMATION_MM", "80.0"))
    recon_slice_thickness_mm = parse(Float64, get(ENV, "XCAT_GE_RECON_SLICE_THICKNESS_MM", "0.625"))
    protocol_additional_al_mm = parse(Float64, get(ENV, "XCAT_GE_PROTOCOL_ADDITIONAL_AL_MM", "4.5"))
    output_dir = get(ENV, "XCAT_GE_SCAN_OUTPUT_DIR", DEFAULT_OUTPUT_DIR)
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00012
begin
    scanner_source_to_isocenter_mm = 625.6
    scanner_source_to_detector_mm = 1100.0
    scanner_detector_cols = 834
    scanner_detector_rows = 256
    scanner_detector_row_size_mm = 0.625
    scanner_detector_col_size_mm = 0.6
    scanner_target_angle_deg = 10.0
    scanner_flat_filter_thickness_mm = 2.5
    scanner_bowtie_filter = Symbol(get(ENV, "XCAT_GE_BOWTIE_FILTER", "ge_revolution_large"))
    scanner_detector_material = Symbol(get(ENV, "XCAT_GE_DETECTOR_MATERIAL", "lumex"))
    scanner_detection_gain = parse(Float64, get(ENV, "XCAT_GE_DETECTION_GAIN", "10.0"))
    scanner_electronic_noise = parse(Float64, get(ENV, "XCAT_GE_ELECTRONIC_NOISE", "0.0"))

    protocol_kvp = parse(Float64, get(ENV, "XCAT_GE_PROTOCOL_KVP", "120.0"))
    protocol_ma = parse(Float64, get(ENV, "XCAT_GE_PROTOCOL_MA", "150.0"))
    protocol_views = parse(Int, get(ENV, "XCAT_GE_PROTOCOL_VIEWS", "984"))
    protocol_rotation_time_s = parse(Float64, get(ENV, "XCAT_GE_PROTOCOL_ROTATION_TIME_S", "1.0"))

    sim_fidelity = Symbol(get(ENV, "XCAT_GE_SIM_FIDELITY", "high"))
    recon_filter = Symbol(get(ENV, "XCAT_GE_RECON_FILTER", "standard"))
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00020
function load_basis_run_manifest(path::AbstractString)
    return JSON.parsefile(path)
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00021
function load_basis_frame(frame_manifest_path::AbstractString)
    manifest = JSON.parsefile(frame_manifest_path)
    raw_path = manifest["raw_path"]
    dims = Tuple(Int.(manifest["dims"]))
    dtype = manifest["dtype"] == "UInt16" ? UInt16 : UInt8
    voxel_size_cm = Tuple(Float64.(manifest["voxel_size_cm"]))
    origin_cm = Tuple(Float64.(manifest["origin_mm"]) ./ 10.0)

    io = open(raw_path, "r")
    data = read!(io, Vector{dtype}(undef, prod(dims)))
    close(io)
    mask = reshape(data, dims)
    materials = JLD2.load(manifest["materials_jld2"], "materials_dict")

    return (manifest=manifest, mask=mask, materials_dict=materials, voxel_size_cm=voxel_size_cm, origin_cm=origin_cm)
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00022
function dynamic_label_ids(materials_dict)
    labels = Int[]
    for (label, material) in materials_dict
        name = lowercase(String(material.name))
        if occursin("_blood", name) || occursin("iodine", name)
            push!(labels, Int(label))
        end
    end
    sort!(unique!(labels))
    return labels
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00023
function crop_mask_to_labels(mask::AbstractArray{<:Unsigned, 3}, labels::AbstractVector{<:Integer}, voxel_size_cm::NTuple{3, Float64}, origin_cm::NTuple{3, Float64}; padding_cm::Float64=3.0)
    isempty(labels) && return (mask=mask, origin_cm=origin_cm, bbox=(axes(mask,1), axes(mask,2), axes(mask,3)))
    wanted = Set(eltype(mask).(labels))
    found = false
    imin = typemax(Int); jmin = typemax(Int); kmin = typemax(Int)
    imax = 0; jmax = 0; kmax = 0
    for k in axes(mask, 3), j in axes(mask, 2), i in axes(mask, 1)
        if mask[i, j, k] in wanted
            found = true
            imin = min(imin, i); jmin = min(jmin, j); kmin = min(kmin, k)
            imax = max(imax, i); jmax = max(jmax, j); kmax = max(kmax, k)
        end
    end
    found || return (mask=mask, origin_cm=origin_cm, bbox=(axes(mask,1), axes(mask,2), axes(mask,3)))

    pad_i = ceil(Int, padding_cm / voxel_size_cm[1])
    pad_j = ceil(Int, padding_cm / voxel_size_cm[2])
    pad_k = ceil(Int, padding_cm / voxel_size_cm[3])

    ir = max(first(axes(mask, 1)), imin - pad_i):min(last(axes(mask, 1)), imax + pad_i)
    jr = max(first(axes(mask, 2)), jmin - pad_j):min(last(axes(mask, 2)), jmax + pad_j)
    kr = max(first(axes(mask, 3)), kmin - pad_k):min(last(axes(mask, 3)), kmax + pad_k)

    cropped = mask[ir, jr, kr]
    new_origin_cm = (
        origin_cm[1] + (first(ir) - 1) * voxel_size_cm[1],
        origin_cm[2] + (first(jr) - 1) * voxel_size_cm[2],
        origin_cm[3] + (first(kr) - 1) * voxel_size_cm[3],
    )
    return (mask=cropped, origin_cm=new_origin_cm, bbox=(ir, jr, kr))
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00024
function downsample_mask_nn(mask::AbstractArray{T, 3}, voxel_size_cm::NTuple{3, Float64}, origin_cm::NTuple{3, Float64}, factor::Int) where {T<:Unsigned}
    factor == 1 && return (mask=mask, voxel_size_cm=voxel_size_cm, origin_cm=origin_cm)
    old_size = size(mask)
    new_size = old_size .÷ factor
    result = Array{T}(undef, new_size...)
    offset = factor ÷ 2
    for k in 1:new_size[3], j in 1:new_size[2], i in 1:new_size[1]
        oi = min((i - 1) * factor + offset + 1, old_size[1])
        oj = min((j - 1) * factor + offset + 1, old_size[2])
        ok = min((k - 1) * factor + offset + 1, old_size[3])
        result[i, j, k] = mask[oi, oj, ok]
    end
    new_voxel = (voxel_size_cm[1] * factor, voxel_size_cm[2] * factor, voxel_size_cm[3] * factor)
    new_origin = (
        origin_cm[1] + offset * voxel_size_cm[1],
        origin_cm[2] + offset * voxel_size_cm[2],
        origin_cm[3] + offset * voxel_size_cm[3],
    )
    return (mask=result, voxel_size_cm=new_voxel, origin_cm=new_origin)
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00025
function build_ge_apex_scanner()
    return BS.Scanner(
        source_to_isocenter=scanner_source_to_isocenter_mm,
        source_to_detector=scanner_source_to_detector_mm,
        detector_rows=scanner_detector_rows,
        detector_cols=scanner_detector_cols,
        detector_row_size=scanner_detector_row_size_mm,
        detector_col_size=scanner_detector_col_size_mm,
        detector_shape=BS.CURVED_DETECTOR,
        focal_spot_width=1.0,
        focal_spot_length=1.0,
        target_angle=scanner_target_angle_deg,
        flat_filter_material=:aluminum,
        flat_filter_thickness=scanner_flat_filter_thickness_mm,
        bowtie_filter=scanner_bowtie_filter,
        detector_material=scanner_detector_material,
        detector_depth=3.0,
        fill_factor_row=0.9,
        fill_factor_col=0.9,
        electronic_noise=scanner_electronic_noise,
        detection_gain=scanner_detection_gain,
    )
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00026
function build_ge_protocol(n_recon_slices::Int)
    filters = protocol_additional_al_mm > 0 ? [("Al", protocol_additional_al_mm)] : Tuple{String, Float64}[]
    return BS.CTProtocol(
        kVp=protocol_kvp,
        mA=protocol_ma,
        views=protocol_views,
        rotation_time=protocol_rotation_time_s,
        collimation_mm=protocol_collimation_mm,
        additional_filters=filters,
    )
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00027
build_sim_options() = BS.SimOptions(fidelity=sim_fidelity)

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00028
function build_recon_options(mask::AbstractArray, voxel_size_cm::NTuple{3, Float64})
    recon_xy = min(recon_xy_cap, min(size(mask, 1), size(mask, 2)))
    target_recon_slices = round(Int, protocol_collimation_mm / recon_slice_thickness_mm)
    recon_slices = min(recon_slices_cap, target_recon_slices, size(mask, 3))
    fov_cm = max(size(mask, 1) * voxel_size_cm[1], size(mask, 2) * voxel_size_cm[2])
    z_cm = recon_slices * recon_slice_thickness_mm / 10.0
    return BS.ReconOptions(
        algorithm=:fdk,
        matrix_size=(recon_xy, recon_xy, recon_slices),
        fov_cm=fov_cm,
        z_cm=z_cm,
        filter=recon_filter,
    )
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00029
function build_water_phantom(scanner::BS.Scanner)
    water_size = (256, 256, 8)
    water_voxel_cm = (0.05, 0.05, scanner.detector_row_size / 10.0)
    water_mask = zeros(UInt8, water_size...)
    cx, cy = water_size[1] ÷ 2, water_size[2] ÷ 2
    r = 100
    for i in 1:water_size[1], j in 1:water_size[2]
        if (i - cx)^2 + (j - cy)^2 <= r^2
            water_mask[i, j, :] .= UInt8(1)
        end
    end
    materials = Dict(0 => BS.XA.Materials.air, 1 => BS.XA.Materials.water)
    phantom = BS.Phantom(to_accel_array(water_mask), materials, water_voxel_cm)
    recon = BS.ReconOptions(
        algorithm=:fdk,
        matrix_size=(256, 256, 8),
        fov_cm=25.0,
        z_cm=8 * scanner.detector_row_size / 10.0,
        filter=:standard,
    )
    return phantom, recon
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00030
function extract_water_mu(vol)
    nx, ny, nz = size(vol)
    cx, cy = nx ÷ 2, ny ÷ 2
    r = max(2, nx ÷ 10)
    z_start = max(1, nz ÷ 4)
    z_end = min(nz, 3 * nz ÷ 4)
    vals = Float64[]
    for k in z_start:z_end, j in (cy - r):(cy + r), i in (cx - r):(cx + r)
        if 1 <= i <= nx && 1 <= j <= ny && (i - cx)^2 + (j - cy)^2 <= r^2
            push!(vals, vol[i, j, k])
        end
    end
    return mean(vals)
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00031
function calibrate_ge_water_mu(scanner::BS.Scanner, protocol::BS.CTProtocol, sim_opts::BS.SimOptions)
    water_phantom, water_recon = build_water_phantom(scanner)
    ws = BS.create_eict_workspace(scanner, protocol, sim_opts, water_recon, water_phantom)
    BS.simulate!(ws, water_phantom, scanner, protocol, sim_opts, water_recon)
    ws_fdk = BS.create_fdk_recon_workspace(ws.sino_noisy_out, ws.geom, water_recon.matrix_size; filter=BS.StandardFilter())
    vol = Array(BS.reconstruct!(ws_fdk, ws.sino_noisy_out, ws.geom, water_recon.matrix_size))
    return extract_water_mu(vol)
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00032
function save_ge_scan_result(output_dir::AbstractString, frame_name::AbstractString, result)
    isdir(output_dir) || mkpath(output_dir)
    stem = joinpath(output_dir, frame_name)
    jld2_path = stem * "_scan.jld2"
    manifest_path = stem * "_scan.json"
    JLD2.jldsave(jld2_path;
        recon_mu=result.recon_mu,
        recon_hu=result.recon_hu,
        mu_water=result.μ_water,
        frame_manifest_path=result.frame_manifest_path,
        crop_bbox_ijk=result.crop_bbox_ijk,
        voxel_size_cm=result.voxel_size_cm,
        origin_cm=result.origin_cm,
        recon_matrix_size=result.recon_opts.matrix_size,
    )
    open(manifest_path, "w") do io
        JSON.print(io, Dict(
            "frame_manifest_path" => result.frame_manifest_path,
            "scan_path" => jld2_path,
            "time_s" => result.time_s,
            "mu_water" => result.μ_water,
            "recon_matrix_size" => collect(result.recon_opts.matrix_size),
            "voxel_size_cm" => collect(result.voxel_size_cm),
            "origin_cm" => collect(result.origin_cm),
            "crop_bbox_ijk" => [collect(result.crop_bbox_ijk[1]), collect(result.crop_bbox_ijk[2]), collect(result.crop_bbox_ijk[3])],
        ))
    end
    return (jld2_path=jld2_path, manifest_path=manifest_path)
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00033
function scan_frame(frame_manifest_path::AbstractString; crop_padding_cm::Float64=3.0, downsample_factor::Int=4, calibrate_water::Bool=true, water_mu_override::Union{Nothing, Float64}=nothing)
    loaded = load_basis_frame(frame_manifest_path)
    dynamic_labels = dynamic_label_ids(loaded.materials_dict)
    cropped = crop_mask_to_labels(loaded.mask, dynamic_labels, loaded.voxel_size_cm, loaded.origin_cm; padding_cm=crop_padding_cm)
    ds = downsample_mask_nn(cropped.mask, loaded.voxel_size_cm, cropped.origin_cm, downsample_factor)

    println("Scanning frame: ", frame_manifest_path)
    println("  original dims = ", size(loaded.mask), ", cropped dims = ", size(cropped.mask), ", downsampled dims = ", size(ds.mask))
    println("  voxel_size_cm = ", ds.voxel_size_cm, ", origin_cm = ", ds.origin_cm)

    mask_accel = to_accel_array(ds.mask)
    phantom = BS.Phantom(mask_accel, Dict(Int(k) => v for (k, v) in loaded.materials_dict), ds.voxel_size_cm; origin=ds.origin_cm)
    scanner = build_ge_apex_scanner()
    recon_opts = build_recon_options(ds.mask, ds.voxel_size_cm)
    protocol = build_ge_protocol(recon_opts.matrix_size[3])
    sim_opts = build_sim_options()

    println("  recon matrix = ", recon_opts.matrix_size, ", fov_cm = ", recon_opts.fov_cm, ", z_cm = ", recon_opts.z_cm)
    println("  protocol views = ", protocol.views, ", kVp = ", protocol.kVp, ", mA = ", protocol.mA, ", fidelity = ", sim_opts.fidelity)

    μ_water = isnothing(water_mu_override) ? (calibrate_water ? calibrate_ge_water_mu(scanner, protocol, sim_opts) : nothing) : water_mu_override
    calibrate_water && println("  calibrated water mu = ", μ_water)

    ws = BS.create_eict_workspace(scanner, protocol, sim_opts, recon_opts, phantom)
    println("  EICT workspace ready")
    BS.simulate!(ws, phantom, scanner, protocol, sim_opts, recon_opts)
    println("  forward simulation complete")
    ws_fdk = BS.create_fdk_recon_workspace(ws.sino_noisy_out, ws.geom, recon_opts.matrix_size; filter=BS.StandardFilter())
    recon_mu = Array(BS.reconstruct!(ws_fdk, ws.sino_noisy_out, ws.geom, recon_opts.matrix_size))
    println("  reconstruction complete")
    recon_hu = isnothing(μ_water) ? nothing : Array(BS.to_hounsfield(recon_mu; μ_water=μ_water))

    ws = nothing
    ws_fdk = nothing
    phantom = nothing
    mask_accel = nothing
    clear_accel!()

    return (
        frame_manifest_path=frame_manifest_path,
        time_s=loaded.manifest["time_s"],
        recon_mu=recon_mu,
        recon_hu=recon_hu,
        μ_water=μ_water,
        scanner=scanner,
        protocol=protocol,
        sim_opts=sim_opts,
        recon_opts=recon_opts,
        voxel_size_cm=ds.voxel_size_cm,
        origin_cm=ds.origin_cm,
        crop_bbox_ijk=cropped.bbox,
        dynamic_labels=dynamic_labels,
    )
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00040
frame_run_manifest = load_basis_run_manifest(frame_run_manifest_path)

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00041
selected_frame_manifests = begin
    frames = frame_run_manifest["frame_manifests"]
    if scan_all_frames
        [String(frame["manifest_path"]) for frame in frames]
    else
        idx = clamp(selected_frame_index, 1, length(frames))
        [String(frames[idx]["manifest_path"])]
    end
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00042
scan_state = if run_scan
    run_stamp = Dates.format(Dates.now(), dateformat"yyyymmddTHHMMSS")
    run_dir = joinpath(output_dir, run_stamp)
    saved = []
    cached_water_mu = Ref{Union{Nothing, Float64}}(nothing)
    for frame_manifest_path in selected_frame_manifests
        result = scan_frame(
            frame_manifest_path;
            crop_padding_cm=crop_padding_cm,
            downsample_factor=downsample_factor,
            calibrate_water=calibrate_water,
            water_mu_override=cached_water_mu[],
        )
        if calibrate_water && isnothing(cached_water_mu[])
            cached_water_mu[] = result.μ_water
        end
        frame_name = splitext(basename(frame_manifest_path))[1]
        saved_paths = save_ge_scan_result(run_dir, frame_name, result)
        push!(saved, Dict(
            "frame_manifest_path" => frame_manifest_path,
            "time_s" => result.time_s,
            "scan_manifest_path" => saved_paths.manifest_path,
            "scan_jld2_path" => saved_paths.jld2_path,
            "mu_water" => result.μ_water,
            "recon_matrix_size" => collect(result.recon_opts.matrix_size),
            "voxel_size_cm" => collect(result.voxel_size_cm),
        ))
    end
    run_manifest_path = joinpath(run_dir, "ge_scan_run.json")
    open(run_manifest_path, "w") do io
        JSON.print(io, Dict(
            "frame_run_manifest_path" => frame_run_manifest_path,
            "scan_manifests" => saved,
            "downsample_factor" => downsample_factor,
            "crop_padding_cm" => crop_padding_cm,
            "calibrate_water" => calibrate_water,
        ))
    end
    (; run_dir, run_manifest_path, saved)
else
    nothing
end

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00050
md"""
## Status

- `run_scan = $(run_scan)`
- `scan_all_frames = $(scan_all_frames)`
- selected frames: $(length(selected_frame_manifests))
- frame manifest set: `$(frame_run_manifest_path)`
- downsample factor: `$(downsample_factor)`
- crop padding: `$(crop_padding_cm) cm`
- water calibration: `$(calibrate_water)`
"""

# ╔═╡ 6c9ec5a6-2d4b-4fd5-8137-0c5e11c00051
if scan_state === nothing
    md"""
    Toggle `run_scan = true` or launch from the shell with:

    ```bash
    XCAT_RUN_GE_SCAN=true XCAT_FRAME_INDEX=1 XCAT_GE_DOWNSAMPLE_FACTOR=8 \
    /home/molloi-lab/.julia/juliaup/julia-1.12.5+0.x64.linux.gnu/bin/julia \
      --project=$(joinpath(dirname(@__DIR__), ".basis_scan_env")) \
      -e 'include("v4/notebooks/xcat_basis_ge_scan.jl")'
    ```
    """
else
    md"""
    **Scan run dir:** `$(scan_state.run_dir)`  
    **Run manifest:** `$(scan_state.run_manifest_path)`
    """
end

# ╔═╡ Cell order:
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00001
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00002
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00003
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00004
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00005
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00006
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00007
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00009
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00010
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00011
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00012
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00020
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00021
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00022
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00023
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00024
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00025
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00026
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00027
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00028
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00029
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00030
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00031
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00032
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00033
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00040
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00041
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00042
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00050
# ╠═6c9ec5a6-2d4b-4fd5-8137-0c5e11c00051
