### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0001
begin
    import Pkg
    Pkg.activate(dirname(@__DIR__))
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0002
using VesselTree

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0003
using Random

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0004
using Dates

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0005
md"""
# XCAT Basis Frame Compiler

This notebook prepares time-resolved phantom frames for `basis_simulator` without modifying the simulator itself.

For each requested time point we compile:
- a full 3D label mask on the original XCAT raw grid
- a `label => XA.Material` dictionary that includes dynamic iodinated-blood mixtures
- per-frame manifests that can be consumed downstream by the GE scanner workflow

The current implementation overlays **only the newly grown distal tree** on top of the baseline XCAT phantom.
"""

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0010
begin
    const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
    const DEFAULT_XCAT_RAW = "/home/molloi-lab/smb_mount/shared_drive/shu_nie/PVAT_Analysis/digital phantoms/vmale50_1600x1400x500_8bit_little_endian_act_1.raw"
    const DEFAULT_XCAT_RAW_LOG = "/home/molloi-lab/smb_mount/shared_drive/shu_nie/PVAT_Analysis/digital phantoms/vmale50_1600x1400x500_8bit_little_endian_log"
    const DEFAULT_XCAT_NONCONTRAST_XLSX = joinpath(dirname(@__DIR__), "..", "basis_simulator", "verification", "data", "xcat", "Material_Spreadsheets", "vmale_50_materials_heart_non_contrast.xlsx")
    const DEFAULT_OUTPUT_DIR = joinpath(dirname(@__DIR__), "output", "basis_frames")
    const MMHG_TO_PA = 133.32239
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0011
begin
    run_pipeline = lowercase(get(ENV, "XCAT_RUN_BASIS_PIPELINE", "false")) in ("1", "true", "yes")
    nrb_path = DEFAULT_XCAT_NRB
    raw_path = DEFAULT_XCAT_RAW
    raw_log_path = DEFAULT_XCAT_RAW_LOG
    materials_xlsx_path = DEFAULT_XCAT_NONCONTRAST_XLSX
    frame_times_s = collect(0.0:1.0:9.0)
    contrast_dt_s = 0.05
    contrast_t_end_s = 25.0
    output_dir = get(ENV, "XCAT_BASIS_OUTPUT_DIR", DEFAULT_OUTPUT_DIR)
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0012
begin
    target_terminals = Dict("LAD" => 8, "LCX" => 72, "RCA" => 24)
    prefix_blend = Dict("LAD" => 0.0, "LCX" => 0.5, "RCA" => 0.3)
    root_target_diam_um = Dict("LAD" => 3500.0, "LCX" => 4500.0, "RCA" => 4500.0)
    territory_fractions = Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35)
    contrast_peak_mg_mL = 15.0
    contrast_t0_s = 0.0
    contrast_tmax_s = 4.0
    contrast_alpha = 8.0
    blood_fraction_bins = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0]
    iodine_bins_mg_mL = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 12.5, 15.0]
    min_fraction = 1e-4
    min_radius_um = 0.0
    supersample = 1
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0020
function build_generated(nrb_path::AbstractString)
    return generate_xcat_kassab_coronary(
        nrb_path,
        kassab_coronary_params();
        phase="dias",
        rng=MersenneTwister(321),
        verbose=true,
        handoff_order=6,
        subdivision_max_order=11,
        apply_geometry=false,
        enforce_full_cutoff=true,
        max_tree_capacity=2_000_000,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        target_terminals=target_terminals,
        territory_fractions=territory_fractions,
    )
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0021
function build_hemodynamic_forest(generated)
    trees = Dict{String, VascularTree}()
    base_params = kassab_coronary_params()
    for name in ["LAD", "LCX", "RCA"]
        tree = deepcopy(generated.forest.trees[name])
        update_radii!(tree, base_params.gamma)
        src = generated.fixed_trees[name]
        target_diam_um = root_target_diam_um[name]
        scale = target_diam_um / (src.segments.radius[Int(src.root_segment_id)] * 2000.0)
        alpha = prefix_blend[name]
        for i in 1:min(src.segments.n, tree.segments.n)
            target_r = src.segments.radius[i] * scale
            tree.segments.radius[i] = (1 - alpha) * tree.segments.radius[i] + alpha * target_r
        end
        trees[name] = tree
    end

    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )
    return CoronaryForest(trees, generated.forest.territory_map, generated.forest.params, pressure_params)
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0022
function build_root_inputs(times)
    return Dict(name => gamma_variate_input(times; amplitude=contrast_peak_mg_mL, t0=contrast_t0_s, tmax=contrast_tmax_s, alpha=contrast_alpha) for name in ("LAD", "LCX", "RCA"))
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0030
workflow_state = if run_pipeline
    generated = build_generated(nrb_path)
    forest = build_hemodynamic_forest(generated)
    contrast_times = collect(0.0:contrast_dt_s:contrast_t_end_s)
    results = simulate_forest_contrast(
        forest,
        contrast_times;
        root_inputs=build_root_inputs(contrast_times),
        recompute_hemodynamics=true,
        storage_type=Float32,
        threaded=true,
    )

    raw = load_xcat_raw_labels(raw_path; log_path=raw_log_path)
    base_materials = load_xcat_materials_from_xlsx(materials_xlsx_path)
    alignment = estimate_xcat_raw_alignment(raw, generated.surfaces)
    fixed_prefix_segments = Dict(name => generated.fixed_trees[name].segments.n for name in keys(generated.fixed_trees))

    run_stamp = format_run_timestamp()
    frame_dir = joinpath(output_dir, run_stamp)
    manifest = generate_basis_frames(
        frame_dir,
        raw,
        base_materials,
        forest,
        results;
        frame_times_s=frame_times_s,
        fixed_prefix_segments=fixed_prefix_segments,
        alignment=alignment,
        blood_fraction_bins=blood_fraction_bins,
        iodine_concentration_bins_mg_mL=iodine_bins_mg_mL,
        min_fraction=min_fraction,
        min_radius_um=min_radius_um,
        supersample=supersample,
        prefix="xcat_basis",
    )

    (; generated, forest, results, raw, base_materials, alignment, fixed_prefix_segments, frame_dir, manifest)
else
    nothing
end

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0040
md"""
## Status

`run_pipeline = $(run_pipeline)`

When enabled, this notebook will:
1. regrow the validated XCAT coronary forest
2. simulate segment-level contrast transport
3. load the baseline XCAT raw + non-contrast material spreadsheet
4. estimate the raw-to-NURBS alignment
5. compile complete phantom frames for `t = $(frame_times_s[1]) : 1 : $(frame_times_s[end]) s`
"""

# ╔═╡ 9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0041
if workflow_state === nothing
    md"""
    Toggle `run_pipeline = true` when you want to write simulator-ready frame artifacts.
    The generated files will be placed under `$(output_dir)`.
    """
else
    md"""
    **Frame directory:** `$(workflow_state.frame_dir)`  
    **Run manifest:** `$(workflow_state.manifest)`  
    **Alignment origin (mm):** `$(workflow_state.alignment.origin_mm)`
    """
end

# ╔═╡ Cell order:
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0001
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0002
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0003
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0004
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0005
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0010
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0011
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0012
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0020
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0021
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0022
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0030
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0040
# ╠═9f7fd35a-b4ef-4d71-86d0-3a5a7a5d0041
