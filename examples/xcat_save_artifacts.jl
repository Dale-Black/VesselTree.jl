using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239
const DEFAULT_TARGET_TERMINALS = Dict("LAD" => 8, "LCX" => 72, "RCA" => 24)
const DEFAULT_PREFIX_BLEND = Dict("LAD" => 0.0, "LCX" => 0.5, "RCA" => 0.3)
const DEFAULT_ROOT_TARGET_DIAM_UM = Dict("LAD" => 3500.0, "LCX" => 4500.0, "RCA" => 4500.0)

const DEFAULT_INPUT_PEAK_CONCENTRATION_MG_ML = 15.0
const DEFAULT_INPUT_T0_S = 0.0
const DEFAULT_INPUT_TMAX_S = 4.0
const DEFAULT_INPUT_ALPHA = 8.0
const DEFAULT_SIMULATION_DT_S = 0.05
const DEFAULT_SIMULATION_TEND_S = 25.0
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "saved_runs")

_env_target_terminals() = Dict(
    "LAD" => parse(Int, get(ENV, "XCAT_TARGET_TERMINALS_LAD", string(DEFAULT_TARGET_TERMINALS["LAD"]))),
    "LCX" => parse(Int, get(ENV, "XCAT_TARGET_TERMINALS_LCX", string(DEFAULT_TARGET_TERMINALS["LCX"]))),
    "RCA" => parse(Int, get(ENV, "XCAT_TARGET_TERMINALS_RCA", string(DEFAULT_TARGET_TERMINALS["RCA"]))),
)

_env_prefix_blend() = Dict(
    "LAD" => parse(Float64, get(ENV, "XCAT_PREFIX_BLEND_LAD", string(DEFAULT_PREFIX_BLEND["LAD"]))),
    "LCX" => parse(Float64, get(ENV, "XCAT_PREFIX_BLEND_LCX", string(DEFAULT_PREFIX_BLEND["LCX"]))),
    "RCA" => parse(Float64, get(ENV, "XCAT_PREFIX_BLEND_RCA", string(DEFAULT_PREFIX_BLEND["RCA"]))),
)

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
        target_terminals=_env_target_terminals(),
        territory_fractions=Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35),
    )
end

function build_hemodynamic_forest(generated)
    blend = _env_prefix_blend()
    trees = Dict{String, VascularTree}()
    base_params = kassab_coronary_params()
    for name in ["LAD", "LCX", "RCA"]
        tree = deepcopy(generated.forest.trees[name])
        update_radii!(tree, base_params.gamma)
        src = generated.fixed_trees[name]
        target_diam_um = DEFAULT_ROOT_TARGET_DIAM_UM[name]
        scale = target_diam_um / (src.segments.radius[Int(src.root_segment_id)] * 2000.0)
        alpha = blend[name]
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

function build_root_inputs(times)
    amplitude = parse(Float64, get(ENV, "CONTRAST_INPUT_PEAK_MG_ML", string(DEFAULT_INPUT_PEAK_CONCENTRATION_MG_ML)))
    t0 = parse(Float64, get(ENV, "CONTRAST_INPUT_T0_S", string(DEFAULT_INPUT_T0_S)))
    tmax = parse(Float64, get(ENV, "CONTRAST_INPUT_TMAX_S", string(DEFAULT_INPUT_TMAX_S)))
    alpha = parse(Float64, get(ENV, "CONTRAST_INPUT_ALPHA", string(DEFAULT_INPUT_ALPHA)))
    return Dict(name => gamma_variate_input(times; amplitude=amplitude, t0=t0, tmax=tmax, alpha=alpha) for name in ("LAD", "LCX", "RCA"))
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    dt = parse(Float64, get(ENV, "CONTRAST_SIMULATION_DT_S", string(DEFAULT_SIMULATION_DT_S)))
    t_end = parse(Float64, get(ENV, "CONTRAST_SIMULATION_TEND_S", string(DEFAULT_SIMULATION_TEND_S)))
    times = collect(0.0:dt:t_end)

    generated = build_generated(nrb_path)
    forest = build_hemodynamic_forest(generated)
    results = simulate_forest_contrast(
        forest,
        times;
        root_inputs=build_root_inputs(times),
        recompute_hemodynamics=true,
        storage_type=Float32,
        threaded=true,
    )

    output_dir = get(ENV, "XCAT_ARTIFACT_OUTPUT_DIR", DEFAULT_OUTPUT_DIR)
    dense_mode = Symbol(lowercase(get(ENV, "XCAT_CONTRAST_DENSE_MODE", "auto")))
    voxel_spacing = parse(Float64, get(ENV, "XCAT_CONTRAST_VOXEL_SPACING_MM", "1.0"))
    voxel_supersample = parse(Int, get(ENV, "XCAT_CONTRAST_VOXEL_SUPERSAMPLE", "2"))

    bundle = save_xcat_run_artifacts(
        output_dir,
        nrb_path,
        forest,
        results;
        voxel_spacing_mm=(voxel_spacing, voxel_spacing, voxel_spacing),
        voxel_padding_mm=(2.0, 2.0, 2.0),
        voxel_supersample=voxel_supersample,
        dense_mode=dense_mode,
    )

    println("
Saved XCAT run artifacts:")
    println("  run dir:      ", bundle.run_dir)
    println("  trees:        ", bundle.tree_snapshot.manifest_path)
    println("  contrast:     ", bundle.contrast_snapshot.manifest_path)
    println("  fused nrb:    ", bundle.fused_nrb_path)
    println("  run manifest: ", bundle.manifest_path)
end

main()
