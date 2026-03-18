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

_env_bool(name::AbstractString, default::Bool) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

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

    viewer = prepare_unified_viewer_data(
        forest,
        results;
        title="XCAT Coronary Unified Viewer",
        domain=generated.domain,
        original_surfaces=generated.vessel_surfaces,
        fixed_trees=generated.fixed_trees,
        domain_max_points=parse(Int, get(ENV, "UNIFIED_VIEWER_DOMAIN_MAX_POINTS", "12000")),
        cavity_max_points=parse(Int, get(ENV, "UNIFIED_VIEWER_CAVITY_MAX_POINTS", "6000")),
        surface_n_u=parse(Int, get(ENV, "UNIFIED_VIEWER_SURFACE_NU", "18")),
        surface_n_v=parse(Int, get(ENV, "UNIFIED_VIEWER_SURFACE_NV", "10")),
        original_surface_max_points=parse(Int, get(ENV, "UNIFIED_VIEWER_SURFACE_MAX_POINTS", "6000")),
        grown_max_segments_per_tree=parse(Int, get(ENV, "UNIFIED_VIEWER_GROWN_MAX_SEGMENTS", "3500")),
        iodine_min_radius_um=parse(Float64, get(ENV, "UNIFIED_VIEWER_MIN_RADIUS_UM", "4.0")),
        iodine_max_segments_per_tree=parse(Int, get(ENV, "UNIFIED_VIEWER_IODINE_MAX_SEGMENTS", "3000")),
        time_stride=parse(Int, get(ENV, "UNIFIED_VIEWER_TIME_STRIDE", "1")),
        radius_scale=parse(Float64, get(ENV, "UNIFIED_VIEWER_RADIUS_SCALE", "10.0")),
        blue_reference_mg_mL=parse(Float64, get(ENV, "UNIFIED_VIEWER_BLUE_REFERENCE", "3.0")),
    )

    output_dir = get(ENV, "UNIFIED_VIEWER_OUTPUT_DIR", joinpath(@__DIR__, "..", "output", "unified_viewer"))
    html_path = export_unified_viewer_html(joinpath(output_dir, "index.html"), viewer; title="XCAT Coronary Unified Viewer")
    println("\nUnified viewer HTML written to: ", html_path)

    if _env_bool("UNIFIED_VIEWER_SERVE", true)
        host = get(ENV, "UNIFIED_VIEWER_HOST", "127.0.0.1")
        port = parse(Int, get(ENV, "UNIFIED_VIEWER_PORT", "8008"))
        server = serve_static_directory(output_dir; host=host, port=port)
        println("Open this in your browser: http://", host, ":", port, "/")
        println("Use the layer toggles to switch domain / XCAT surfaces / fixed trunks / grown tree / iodine.")
        if _env_bool("UNIFIED_VIEWER_BLOCK", true)
            println("Server is running. Press Ctrl+C to stop.")
            while true
                sleep(3600.0)
            end
        else
            sleep(parse(Float64, get(ENV, "UNIFIED_VIEWER_SERVE_SECONDS", "2.0")))
            close(server)
        end
    end
end

main()
