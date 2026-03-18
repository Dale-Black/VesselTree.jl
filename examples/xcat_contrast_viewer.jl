using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239

const DEFAULT_INPUT_PEAK_CONCENTRATION_MG_ML = 15.0
const DEFAULT_INPUT_T0_S = 0.0
const DEFAULT_INPUT_TMAX_S = 4.0
const DEFAULT_INPUT_ALPHA = 8.0
const DEFAULT_SIMULATION_DT_S = 0.05
const DEFAULT_SIMULATION_TEND_S = 25.0
const DEFAULT_CENTERLINE_RESAMPLE = true
const DEFAULT_CENTERLINE_MIN_SPACING_MM = 0.75
const DEFAULT_CENTERLINE_MAX_SPACING_MM = 1.5
const DEFAULT_CENTERLINE_RADIUS_SPACING_FACTOR = 0.75
const DEFAULT_ENFORCE_PHASE1_MORPHOMETRY = false
const DEFAULT_PHASE1_LENGTH_SIGMA = 5.0
const DEFAULT_TERRITORY_TREE_DISTANCE_WEIGHT = 0.0
const DEFAULT_TARGET_TERMINALS = Dict("LAD" => 18, "LCX" => 14, "RCA" => 18)

_env_bool(name::AbstractString, default::Bool) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

_env_target_terminals() = Dict(
    "LAD" => parse(Int, get(ENV, "XCAT_TARGET_TERMINALS_LAD", string(DEFAULT_TARGET_TERMINALS["LAD"]))),
    "LCX" => parse(Int, get(ENV, "XCAT_TARGET_TERMINALS_LCX", string(DEFAULT_TARGET_TERMINALS["LCX"]))),
    "RCA" => parse(Int, get(ENV, "XCAT_TARGET_TERMINALS_RCA", string(DEFAULT_TARGET_TERMINALS["RCA"]))),
)


_env_phase1_morphometry_by_tree() = begin
    config = Dict{String, Bool}()
    for name in ("LAD", "LCX", "RCA")
        key = "XCAT_ENFORCE_PHASE1_MORPHOMETRY_" * name
        haskey(ENV, key) || continue
        config[name] = _env_bool(key, false)
    end
    config
end

_env_phase1_length_sigma_by_tree() = begin
    config = Dict{String, Float64}()
    for name in ("LAD", "LCX", "RCA")
        key = "XCAT_PHASE1_LENGTH_SIGMA_" * name
        haskey(ENV, key) || continue
        config[name] = parse(Float64, ENV[key])
    end
    config
end

function build_forest(nrb_path::AbstractString)
    rng = MersenneTwister(321)
    generated = generate_xcat_kassab_coronary(
        nrb_path,
        kassab_coronary_params();
        phase="dias",
        rng=rng,
        verbose=true,
        handoff_order=6,
        subdivision_max_order=6,
        apply_geometry=false,
        enforce_phase1_morphometry=_env_bool("XCAT_ENFORCE_PHASE1_MORPHOMETRY", DEFAULT_ENFORCE_PHASE1_MORPHOMETRY),
        phase1_length_sigma=parse(Float64, get(ENV, "XCAT_PHASE1_LENGTH_SIGMA", string(DEFAULT_PHASE1_LENGTH_SIGMA))),
        phase1_morphometry_by_tree=_env_phase1_morphometry_by_tree(),
        phase1_length_sigma_by_tree=_env_phase1_length_sigma_by_tree(),
        territory_tree_distance_weight=parse(Float64, get(ENV, "XCAT_TERRITORY_TREE_DISTANCE_WEIGHT", string(DEFAULT_TERRITORY_TREE_DISTANCE_WEIGHT))),
        centerline_resample=_env_bool("XCAT_CENTERLINE_RESAMPLE", DEFAULT_CENTERLINE_RESAMPLE),
        centerline_min_spacing_mm=parse(Float64, get(ENV, "XCAT_CENTERLINE_MIN_SPACING_MM", string(DEFAULT_CENTERLINE_MIN_SPACING_MM))),
        centerline_max_spacing_mm=parse(Float64, get(ENV, "XCAT_CENTERLINE_MAX_SPACING_MM", string(DEFAULT_CENTERLINE_MAX_SPACING_MM))),
        centerline_radius_spacing_factor=parse(Float64, get(ENV, "XCAT_CENTERLINE_RADIUS_SPACING_FACTOR", string(DEFAULT_CENTERLINE_RADIUS_SPACING_FACTOR))),
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        target_terminals=_env_target_terminals(),
        territory_fractions=Dict(
            "LAD" => 0.40,
            "LCX" => 0.25,
            "RCA" => 0.35,
        ),
    )

    forest = generated.forest
    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )

    return CoronaryForest(forest.trees, forest.territory_map, forest.params, pressure_params)
end

function build_root_inputs(times)
    amplitude = parse(Float64, get(ENV, "CONTRAST_INPUT_PEAK_MG_ML", string(DEFAULT_INPUT_PEAK_CONCENTRATION_MG_ML)))
    t0 = parse(Float64, get(ENV, "CONTRAST_INPUT_T0_S", string(DEFAULT_INPUT_T0_S)))
    tmax = parse(Float64, get(ENV, "CONTRAST_INPUT_TMAX_S", string(DEFAULT_INPUT_TMAX_S)))
    alpha = parse(Float64, get(ENV, "CONTRAST_INPUT_ALPHA", string(DEFAULT_INPUT_ALPHA)))

    return Dict(
        "LAD" => gamma_variate_input(times; amplitude=amplitude, t0=t0, tmax=tmax, alpha=alpha),
        "LCX" => gamma_variate_input(times; amplitude=amplitude, t0=t0, tmax=tmax, alpha=alpha),
        "RCA" => gamma_variate_input(times; amplitude=amplitude, t0=t0, tmax=tmax, alpha=alpha),
    )
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    dt = parse(Float64, get(ENV, "CONTRAST_SIMULATION_DT_S", string(DEFAULT_SIMULATION_DT_S)))
    t_end = parse(Float64, get(ENV, "CONTRAST_SIMULATION_TEND_S", string(DEFAULT_SIMULATION_TEND_S)))
    times = collect(0.0:dt:t_end)

    forest = build_forest(nrb_path)
    results = simulate_forest_contrast(
        forest,
        times;
        root_inputs=build_root_inputs(times),
        recompute_hemodynamics=true,
        storage_type=Float32,
        threaded=true,
    )

    viewer = prepare_contrast_viewer_data(
        forest,
        results;
        min_radius_um=parse(Float64, get(ENV, "CONTRAST_VIEWER_MIN_RADIUS_UM", "4.0")),
        max_segments_per_tree=parse(Int, get(ENV, "CONTRAST_VIEWER_MAX_SEGMENTS", "2500")),
        time_stride=parse(Int, get(ENV, "CONTRAST_VIEWER_TIME_STRIDE", "1")),
        radius_scale=parse(Float64, get(ENV, "CONTRAST_VIEWER_RADIUS_SCALE", "10.0")),
        blue_reference_mg_mL=parse(Float64, get(ENV, "CONTRAST_VIEWER_BLUE_REFERENCE", "3.0")),
    )

    output_dir = get(ENV, "CONTRAST_VIEWER_OUTPUT_DIR", joinpath(@__DIR__, "..", "output", "contrast_viewer"))
    html_path = export_contrast_viewer_html(
        joinpath(output_dir, "index.html"),
        viewer;
        title="XCAT Coronary Contrast Viewer",
    )

    println("\nViewer HTML written to: ", html_path)

    if _env_bool("CONTRAST_VIEWER_SERVE", true)
        host = get(ENV, "CONTRAST_VIEWER_HOST", "127.0.0.1")
        port = parse(Int, get(ENV, "CONTRAST_VIEWER_PORT", "8008"))
        server = serve_static_directory(output_dir; host=host, port=port)
        println("Open this in your browser: http://", host, ":", port, "/")
        println("Use the slider or Play/Pause controls to inspect the iodine heat map over time.")

        if _env_bool("CONTRAST_VIEWER_BLOCK", true)
            println("Server is running. Press Ctrl+C to stop.")
            while true
                sleep(3600.0)
            end
        else
            sleep(parse(Float64, get(ENV, "CONTRAST_VIEWER_SERVE_SECONDS", "2.0")))
            close(server)
        end
    end
end

main()
