using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239
const DEFAULT_ENFORCE_PHASE1_MORPHOMETRY = false
const DEFAULT_PHASE1_LENGTH_SIGMA = 5.0
const DEFAULT_TERRITORY_TREE_DISTANCE_WEIGHT = 0.0
const DEFAULT_TARGET_TERMINALS = Dict("LAD" => 18, "LCX" => 14, "RCA" => 18)
const DEFAULT_DT_S = 0.05
const DEFAULT_TEND_S = 25.0
const DEFAULT_INPUT_PEAK_CONCENTRATION_MG_ML = 15.0
const DEFAULT_INPUT_T0_S = 0.0
const DEFAULT_INPUT_TMAX_S = 4.0
const DEFAULT_INPUT_ALPHA = 8.0

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
        verbose=_env_bool("XCAT_SANITY_VERBOSE", true),
        handoff_order=6,
        subdivision_max_order=6,
        apply_geometry=false,
        enforce_phase1_morphometry=_env_bool("XCAT_ENFORCE_PHASE1_MORPHOMETRY", DEFAULT_ENFORCE_PHASE1_MORPHOMETRY),
        phase1_length_sigma=parse(Float64, get(ENV, "XCAT_PHASE1_LENGTH_SIGMA", string(DEFAULT_PHASE1_LENGTH_SIGMA))),
        phase1_morphometry_by_tree=_env_phase1_morphometry_by_tree(),
        phase1_length_sigma_by_tree=_env_phase1_length_sigma_by_tree(),
        territory_tree_distance_weight=parse(Float64, get(ENV, "XCAT_TERRITORY_TREE_DISTANCE_WEIGHT", string(DEFAULT_TERRITORY_TREE_DISTANCE_WEIGHT))),
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        target_terminals=_env_target_terminals(),
        territory_fractions=Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35),
    )

    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )

    forest = generated.forest
    return CoronaryForest(forest.trees, forest.territory_map, forest.params, pressure_params)
end

function cumulative_path_transit(tree::VascularTree, transit_time_s::AbstractVector, seg_id::Int)
    topo = tree.topology
    total = 0.0
    current = seg_id
    while current > 0
        total += transit_time_s[current]
        current = Int(topo.parent_id[current])
    end
    return total
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    dt = parse(Float64, get(ENV, "CONTRAST_SIMULATION_DT_S", string(DEFAULT_DT_S)))
    t_end = parse(Float64, get(ENV, "CONTRAST_SIMULATION_TEND_S", string(DEFAULT_TEND_S)))
    times = collect(0.0:dt:t_end)
    amplitude = parse(Float64, get(ENV, "CONTRAST_INPUT_PEAK_MG_ML", string(DEFAULT_INPUT_PEAK_CONCENTRATION_MG_ML)))
    t0 = parse(Float64, get(ENV, "CONTRAST_INPUT_T0_S", string(DEFAULT_INPUT_T0_S)))
    tmax = parse(Float64, get(ENV, "CONTRAST_INPUT_TMAX_S", string(DEFAULT_INPUT_TMAX_S)))
    alpha = parse(Float64, get(ENV, "CONTRAST_INPUT_ALPHA", string(DEFAULT_INPUT_ALPHA)))

    forest = build_forest(nrb_path)
    root_inputs = Dict(
        name => gamma_variate_input(times; amplitude=amplitude, t0=t0, tmax=tmax, alpha=alpha)
        for name in keys(forest.trees)
    )
    results = simulate_forest_contrast(forest, times; root_inputs=root_inputs, threaded=true)

    println("
Contrast sanity summary")
    println("  targets = ", _env_target_terminals())
    println("  phase1 morphometry gate = ", _env_bool("XCAT_ENFORCE_PHASE1_MORPHOMETRY", DEFAULT_ENFORCE_PHASE1_MORPHOMETRY))
    println("  phase1 sigma = ", parse(Float64, get(ENV, "XCAT_PHASE1_LENGTH_SIGMA", string(DEFAULT_PHASE1_LENGTH_SIGMA))))

    for name in sort(collect(keys(forest.trees)))
        tree = forest.trees[name]
        params = forest.tree_params[name]
        result = results[name]
        maxc = vec(maximum(result.concentration; dims=2))
        path_tau = [cumulative_path_transit(tree, result.transit_time_s, i) for i in 1:tree.segments.n]

        base_tree = deepcopy(tree)
        baseline_flow = compute_root_flow_mLmin(base_tree, params)
        hype_tree = deepcopy(tree)
        apply_vasodilation!(hype_tree; factor=1.6, min_diameter_um=0.0, max_diameter_um=400.0)
        hyperemic_flow = compute_root_flow_mLmin(hype_tree, params)

        println("
$(name)")
        println("  seg = ", tree.segments.n, ", terminals = ", tree.n_terminals)
        println("  baseline flow mL/min = ", round(baseline_flow, digits=3))
        println("  hyperemic flow mL/min = ", round(hyperemic_flow, digits=3))
        println("  median path_tau s = ", round(median(path_tau), digits=3))
        println("  p90 path_tau s = ", round(quantile(path_tau, 0.9), digits=3))
        println("  p99 path_tau s = ", round(quantile(path_tau, 0.99), digits=3))
        println("  max path_tau s = ", round(maximum(path_tau), digits=3))
        println("  never >1e-6 mg/mL = ", count(<=(1e-6), maxc))
        println("  never >1e-3 mg/mL = ", count(<=(1e-3), maxc))
        println("  median max concentration mg/mL = ", round(median(maxc), digits=3))
        println("  min positive max concentration mg/mL = ", round(minimum(filter(>(0.0), maxc); init=0.0), digits=6))
    end
end

main()
