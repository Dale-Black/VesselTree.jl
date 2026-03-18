using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239
const DEFAULT_TREE_NAME = "LAD"
const DEFAULT_PHASE = "dias"
const DEFAULT_TARGET_TERMINALS = 30
const DEFAULT_ENFORCE_PHASE1_MORPHOMETRY = false
const DEFAULT_PHASE1_LENGTH_SIGMA = 5.0

_env_bool(name::AbstractString, default::Bool) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

function single_tree_spec(name::String; phase::String=DEFAULT_PHASE, target_terminals::Int=DEFAULT_TARGET_TERMINALS)
    prefix = string(phase, "_")
    upper = uppercase(name)
    tree_spec = if upper == "LAD"
        NRBTreeSpec("LAD", [prefix * "lad1", prefix * "lad2", prefix * "lad3"]; target_terminals=target_terminals, territory_fraction=1.0)
    elseif upper == "LCX"
        NRBTreeSpec("LCX", [prefix * "lcx"]; target_terminals=target_terminals, territory_fraction=1.0)
    elseif upper == "RCA"
        NRBTreeSpec("RCA", [prefix * "rca1", prefix * "rca2"]; target_terminals=target_terminals, territory_fraction=1.0)
    else
        error("Unsupported tree name: $name")
    end
    return NRBOrganSpec(
        "XCAT $(upper) only",
        prefix * "pericardium",
        xcat_default_cavity_names(; phase=phase),
        [tree_spec],
        prefix * "aorta",
    )
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
    tree_name = get(ENV, "XCAT_SINGLE_TREE", DEFAULT_TREE_NAME)
    phase = get(ENV, "XCAT_PHASE", DEFAULT_PHASE)
    target_terminals = parse(Int, get(ENV, "XCAT_TARGET_TERMINALS", string(DEFAULT_TARGET_TERMINALS)))
    spec = single_tree_spec(tree_name; phase=phase, target_terminals=target_terminals)
    params = VesselTree._params_for_tree(tree_name, kassab_coronary_params())

    result = generate_nrb_kassab_forest(
        nrb_path,
        spec,
        params;
        rng=MersenneTwister(321),
        verbose=true,
        handoff_order=6,
        subdivision_max_order=6,
        apply_geometry=false,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        enforce_phase1_morphometry=_env_bool("XCAT_ENFORCE_PHASE1_MORPHOMETRY", DEFAULT_ENFORCE_PHASE1_MORPHOMETRY),
        phase1_length_sigma=parse(Float64, get(ENV, "XCAT_PHASE1_LENGTH_SIGMA", string(DEFAULT_PHASE1_LENGTH_SIGMA))),
    )

    pressure_params = with_hemodynamics(VesselTree._params_for_tree(tree_name, kassab_coronary_params()); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA)
    forest = CoronaryForest(result.forest.trees, result.forest.territory_map, result.forest.params, Dict(tree_name => pressure_params))
    tree = forest.trees[tree_name]

    times = collect(0.0:0.05:25.0)
    root_input = gamma_variate_input(times; amplitude=15.0, t0=0.0, tmax=4.0, alpha=8.0)
    contrast = simulate_contrast_transport(tree, times; root_input=root_input, params=pressure_params, recompute_hemodynamics=true, return_outlet=true)
    maxc = vec(maximum(contrast.concentration; dims=2))
    path_tau = [cumulative_path_transit(tree, contrast.transit_time_s, i) for i in 1:tree.segments.n]

    base_tree = deepcopy(tree)
    baseline_flow = compute_root_flow_mLmin(base_tree, pressure_params)
    hype_tree = deepcopy(tree)
    apply_vasodilation!(hype_tree; factor=1.6, min_diameter_um=0.0, max_diameter_um=400.0)
    hyperemic_flow = compute_root_flow_mLmin(hype_tree, pressure_params)

    println("
Single-tree sanity: ", tree_name)
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
end

main()
