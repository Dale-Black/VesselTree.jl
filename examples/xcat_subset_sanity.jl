using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239
const DEFAULT_PHASE = "dias"
const DEFAULT_TREE_SET = ["LAD", "LCX", "RCA"]
const DEFAULT_TARGET_TERMINALS = Dict("LAD" => 18, "LCX" => 14, "RCA" => 18)
const DEFAULT_TERRITORY_FRACTIONS = Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35)

_env_bool(name::AbstractString, default::Bool) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

function _tree_spec(name::String, phase::String, target_terminals::Int, territory_fraction::Float64)
    prefix = string(phase, "_")
    upper = uppercase(name)
    if upper == "LAD"
        return NRBTreeSpec("LAD", [prefix * "lad1", prefix * "lad2", prefix * "lad3"]; target_terminals=target_terminals, territory_fraction=territory_fraction)
    elseif upper == "LCX"
        return NRBTreeSpec("LCX", [prefix * "lcx"]; target_terminals=target_terminals, territory_fraction=territory_fraction)
    elseif upper == "RCA"
        return NRBTreeSpec("RCA", [prefix * "rca1", prefix * "rca2"]; target_terminals=target_terminals, territory_fraction=territory_fraction)
    else
        error("Unsupported tree name: $name")
    end
end

function subset_spec(names::Vector{String}; phase::String=DEFAULT_PHASE, use_original_fractions::Bool=false)
    fractions = if use_original_fractions
        raw = [DEFAULT_TERRITORY_FRACTIONS[name] for name in names]
        total = sum(raw)
        [v / total for v in raw]
    else
        fill(1.0 / length(names), length(names))
    end
    tree_specs = [_tree_spec(name, phase, DEFAULT_TARGET_TERMINALS[name], fractions[i]) for (i, name) in enumerate(names)]
    return NRBOrganSpec(
        "XCAT subset",
        string(phase, "_pericardium"),
        xcat_default_cavity_names(; phase=phase),
        tree_specs,
        string(phase, "_aorta"),
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
    names = split(get(ENV, "XCAT_TREE_SET", join(DEFAULT_TREE_SET, ",")), ',')
    names = [uppercase(strip(name)) for name in names if !isempty(strip(name))]
    spec = subset_spec(names; use_original_fractions=_env_bool("XCAT_SUBSET_USE_ORIGINAL_FRACTIONS", false))

    tree_params = Dict(name => with_hemodynamics(VesselTree._params_for_tree(name, kassab_coronary_params()); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA) for name in names)
    result = generate_nrb_kassab_forest(
        nrb_path,
        spec,
        kassab_coronary_params();
        rng=MersenneTwister(321),
        verbose=true,
        handoff_order=6,
        subdivision_max_order=6,
        apply_geometry=false,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        tree_params=tree_params,
    )
    forest = CoronaryForest(result.forest.trees, result.forest.territory_map, result.forest.params, tree_params)

    times = collect(0.0:0.05:25.0)
    root_inputs = Dict(name => gamma_variate_input(times; amplitude=15.0, t0=0.0, tmax=4.0, alpha=8.0) for name in names)
    contrast = simulate_forest_contrast(forest, times; root_inputs=root_inputs, threaded=true)

    println("
Subset sanity: ", names)
    for name in sort(names)
        tree = forest.trees[name]
        params = forest.tree_params[name]
        result = contrast[name]
        maxc = vec(maximum(result.concentration; dims=2))
        path_tau = [cumulative_path_transit(tree, result.transit_time_s, i) for i in 1:tree.segments.n]
        println("
$(name)")
        println("  seg = ", tree.segments.n, ", terminals = ", tree.n_terminals)
        println("  median path_tau s = ", round(median(path_tau), digits=3))
        println("  p99 path_tau s = ", round(quantile(path_tau, 0.99), digits=3))
        println("  max path_tau s = ", round(maximum(path_tau), digits=3))
        println("  never >1e-6 mg/mL = ", count(<=(1e-6), maxc))
        println("  never >1e-3 mg/mL = ", count(<=(1e-3), maxc))
        println("  baseline flow mL/min = ", round(compute_root_flow_mLmin(deepcopy(tree), params), digits=3))
    end
end

main()
