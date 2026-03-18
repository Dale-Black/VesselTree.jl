using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239
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
        verbose=false,
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

function ancestry(tree::VascularTree, seg_id::Int)
    topo = tree.topology
    chain = Int[]
    current = seg_id
    while current > 0
        push!(chain, current)
        current = Int(topo.parent_id[current])
    end
    reverse!(chain)
    return chain
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    times = collect(0.0:0.05:25.0)
    root_input = gamma_variate_input(times; amplitude=15.0, t0=0.0, tmax=4.0, alpha=8.0)

    forest = build_forest(nrb_path)
    tree = deepcopy(forest.trees["LAD"])
    params = forest.tree_params["LAD"]
    result = simulate_contrast_transport(tree, times; root_input=root_input, params=params, recompute_hemodynamics=true, return_outlet=true)

    compute_resistances!(tree, params.blood_viscosity)
    compute_flows!(tree, params)
    compute_pressures!(tree, params)

    seg = tree.segments
    topo = tree.topology
    maxc = vec(maximum(result.concentration; dims=2))
    path_tau = [cumulative_path_transit(tree, result.transit_time_s, i) for i in 1:seg.n]

    println("LAD summary")
    println("  seg=", seg.n)
    println("  median maxc=", median(maxc))
    println("  median path_tau=", median(path_tau))
    println("  p90 path_tau=", quantile(path_tau, 0.9))
    println("  p99 path_tau=", quantile(path_tau, 0.99))
    println("  never >1e-6 = ", count(<=(1e-6), maxc))
    println("  never >1e-3 = ", count(<=(1e-3), maxc))

    idx = sortperm(path_tau; rev=true)
    println("
Top 12 path_tau segments")
    for seg_id in idx[1:min(12, length(idx))]
        println(
            "seg=", seg_id,
            " order=", topo.strahler_order[seg_id],
            " terminal=", topo.is_terminal[seg_id],
            " path_tau=", round(path_tau[seg_id], digits=3),
            " local_tau=", round(result.transit_time_s[seg_id], digits=3),
            " len_mm=", round(seg.seg_length[seg_id], digits=3),
            " diam_um=", round(seg.radius[seg_id] * 2000.0, digits=3),
            " flow=", seg.flow[seg_id],
            " maxc=", maxc[seg_id],
        )
    end

    worst = idx[1]
    println("
Worst ancestry: ", worst)
    for seg_id in ancestry(tree, worst)
        println(
            "seg=", seg_id,
            " parent=", topo.parent_id[seg_id],
            " children=", get_children(topo, Int32(seg_id)),
            " order=", topo.strahler_order[seg_id],
            " len_mm=", round(seg.seg_length[seg_id], digits=3),
            " diam_um=", round(seg.radius[seg_id] * 2000.0, digits=3),
            " flow=", seg.flow[seg_id],
            " tau=", round(result.transit_time_s[seg_id], digits=3),
            " path_tau=", round(path_tau[seg_id], digits=3),
            " maxc=", maxc[seg_id],
        )
    end
end

main()
