using Random
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
        verbose=true,
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
    forest = generated.forest
    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )
    return CoronaryForest(forest.trees, forest.territory_map, forest.params, pressure_params)
end

function nearest_segment(tree::VascularTree, point)
    best_i = 1
    best_d = Inf
    seg = tree.segments
    for i in 1:seg.n
        d = point_segment_distance(point[1], point[2], point[3], seg.proximal_x[i], seg.proximal_y[i], seg.proximal_z[i], seg.distal_x[i], seg.distal_y[i], seg.distal_z[i])
        if d < best_d
            best_d = d
            best_i = i
        end
    end
    return best_i, best_d
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

function first_positive_time(times, curve; threshold=1e-6)
    for i in eachindex(times)
        if curve[i] > threshold
            return times[i]
        end
    end
    return Inf
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    point = (136.46, 253.99, 290.35)
    times = collect(0.0:0.05:25.0)
    root_input = gamma_variate_input(times; amplitude=15.0, t0=0.0, tmax=4.0, alpha=8.0)

    forest = build_forest(nrb_path)
    tree = forest.trees["LCX"]
    params = forest.tree_params["LCX"]
    result = simulate_contrast_transport(tree, times; root_input=root_input, params=params, recompute_hemodynamics=true, return_outlet=true)
    work_tree = deepcopy(tree)
    compute_resistances!(work_tree, params.blood_viscosity)
    compute_flows!(work_tree, params)
    compute_pressures!(work_tree, params)

    seg_id, dist = nearest_segment(work_tree, point)
    seg = work_tree.segments
    conc = result.concentration[seg_id, :]
    println("nearest segment = ", seg_id)
    println("distance_mm = ", round(dist, digits=4))
    println("diameter_um = ", round(seg.radius[seg_id] * 2000.0, digits=4))
    println("length_mm = ", round(seg.seg_length[seg_id], digits=6))
    println("flow_m3_s = ", seg.flow[seg_id])
    println("local_tau_s = ", round(result.transit_time_s[seg_id], digits=4))
    println("path_tau_s = ", round(cumulative_path_transit(tree, result.transit_time_s, seg_id), digits=4))
    println("first_positive_1e-6 = ", first_positive_time(times, conc; threshold=1e-6))
    println("first_positive_1e-3 = ", first_positive_time(times, conc; threshold=1e-3))
    println("final_conc_mg_mL = ", conc[end])
    println("max_conc_mg_mL = ", maximum(conc))
end

main()
