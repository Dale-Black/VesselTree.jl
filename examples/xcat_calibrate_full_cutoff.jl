using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239

function env_int(name::String, default::Int)
    value = get(ENV, name, nothing)
    isnothing(value) && return default
    return parse(Int, value)
end

function env_bool(name::String, default::Bool)
    value = get(ENV, name, nothing)
    isnothing(value) && return default
    value = lowercase(strip(value))
    return value in ("1", "true", "yes", "on")
end

function actual_territory_fractions(result)
    pts = getfield(result.domain, :interior_points)
    tmap = result.forest.territory_map
    counts = Dict(name => 0 for name in tmap.tree_names)
    for i in axes(pts, 1)
        owner = query_territory(tmap, pts[i, 1], pts[i, 2], pts[i, 3])
        counts[owner] += 1
    end
    total = sum(values(counts))
    total == 0 && return Dict(name => 0.0 for name in keys(counts))
    return Dict(name => counts[name] / total for name in sort(collect(keys(counts))))
end

function proximal_root_stats(result)
    stats = Dict{String, NamedTuple}()
    for name in ["LAD", "LCX", "RCA"]
        xtree = result.xcat_trees[name]
        root = xtree.segments[xtree.root_segment]
        stats[name] = (
            first_radius_um=root.radii[1] * 1000.0,
            max_radius_um=maximum(root.radii) * 1000.0,
            imported_root_segment_radius_um=result.fixed_trees[name].segments.radius[Int(result.fixed_trees[name].root_segment_id)] * 1000.0,
        )
    end
    return stats
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB

    targets = Dict(
        "LAD" => env_int("LAD_TARGET", 12),
        "LCX" => env_int("LCX_TARGET", 8),
        "RCA" => env_int("RCA_TARGET", 10),
    )

    territory_fractions = Dict(
        "LAD" => parse(Float64, get(ENV, "LAD_TERRITORY", "0.40")),
        "LCX" => parse(Float64, get(ENV, "LCX_TERRITORY", "0.25")),
        "RCA" => parse(Float64, get(ENV, "RCA_TERRITORY", "0.35")),
    )

    rng = MersenneTwister(env_int("XCAT_SEED", 321))
    params = kassab_coronary_params()
    verbose = env_bool("XCAT_VERBOSE", true)
    target_interior = env_int("XCAT_TARGET_INTERIOR", 60_000)
    max_candidates = env_int("XCAT_MAX_CANDIDATES", 300_000)
    batch_size = env_int("XCAT_BATCH_SIZE", 8192)
    tree_capacity = env_int("XCAT_TREE_CAPACITY", 10_000)
    handoff_order = env_int("XCAT_HANDOFF_ORDER", 6)
    subdivision_max_order = env_int("XCAT_SUBDIVISION_MAX_ORDER", 11)
    max_tree_capacity = env_int("XCAT_MAX_TREE_CAPACITY", 2_000_000)
    apply_geometry = env_bool("XCAT_APPLY_GEOMETRY", false)
    enforce_full_cutoff = env_bool("XCAT_ENFORCE_FULL_CUTOFF", true)

    result = generate_xcat_kassab_coronary(
        nrb_path,
        params;
        phase="dias",
        rng=rng,
        verbose=verbose,
        handoff_order=handoff_order,
        subdivision_max_order=subdivision_max_order,
        apply_geometry=apply_geometry,
        enforce_full_cutoff=enforce_full_cutoff,
        max_tree_capacity=max_tree_capacity,
        target_interior=target_interior,
        max_candidates=max_candidates,
        batch_size=batch_size,
        tree_capacity=tree_capacity,
        target_terminals=targets,
        territory_fractions=territory_fractions,
    )

    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )

    println("\nCalibration summary")
    println("  targets = $(targets)")
    println("  territories = $(territory_fractions)")
    println("  actual territories = $(actual_territory_fractions(result))")
    println("  proximal root stats = $(proximal_root_stats(result))")
    println("  handoff_order = $(handoff_order), subdivision_max_order = $(subdivision_max_order)")
    println("  max_tree_capacity = $(max_tree_capacity), apply_geometry = $(apply_geometry)")
    hyperemic_targets = Dict("LAD" => 242.0, "LCX" => 116.0, "RCA" => 214.0)
    hyperemic_score = 0.0
    for name in ["LAD", "LCX", "RCA"]
        tree = result.forest.trees[name]
        baseline_tree = deepcopy(tree)
        hyperemic_tree = deepcopy(tree)
        apply_vasodilation!(hyperemic_tree; factor=1.6, min_diameter_um=0.0, max_diameter_um=400.0)

        baseline_flow = compute_root_flow_mLmin(baseline_tree, pressure_params[name])
        hyperemic_flow = compute_root_flow_mLmin(hyperemic_tree, pressure_params[name])
        hyperemic_score += abs(hyperemic_flow - hyperemic_targets[name])
        println(
            rpad(name, 6),
            "segments=", tree.segments.n,
            " terminals=", tree.n_terminals,
            " baseline=", round(baseline_flow, digits=3), " mL/min",
            " hyperemic=", round(hyperemic_flow, digits=3), " mL/min",
        )
    end
    println("  hyperemic absolute-error score = ", round(hyperemic_score, digits=3))
end

main()
