using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const DEFAULT_PHASE = "dias"
const DEFAULT_TARGET_TERMINALS = Dict("LAD" => 18, "LCX" => 14, "RCA" => 18)
const DEFAULT_TERRITORY_FRACTIONS = Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35)

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    phase = get(ENV, "XCAT_PHASE", DEFAULT_PHASE)
    names = split(get(ENV, "XCAT_TREE_SET", "LAD,LCX,RCA"), ',')
    names = [uppercase(strip(name)) for name in names if !isempty(strip(name))]
    use_original = lowercase(get(ENV, "XCAT_SUBSET_USE_ORIGINAL_FRACTIONS", "true")) in ("1","true","yes","on")

    fractions = if use_original
        raw = [DEFAULT_TERRITORY_FRACTIONS[name] for name in names]
        total = sum(raw)
        Dict(name => raw[i] / total for (i, name) in enumerate(names))
    else
        Dict(name => 1.0 / length(names) for name in names)
    end
    targets = Dict(name => parse(Int, get(ENV, "XCAT_TARGET_TERMINALS_" * name, string(DEFAULT_TARGET_TERMINALS[name]))) for name in names)

    spec = let prefix = string(phase, "_")
        tree_specs = NRBTreeSpec[]
        for name in names
            if name == "LAD"
                push!(tree_specs, NRBTreeSpec("LAD", [prefix * "lad1", prefix * "lad2", prefix * "lad3"]; target_terminals=targets[name], territory_fraction=fractions[name]))
            elseif name == "LCX"
                push!(tree_specs, NRBTreeSpec("LCX", [prefix * "lcx"]; target_terminals=targets[name], territory_fraction=fractions[name]))
            elseif name == "RCA"
                push!(tree_specs, NRBTreeSpec("RCA", [prefix * "rca1", prefix * "rca2"]; target_terminals=targets[name], territory_fraction=fractions[name]))
            end
        end
        NRBOrganSpec("phase1 diagnostic", prefix * "pericardium", xcat_default_cavity_names(; phase=phase), tree_specs, prefix * "aorta")
    end

    tree_params = Dict(name => VesselTree._params_for_tree(name, kassab_coronary_params()) for name in names)
    result = generate_nrb_continuation_forest(
        nrb_path,
        spec,
        kassab_coronary_params();
        rng=MersenneTwister(321),
        verbose=true,
        kassab=true,
        handoff_order=6,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        tree_params=tree_params,
    )

    println("
Phase-1 terminal order summary")
    for name in sort(names)
        tree = result.forest.trees[name]
        params = tree_params[name]
        assign_strahler_orders!(tree, params)
        topo = tree.topology
        orders = Int[]
        for i in 1:tree.segments.n
            topo.is_terminal[i] || continue
            push!(orders, Int(topo.strahler_order[i]))
        end
        counts = Dict{Int, Int}()
        for ord in orders
            counts[ord] = get(counts, ord, 0) + 1
        end
        println("
$(name)")
        println("  terminals = ", length(orders))
        println("  median order = ", median(orders))
        println("  min/max order = ", minimum(orders), " / ", maximum(orders))
        println("  histogram = ", sort(collect(counts); by=first))
    end
end

main()
