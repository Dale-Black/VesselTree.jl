using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239

function length_stats_by_order(tree::VascularTree, order::Int)
    topo = tree.topology
    seg = tree.segments
    lengths = Float64[]
    for i in 1:seg.n
        if Int(topo.strahler_order[i]) == order
            push!(lengths, seg.seg_length[i])
        end
    end
    isempty(lengths) && return nothing
    return (n=length(lengths), min=minimum(lengths), med=median(lengths))
end

function summarize_result(label, result)
    println("\n=== ", label, " ===")
    for name in ["LAD", "LCX", "RCA"]
        tree = result.forest.trees[name]
        params = name == "LAD" ? kassab_lad_params() : name == "LCX" ? kassab_lcx_params() : kassab_rca_params()
        hparams = with_hemodynamics(params; root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA)
        base = deepcopy(tree)
        hype = deepcopy(tree)
        apply_vasodilation!(hype; factor=1.6, min_diameter_um=0.0, max_diameter_um=400.0)
        baseline = compute_root_flow_mLmin(base, hparams)
        hyperemic = compute_root_flow_mLmin(hype, hparams)
        println(name, ": seg=", tree.segments.n, ", term=", tree.n_terminals, ", baseline=", round(baseline, digits=3), ", hyperemic=", round(hyperemic, digits=3))
        if name == "LCX"
            for ord in (8, 9, 10)
                stats = length_stats_by_order(tree, ord)
                if stats !== nothing
                    println("  LCX order ", ord, ": n=", stats.n, ", min_mm=", round(stats.min, digits=4), ", med_mm=", round(stats.med, digits=4))
                end
            end
        end
    end
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    params = kassab_coronary_params()
    common = (
        phase="dias",
        rng=MersenneTwister(321),
        verbose=true,
        handoff_order=6,
        subdivision_max_order=6,
        apply_geometry=false,
        enforce_full_cutoff=false,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        max_tree_capacity=2_000_000,
        target_terminals=Dict("LAD" => 18, "LCX" => 14, "RCA" => 18),
        territory_fractions=Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35),
    )

    raw = generate_xcat_kassab_coronary(nrb_path, params; common..., centerline_resample=false)
    summarize_result("no resampling", raw)

    smooth = generate_xcat_kassab_coronary(
        nrb_path,
        params;
        common...,
        centerline_resample=true,
        centerline_min_spacing_mm=0.75,
        centerline_max_spacing_mm=1.5,
        centerline_radius_spacing_factor=0.75,
    )
    summarize_result("with resampling", smooth)
end

main()
