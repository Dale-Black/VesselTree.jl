using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    rng = MersenneTwister(321)
    params = kassab_coronary_params()

    # Verified small-scale full-cutoff preset:
    # continuation stays moderate, but subdivision resolves all remaining orders
    # down to order 0 / cutoff scale.
    result = generate_xcat_kassab_coronary(
        nrb_path,
        params;
        phase="dias",
        rng=rng,
        verbose=true,
        handoff_order=6,
        subdivision_max_order=11,
        apply_geometry=false,
        enforce_full_cutoff=true,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        target_terminals=Dict(
            "LAD" => 12,
            "LCX" => 8,
            "RCA" => 10,
        ),
        territory_fractions=Dict(
            "LAD" => 0.40,
            "LCX" => 0.25,
            "RCA" => 0.35,
        ),
    )

    println("Full-cutoff summary")
    for name in ["LAD", "LCX", "RCA"]
        tree = result.forest.trees[name]
        tp = name == "LAD" ? kassab_lad_params() : name == "LCX" ? kassab_lcx_params() : kassab_rca_params()
        assign_strahler_orders!(tree, tp)
        terminal_orders = [Int(tree.topology.strahler_order[i]) for i in 1:tree.segments.n if tree.topology.is_terminal[i]]
        diam_um = [tree.segments.radius[i] * 2000.0 for i in 1:tree.segments.n if tree.topology.is_terminal[i]]
        println(
            rpad(name, 6),
            " segments=", tree.segments.n,
            " terminals=", tree.n_terminals,
            " max_terminal_order=", maximum(terminal_orders),
            " terminal_diam_um=[", round(minimum(diam_um), digits=2), ", ", round(maximum(diam_um), digits=2), "]",
        )
    end
end

main()
