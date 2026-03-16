using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239

function clone_tree(tree::VascularTree)
    return deepcopy(tree)
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(@__DIR__), "output", "xcat_wenbo")
    mkpath(output_dir)

    rng = MersenneTwister(123)
    params = kassab_coronary_params()
    result = generate_xcat_kassab_coronary(
        nrb_path,
        params;
        phase="dias",
        rng=rng,
        verbose=false,
        handoff_order=6,
        subdivision_max_order=6,
        target_interior=80_000,
        max_candidates=500_000,
        batch_size=8192,
        tree_capacity=20_000,
        target_terminals=Dict(
            "LAD" => 140,
            "LCX" => 90,
            "RCA" => 110,
        ),
        territory_fractions=Dict(
            "LAD" => 0.40,
            "LCX" => 0.25,
            "RCA" => 0.35,
        ),
    )

    export_paths = export_forest_wenbo_txt(result.forest, output_dir)
    println("Exported Wenbo txt files:")
    for path in export_paths
        println("  ", path)
    end

    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )

    println()
    println("Hemodynamic summary (built-in Poiseuille solver, 100 mmHg inlet, 15 mmHg outlet)")
    for name in ["LAD", "LCX", "RCA"]
        baseline_tree = clone_tree(result.forest.trees[name])
        baseline_flow = compute_root_flow_mLmin(baseline_tree, pressure_params[name])

        hyperemic_tree = clone_tree(result.forest.trees[name])
        apply_vasodilation!(hyperemic_tree; factor=1.6, min_diameter_um=0.0, max_diameter_um=400.0)
        hyperemic_flow = compute_root_flow_mLmin(hyperemic_tree, pressure_params[name])

        println(
            rpad(name, 6),
            " baseline=", round(baseline_flow, digits=3), " mL/min",
            "   vasodilated=", round(hyperemic_flow, digits=3), " mL/min",
        )
    end
end

main()
