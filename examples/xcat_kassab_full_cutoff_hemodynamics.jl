using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB

    rng = MersenneTwister(321)
    params = kassab_coronary_params()
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

    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )

    println("Full-cutoff hemodynamic summary")
    for name in ["LAD", "LCX", "RCA"]
        baseline_tree = deepcopy(result.forest.trees[name])
        hyperemic_tree = deepcopy(result.forest.trees[name])
        apply_vasodilation!(hyperemic_tree; factor=1.6, min_diameter_um=0.0, max_diameter_um=400.0)

        baseline_flow = compute_root_flow_mLmin(baseline_tree, pressure_params[name])
        hyperemic_flow = compute_root_flow_mLmin(hyperemic_tree, pressure_params[name])

        println(
            rpad(name, 6),
            " baseline=", round(baseline_flow, digits=3), " mL/min",
            "   vasodilated=", round(hyperemic_flow, digits=3), " mL/min",
        )
    end
end

main()
