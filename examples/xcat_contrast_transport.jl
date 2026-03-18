using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"
const MMHG_TO_PA = 133.32239

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    rng = MersenneTwister(321)

    forest = generate_xcat_kassab_coronary(
        nrb_path,
        kassab_coronary_params();
        phase="dias",
        rng=rng,
        verbose=true,
        handoff_order=6,
        subdivision_max_order=6,
        apply_geometry=false,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        target_terminals=Dict(
            "LAD" => 18,
            "LCX" => 14,
            "RCA" => 18,
        ),
        territory_fractions=Dict(
            "LAD" => 0.40,
            "LCX" => 0.25,
            "RCA" => 0.35,
        ),
    ).forest

    times = collect(0.0:0.05:25.0)
    root_inputs = Dict(
        "LAD" => gamma_variate_input(times; amplitude=15.0, t0=0.0, tmax=4.0, alpha=8.0),
        "LCX" => gamma_variate_input(times; amplitude=15.0, t0=0.0, tmax=4.0, alpha=8.0),
        "RCA" => gamma_variate_input(times; amplitude=15.0, t0=0.0, tmax=4.0, alpha=8.0),
    )

    pressure_params = Dict(
        "LAD" => with_hemodynamics(kassab_lad_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "LCX" => with_hemodynamics(kassab_lcx_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
        "RCA" => with_hemodynamics(kassab_rca_params(); root_pressure=100.0 * MMHG_TO_PA, terminal_pressure=15.0 * MMHG_TO_PA),
    )

    forest_with_params = CoronaryForest(forest.trees, forest.territory_map, forest.params, pressure_params)
    results = simulate_forest_contrast(
        forest_with_params,
        times;
        root_inputs=root_inputs,
        recompute_hemodynamics=true,
        storage_type=Float32,
        threaded=true,
    )

    println("\nContrast transport summary")
    for name in ["LAD", "LCX", "RCA"]
        conc = results[name].concentration
        peak_value = maximum(conc)
        peak_idx = argmax(conc)
        seg_idx, time_idx = Tuple(CartesianIndices(conc)[peak_idx])
        println(
            rpad(name, 6),
            "peak_conc=", round(peak_value, digits=3), " mg/mL",
            "  peak_time=", round(times[time_idx], digits=2), " s",
            "  peak_segment=", seg_idx,
        )
    end
end

main()
