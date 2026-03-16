using Random
using VesselTree

# Template for organs beyond the heart.
# Replace the placeholder surface names with the names present in your .nrb file.

const DEFAULT_NRB = "/path/to/your_model.nrb"

function example_organ_spec()
    return NRBOrganSpec(
        "Example Organ",
        "example_outer_surface",
        ["example_cavity_1", "example_cavity_2"],
        [
            NRBTreeSpec(
                "MAIN_TREE_A",
                ["example_tree_a_root", "example_tree_a_mid", "example_tree_a_distal"];
                target_terminals=40,
                territory_fraction=0.55,
            ),
            NRBTreeSpec(
                "MAIN_TREE_B",
                ["example_tree_b_root", "example_tree_b_distal"];
                target_terminals=30,
                territory_fraction=0.45,
            ),
        ],
        "example_reference_surface",
    )
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_NRB
    spec = example_organ_spec()
    rng = MersenneTwister(42)
    params = kassab_coronary_params()

    result = generate_nrb_continuation_forest(
        nrb_path,
        spec,
        params;
        rng=rng,
        verbose=true,
        target_interior=40_000,
        max_candidates=200_000,
        batch_size=4096,
        tree_capacity=5000,
    )

    println("Generated continuation forest for $(spec.name)")
    for (name, tree) in sort(collect(result.forest.trees), by=x -> x[1])
        println("  $(name): $(tree.segments.n) segments, $(tree.n_terminals) terminals")
    end
end

main()
