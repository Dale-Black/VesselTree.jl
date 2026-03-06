# basic_tree.jl — Generate a single coronary tree and validate against Kassab
#
# Usage: julia --project=. examples/basic_tree.jl

using VesselTree
using Random

println("=== VesselTree.jl: Basic Tree Example ===\n")

# Parameters
params = kassab_coronary_params()
rng = MersenneTwister(42)
domain = SphereDomain((0.0, 0.0, 0.0), 5.0)

# Generate a single tree with 100 terminals, handoff at order 5
configs = [
    TreeConfig("LAD", (0.0, 0.0, 0.0), 1.5, (0.0, -1.0, 0.0), 100, 1.0),
]

println("Generating tree (100 CCO terminals, handoff_order=5)...")
t0 = time()
forest = generate_kassab_coronary(
    domain, params;
    rng=rng, verbose=true, handoff_order=5,
    tree_configs=configs,
)
elapsed = time() - t0
println("\nCompleted in $(round(elapsed, digits=1))s\n")

tree = forest.trees["LAD"]

# Validation
println("--- Validation Report ---")
report = validate_tree(tree, params)
print_report(report)

# Report card
println("\n--- Kassab Report Card ---")
card = generate_report_card(tree, params)
print_report_card(stdout, card)

# Export
tmpdir = mktempdir()
vtp_path = export_centerlines_vtp(tree, joinpath(tmpdir, "basic_tree"))
println("\nVTP exported to: $vtp_path")

csv_path = export_csv(tree, joinpath(tmpdir, "basic_tree.csv"))
println("CSV exported to: $csv_path")

stl_path = export_stl(tree, joinpath(tmpdir, "basic_tree.stl"))
println("STL exported to: $stl_path")

json_path = export_graph_json(tree, joinpath(tmpdir, "basic_tree.json"))
println("JSON graph exported to: $json_path")

# Save/load roundtrip
jld2_path = joinpath(tmpdir, "basic_tree.jld2")
save_tree(jld2_path, tree)
loaded = load_tree(jld2_path)
println("\nJLD2 roundtrip: $(loaded.segments.n) segments (matches: $(loaded.segments.n == tree.segments.n))")

println("\n=== Done ===")
