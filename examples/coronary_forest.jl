# coronary_forest.jl — Generate a complete LAD+LCX+RCA coronary forest
#
# Usage: julia --project=. examples/coronary_forest.jl

using VesselTree
using Random

println("=== VesselTree.jl: Coronary Forest Example ===\n")

params = kassab_coronary_params()
rng = MersenneTwister(123)
domain = SphereDomain((0.0, 0.0, 0.0), 8.0)

# Default coronary configs: LAD (40%), LCX (25%), RCA (35%)
# Using smaller terminal counts for demo speed
configs = [
    TreeConfig("LAD", (-2.0, 1.0, 0.0), 1.5, (0.0, -1.0, -0.5), 50, 0.40),
    TreeConfig("LCX", (-2.0, -1.0, 0.0), 1.2, (0.0, -1.0, 0.5), 30, 0.25),
    TreeConfig("RCA", (2.0, 0.0, 0.0), 1.3, (0.0, -1.0, 0.0), 40, 0.35),
]

println("Generating 3-tree coronary forest...")
println("  LAD: 50 terminals")
println("  LCX: 30 terminals")
println("  RCA: 40 terminals")
println()

t0 = time()
forest = generate_kassab_coronary(
    domain, params;
    rng=rng, verbose=true, handoff_order=4,
    tree_configs=configs,
)
elapsed = time() - t0
println("\nTotal generation time: $(round(elapsed, digits=1))s\n")

# Per-tree validation
for name in sort(collect(keys(forest.trees)))
    tree = forest.trees[name]
    println("--- $name: $(tree.segments.n) segments ---")
    report = validate_tree(tree, params)
    println("  Orders populated: $(length(report.diameter_ks_pvalues))")
    println("  Terminals: $(tree.n_terminals)")
    println("  Bifurcations: $(tree.n_bifurcations)")
    println("  Trifurcations: $(tree.n_trifurcations)")
end

# Export forest
tmpdir = mktempdir()
paths = export_forest_vtp(forest, tmpdir)
println("\nForest VTP files:")
for p in paths
    println("  $p")
end

println("\n=== Done ===")
