# custom_vasculature.jl — Create custom parameters for non-coronary vasculature
#
# Usage: julia --project=. examples/custom_vasculature.jl

using VesselTree
using Random

println("=== VesselTree.jl: Custom Vasculature Example ===\n")

# Start with coronary params as a template
base = kassab_coronary_params()

# Create cerebral-like parameters with fewer orders and smaller vessels
# MorphometricParams fields (positional):
#   gamma, diameter_mean, diameter_sd, length_mean, length_sd,
#   diameter_bounds, connectivity_matrix,
#   asymmetry_alpha, asymmetry_beta,
#   trifurcation_chi_th, sprouting_rho_th,
#   blood_viscosity, root_pressure, terminal_pressure,
#   vessel_cutoff_um, n_orders

n_ord = 8
diam_mean = [8.0, 15.0, 30.0, 60.0, 120.0, 250.0, 500.0, 1000.0]
diam_sd = [1.0, 3.0, 6.0, 12.0, 25.0, 50.0, 100.0, 200.0]

# Compute diameter bounds from means (geometric midpoints)
diam_bounds = Float64[0.0]
for i in 1:(n_ord - 1)
    push!(diam_bounds, sqrt(diam_mean[i] * diam_mean[i + 1]))
end
push!(diam_bounds, diam_mean[end] * 2.0)

custom_params = MorphometricParams(
    base.gamma,
    diam_mean,
    diam_sd,
    [20.0, 40.0, 100.0, 250.0, 600.0, 1500.0, 4000.0, 10000.0],
    [5.0, 10.0, 25.0, 60.0, 150.0, 400.0, 1000.0, 2500.0],
    diam_bounds,
    base.connectivity_matrix[1:n_ord, 1:n_ord],
    base.asymmetry_alpha,
    base.asymmetry_beta,
    base.trifurcation_chi_th,
    base.sprouting_rho_th,
    base.blood_viscosity,
    base.root_pressure,
    base.terminal_pressure,
    base.vessel_cutoff_um,
    n_ord,
)

rng = MersenneTwister(77)
domain = SphereDomain((0.0, 0.0, 0.0), 4.0)

configs = [
    TreeConfig("MCA", (0.0, 0.0, 0.0), 0.5, (1.0, -1.0, 0.0), 30, 1.0),
]

println("Generating cerebral-like tree (8 orders, 30 terminals)...")
t0 = time()
forest = generate_kassab_coronary(
    domain, custom_params;
    rng=rng, verbose=true, handoff_order=4,
    tree_configs=configs,
)
elapsed = time() - t0

tree = forest.trees["MCA"]
println("\nGenerated $(tree.segments.n) segments in $(round(elapsed, digits=1))s")

# Show diameter range
seg = tree.segments
n = tree.segments.n
min_d = minimum(seg.radius[i] * 2 * 1000 for i in 1:n)
max_d = maximum(seg.radius[i] * 2 * 1000 for i in 1:n)
println("Diameter range: $(round(min_d, digits=1)) - $(round(max_d, digits=1)) um")

# Export
tmpdir = mktempdir()
vtp_path = export_centerlines_vtp(tree, joinpath(tmpdir, "cerebral"))
println("VTP exported to: $vtp_path")

println("\n=== Done ===")
