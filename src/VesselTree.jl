module VesselTree

import AcceleratedKernels as AK
using LinearAlgebra
using Random
using Distributions
using StaticArrays

# Scientific constants
const MURRAY_GAMMA = 7 / 3           # Huo-Kassab 2007 (NOT 3.0)
const VESSEL_CUTOFF_UM = 8.0         # Capillary diameter (um)
const BLOOD_VISCOSITY = 0.0035       # Pa*s (3.5 cP)
const ROOT_PRESSURE = 13332.0        # Pa (100 mmHg)
const TERMINAL_PRESSURE = 3999.6     # Pa (30 mmHg)
const TRIFURCATION_CHI_TH = 0.83     # Barabasi 2026
const SPROUTING_RHO_TH = 0.83        # Barabasi 2026
const ASYMMETRY_ALPHA = 2.5          # Beta distribution parameter
const ASYMMETRY_BETA = 0.8           # Beta distribution parameter

include("types.jl")
include("parameters.jl")
include("domain.jl")
include("distance.jl")
include("intersection.jl")
include("kamiya.jl")
include("murray.jl")
include("spatial.jl")
include("growth.jl")
include("hemodynamics.jl")
include("kassab.jl")
include("barabasi.jl")

export VascularTree, SegmentData, TreeTopology
export add_segment!, get_children, n_segments
export MorphometricParams, kassab_coronary_params, classify_order
export AbstractDomain, SphereDomain, BoxDomain, EllipsoidDomain
export in_domain, sample_point, signed_distance
export point_segment_distance, point_segment_distance_sq
export compute_all_distances!, find_nearest_segments
export segments_intersect, segments_min_distance_sq
export check_intersections!, has_any_intersection, check_domain_crossing
export surface_cost_at_t, optimize_bifurcation_point
export evaluate_all_connections!, select_best_connection
export compute_radii_symmetric, compute_radii_asymmetric
export update_radii!, grow_tree!, add_bifurcation!
export SpatialGrid, build_grid, query_nearby
export compute_resistances!, compute_flows!, compute_pressures!, validate_hemodynamics
export assign_terminal_flows!, recompute_radii_from_flow!
export assign_strahler_orders!, sample_asymmetry, compute_daughter_radii, apply_kassab_radii!
export build_empirical_connectivity, validate_connectivity, constrain_connectivity!
export compute_chi, compute_rho, classify_junction, compute_junction_angles, apply_junction_geometry!
export compute_trifurcation_angles, check_trifurcation_merge, merge_to_trifurcation!, apply_trifurcation_geometry!

end # module VesselTree
