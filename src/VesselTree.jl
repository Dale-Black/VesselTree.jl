module VesselTree

import AcceleratedKernels as AK
using LinearAlgebra
using Random
using Distributions
using StaticArrays
using JLD2
using WriteVTK
using DelimitedFiles

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
include("elements.jl")
include("barabasi.jl")
include("validation.jl")
include("subdivision.jl")
include("forest.jl")
include("contrast_transport.jl")
include("contrast_viewer.jl")
include("export.jl")
include("xcat_nrb.jl")
include("unified_viewer.jl")
include("xcat_integration.jl")
include("nrb_growth.jl")

export VascularTree, SegmentData, TreeTopology
export add_segment!, get_children, n_segments
export MorphometricParams, kassab_coronary_params, kassab_rca_params, kassab_lad_params, kassab_lcx_params, classify_order, with_hemodynamics
export AbstractDomain, SphereDomain, BoxDomain, EllipsoidDomain, EllipsoidShellDomain, CSVVolumeDomain, CSVShellDomain
export in_domain, sample_point, signed_distance, project_to_domain, default_coronary_domain, csv_volume_domain, csv_shell_domain, default_coronary_volume_domain, domain_bounds
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
export ContrastTransportResult, gamma_variate_input, root_pulse_input
export simulate_contrast_transport, simulate_forest_contrast, segment_mass_mg
export ContrastViewerTreeData, ContrastViewerData
export prepare_contrast_viewer_data, export_contrast_viewer_html, serve_static_directory
export UnifiedViewerPointLayerData, UnifiedViewerLineLayerData, UnifiedViewerData
export prepare_unified_viewer_data, export_unified_viewer_html
export assign_strahler_orders!, assign_strahler_orders_simple!, assign_diameter_defined_strahler!
export sample_asymmetry, compute_daughter_radii, apply_kassab_radii!, apply_full_kassab_radii!
export ElementData, group_into_elements, build_element_connectivity
export compute_element_statistics, compute_se_ratios
export build_empirical_connectivity, validate_connectivity, constrain_connectivity!
export compute_chi, compute_rho, classify_junction, compute_junction_angles, apply_junction_geometry!
export compute_trifurcation_angles, check_trifurcation_merge, merge_to_trifurcation!, apply_trifurcation_geometry!
export ValidationReport, validate_tree, print_report
export ks_test_onesample, validate_diameters_per_order, validate_lengths_per_order
export validate_segment_counts, validate_element_counts, validate_se_ratios
export validate_asymmetry_ks, validate_connectivity_matrix
export generate_report_card, print_report_card
export TreeConfig, TerritoryMap, initialize_territories, query_territory
export sample_in_territory, check_inter_tree_collision, coronary_tree_configs
export fixed_tree_configs, continue_coronary_forest
export estimate_total_segments, estimate_subdivision_capacity, subdivide_terminals!
export CoronaryForest, generate_coronary_forest, generate_kassab_coronary, validate_forest
export save_tree, load_tree, export_centerlines_vtp, export_forest_vtp
export export_stl, export_graph_json, export_csv
export export_wenbo_txt, export_forest_wenbo_txt
export XCATNurbsSurface, parse_xcat_nrb, xcat_bounds, xcat_center, xcat_object_dict, xcat_select_objects, xcat_summary_rows
export xcat_uv_counts, xcat_degrees, xcat_surface_point, xcat_surface_normal, xcat_sample_surface, xcat_export_sampled_surface_csv
export xcat_sampled_surface_rows, xcat_default_cavity_names, xcat_myocardial_shell_domain
export XCATCenterline, xcat_surface_axis, xcat_centerline_from_surface, xcat_resample_centerline, xcat_centerline_length_mm, xcat_centerline_summary_rows, xcat_export_centerline_csv
export xcat_reverse_centerline, xcat_merge_centerline_chain, xcat_build_coronary_trunks
export XCATTreeConnection, XCATCenterlineTree, xcat_build_coronary_trees, xcat_tree_summary_rows
export xcat_tree_to_vascular_tree, xcat_trees_to_vascular_trees
export xcat_default_coronary_surface_names, xcat_import_coronary_trees, generate_xcat_coronary_forest, generate_xcat_kassab_coronary
export NRBTreeSpec, NRBOrganSpec, xcat_heart_organ_spec
export nrb_vessel_surface_names, nrb_tree_target_terminals, nrb_tree_territory_fractions
export nrb_shell_domain, nrb_import_fixed_trees, generate_nrb_continuation_forest, generate_nrb_kassab_forest
export apply_vasodilation!, compute_root_flow_mLmin

# Visualization stubs (implemented by VesselTreeMakieExt when CairoMakie is loaded)
function plot_tree end
function plot_tree_2d end
function plot_validation_report end
function plot_forest end
export plot_tree, plot_tree_2d, plot_validation_report, plot_forest

end # module VesselTree
