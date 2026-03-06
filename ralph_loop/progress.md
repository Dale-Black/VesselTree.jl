# VesselTree.jl — Progress Log

## Baseline (2026-03-05)

### Starting State
- Source files: 0 (greenfield project)
- Test count: 0 (fresh start)
- Dependencies: AcceleratedKernels.jl (mandatory), StaticArrays, Distributions, NearestNeighbors, WriteVTK, JLD2

### Regression Gate Baseline
- Julia tests: 0 pass (fresh)
- AK compliance: N/A (no source yet)
- Package loads: untested

### Story Status
- Total: 28 stories (VESSEL-1000 to VESSEL-1027)
- Done: 0
- Open: 28
- Milestones: P0 (Discovery) through P8 (Integration)

### Key Scientific Parameters
- Murray exponent: 7/3 = 2.33 (Huo-Kassab 2007, NOT 3.0)
- Vessel cutoff: 8 um (capillary, Kassab Order 0)
- Trifurcation threshold chi: 0.83 (Barabasi 2026)
- Sprouting threshold rho: 0.83 (blood vessels)
- Asymmetry median: 0.76 (Kassab 1993, Beta(2.5, 0.8))

## Log

### 2026-03-05: VESSEL-1000 [PASS]
- Attempted: Deep-dive into AcceleratedKernels.jl v0.4.3 API for vascular tree generation patterns
- Result: All 10 research questions answered with locally verified code snippets
- Regression gate: N/A (discovery story, no source code)
- Learning:
  - No `argmin`/`argmax` in AK — implement via `mapreduce` with `(value, index)` tuples and explicit `neutral`
  - `merge_sort_by_key!` is GPU-only — use `sortperm` on CPU
  - CPU overhead is ~20-70us per AK call; crossover at ~100K elements with 4 threads
  - At 2M segments (our target): 1.83x speedup with 4 CPU threads; GPU would be dramatically better
  - SoA layout with `@view field[1:n]` works perfectly for growing arrays
  - Closures capture external arrays correctly but must be in function scope (not module scope) for GPU
  - Errors propagate as `TaskFailedException` on CPU; guard domain errors (sqrt, log) explicitly
  - `min_elems` parameter is critical for tuning small-array performance
  - `AK.any(identity, bool_array)` is the pattern for intersection detection
  - Pre-allocate index arrays for mapreduce-based argmin to avoid GC pressure
- Next: VESSEL-1001 (OpenCCO algorithm deep-dive) or other P0 discovery stories

### 2026-03-05: VESSEL-1001 [PASS]
- Attempted: Deep-dive into OpenCCO C++ codebase and IPOL paper for CCO algorithm
- Result: Complete algorithm documented with all equations, pseudocode, data structures, and scaling analysis
- Regression gate: N/A (discovery story, no source code)
- Learning:
  - CCO uses volume minimization; we'll use surface minimization (Barabasi)
  - OpenCCO defaults to gamma=3.0; we use 7/3 (Huo-Kassab)
  - Kamiya optimization solves a 2x2 nonlinear system — can use simple Newton instead of Ceres
  - Key bottleneck is O(k^2) scaling from brute-force intersection checking — spatial indexing essential
  - OpenCCO stores only distal points (proximal = parent's distal); we store both for AK parallelism
  - Deepcopy for each candidate is expensive; use undo/rollback instead
  - Distance threshold shrinks by 0.9x if no valid point found after max_trials
  - Degenerate check: 2*r <= l (diameter must not exceed segment length)
  - OpenCCO demonstrates up to 20K terminals; our target is 2M (100x more)
  - Length factor approach: scale distances rather than moving points (avoids O(k) coordinate updates)
- Next: VESSEL-1002 (Kassab morphometry) or VESSEL-1003 (Barabasi surface optimization)

### 2026-03-05: VESSEL-1002 [PASS]
- Attempted: Extract and document all quantitative data from Kassab 1993 and related papers
- Result: Complete morphometric reference compiled from Kassab 1993, Wischgoll 2009, Kaimovitz 2005/2010, Huo-Kassab 2007/2009/2012
- Regression gate: N/A (discovery story, no source code)
- Learning:
  - 12 orders (0-11): capillaries (8um) to main stem (4500um)
  - Diameter boundaries follow ~1.8x geometric progression
  - Asymmetry ratio well-modeled by Beta(2.5, 0.8), median=0.76
  - ~6M segments per major coronary artery; ~50% are terminals
  - L/D ratio increases with order: 2.5 (capillary) to 33.3 (stem)
  - Gamma = 7/3 from tree-level energy minimization (NOT single-bifurcation Murray)
  - Angle distribution is bimodal (sprouting ~90 deg + branching ~45 deg)
  - Connectivity matrix shows highly asymmetric branching at higher orders
  - Human data more symmetric than porcine (median S=0.59 vs <0.40 for large vessels)
  - Validation: KS test per order, Murray's law check, angle bimodality
- Next: VESSEL-1003 (Barabasi surface optimization)

### 2026-03-05: VESSEL-1003 [PASS]
- Attempted: Study Barabasi 2026 Nature paper and min-surf-netw code for surface optimization
- Result: Complete documentation of chi/rho parameters, sprouting/branching regimes, trifurcation stability, junction surface computation, and implementation strategy
- Regression gate: N/A (discovery story, no source code)
- Learning:
  - chi = w/r (circumference/spacing); trifurcation transition at chi=0.83
  - rho = r_small/r_large; sprouting (perpendicular) when rho < 0.83, branching (Y-shape) when rho >= 0.83
  - Blood vessels: 8.3% trifurcations, 25.6% sprouting, 91% bifurcations
  - min-surf-netw is Mathematica (NOT Python as assumed in PRD) using FEM mesh optimization
  - Full PDE solver is overkill for our use — use threshold rules + Murray angles instead
  - Surface cost (r*l) vs volume cost (r^2*l) — gamma=7/3 is consistent with surface dominance
  - Bimodal angle distribution emerges naturally from sprouting+branching separation
  - Lambda order parameter: lambda=l/w, lambda->0 at trifurcation merge
  - Junction surfaces modeled via "pants decomposition" of tube cross-sections
  - All real networks are ~25% longer than Steiner-optimal (surface cost penalty)
- Next: VESSEL-1004 (Project setup: Package skeleton with AK.jl) — all P0 discovery complete

### 2026-03-05: VESSEL-1004 [PASS]
- Attempted: Create Julia package skeleton with AcceleratedKernels.jl and scientific constants
- Result: Package loads, 21 tests pass, CLAUDE.md created
- Regression gate: 21 tests pass, 0 fail, 0 error
- Files created:
  - Project.toml — deps: AcceleratedKernels 0.4, Distributions 0.25, StaticArrays 1, LinearAlgebra, Random
  - src/VesselTree.jl — module with 9 scientific constants (MURRAY_GAMMA=7/3, VESSEL_CUTOFF_UM=8.0, etc.)
  - test/runtests.jl — 3 testsets: module loading (3), scientific constants (12), AK availability (6)
  - CLAUDE.md — project context, conventions, AK gotchas, scientific constants
- Learning:
  - AK does not support BitVector — must convert to Vector{Bool} before passing to AK.any/all
  - Correct AcceleratedKernels UUID is 6a4ca0a5-0e36-4168-a932-d9be78d558f1 (not the one in General registry search)
  - `import AcceleratedKernels as AK` pattern works well for namespacing
- Next: VESSEL-1005 (Core types: SoA segment storage + tree topology)

### 2026-03-05: VESSEL-1005 [PASS]
- Attempted: Implement SoA segment storage, tree topology, and composite VascularTree type
- Result: 77 new tests pass (98 total), all acceptance criteria met
- Regression gate: 98 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/types.jl — SegmentData{V}, TreeTopology, VascularTree, add_segment!, get_children, n_segments
  - test/test_types.jl — 12 testsets with 77 tests
  - src/VesselTree.jl — include types.jl, export public API
  - test/runtests.jl — include test_types.jl
- Learning:
  - Terminal counting: first child of parent is net-zero (parent loses terminal, child gains it); only second child is net +1
  - SoA views (`@view field[1:n]`) work with AK.maximum and other operations
  - Capacity-based pre-allocation works cleanly: error on overflow, views for active region
  - junction_type as Symbol (:none, :bifurcation, :trifurcation) keeps topology CPU-side
- Next: VESSEL-1006 (Morphometric parameters + domain types)

### 2026-03-05: VESSEL-1006 [PASS]
- Attempted: Implement morphometric parameter structs and perfusion domain types
- Result: 425 new tests pass (523 total), all acceptance criteria met, P1 milestone complete
- Regression gate: 523 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/parameters.jl — MorphometricParams struct, kassab_coronary_params(), classify_order()
  - src/domain.jl — AbstractDomain, SphereDomain, BoxDomain, EllipsoidDomain with in_domain/sample_point/signed_distance
  - test/test_parameters.jl — 94 tests covering tables, connectivity matrix, thresholds, classify_order
  - test/test_domain.jl — 331 tests (3 domain types × boundary/inside/outside/sampling)
  - src/VesselTree.jl — includes + exports for parameters and domain modules
- Learning:
  - Connectivity matrix is 12x12 (daughter_order+1 × parent_order+1), column sums ~2-3
  - EllipsoidDomain signed_distance is approximate (exact requires iterative Newton solve)
  - Rejection sampling for sphere/ellipsoid sample_point is adequate for moderate aspect ratios
  - classify_order scans from highest order downward for correct boundary assignment
- Next: VESSEL-1007 (AK-accelerated distance kernels) — P2 CCO Engine begins

### 2026-03-05: VESSEL-1007 [PASS]
- Attempted: Implement AK-accelerated point-to-segment distance kernels
- Result: 1022 new tests pass (1545 total), all acceptance criteria met
- Regression gate: 1545 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/distance.jl — point_segment_distance/sq (scalar), compute_all_distances! (AK.foreachindex), find_nearest_segments (AK.sortperm)
  - test/test_distance.jl — 13 testsets with 1022 tests (1000 from AK vs scalar comparison on 1000 segments)
  - src/VesselTree.jl — include + exports
- Learning:
  - AK.foreachindex with @view for active region works cleanly for distance computation
  - AK.sortperm returns correct permutation for find_nearest_segments
  - Zero-length segment handling: degenerate case returns distance to the single point
  - All 6 geometric cases (on-segment, proximal, distal, perpendicular, behind, beyond) verified
- Next: VESSEL-1008 (AK-accelerated intersection testing)

### 2026-03-05: VESSEL-1008 [PASS]
- Attempted: Implement AK-accelerated segment-segment intersection testing
- Result: 1018 new tests pass (2563 total), all acceptance criteria met
- Regression gate: 2563 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/intersection.jl — segments_min_distance_sq, segments_intersect, check_intersections! (AK.foreachindex), has_any_intersection (AK.any), check_domain_crossing
  - test/test_intersection.jl — 14 testsets with 1018 tests
- Learning:
  - Closest-point-on-two-segments algorithm handles all degenerate cases (parallel, collinear, zero-length)
  - AK.any(identity, bool_view) works correctly for intersection reduction
  - check_domain_crossing samples interior points (not endpoints) to detect boundary crossing
- Next: VESSEL-1009 (Kamiya bifurcation point optimization)

### 2026-03-05: VESSEL-1009 [PASS]
- Attempted: Implement Kamiya bifurcation point optimization with surface cost
- Result: 224 new tests pass (2787 total), all acceptance criteria met
- Regression gate: 2787 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/kamiya.jl — surface_cost_at_t, optimize_bifurcation_point (golden section), evaluate_all_connections! (AK.foreachindex), select_best_connection (AK.minimum + scan)
  - test/test_kamiya.jl — 11 testsets with 224 tests
- Learning:
  - AK.mapreduce doesn't support UnitRange — must use concrete arrays or AK.minimum + manual argmin
  - Surface cost minimum isn't at segment midpoint even for perpendicular terminal — parent radius dominates
  - Golden section search converges in ~20 iterations to high precision
  - Murray's law: r_parent^gamma = r_left^gamma + r_right^gamma verified for both symmetric and asymmetric splits
- Next: VESSEL-1010 (Core CCO growth loop + Murray's law propagation)

### 2026-03-05: VESSEL-1010 [PASS]
- Attempted: Implement core CCO growth loop with Murray's law radius propagation
- Result: 45 new tests pass (2832 total), all acceptance criteria met
- Regression gate: 2832 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/murray.jl — update_radii! (bottom-up Murray's law propagation by generation)
  - src/growth.jl — grow_tree!, add_bifurcation!, sample_terminal_candidate
  - test/test_growth.jl — 9 testsets with 45 tests
- Learning:
  - add_bifurcation! must transfer original children to continuation segment when splitting
  - Adaptive distance threshold d_thresh = domain_size / (10 * (n+1)^(1/3)) works for space-filling
  - Intersection check with min_dist = d_thresh * 0.1 prevents self-intersection
  - Murray's law holds exactly at every bifurcation after update_radii! (verified)
  - Growth rate depends on domain size vs intersection threshold balance
- Next: VESSEL-1011 (Spatial indexing for fast segment queries)

### 2026-03-05: VESSEL-1011 [PASS]
- Attempted: Implement spatial grid for fast segment queries
- Result: 32 new tests pass (2864 total), P2 CCO Engine milestone complete
- Regression gate: 2864 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/spatial.jl — SpatialGrid, build_grid, insert!, query_nearby
  - test/test_spatial.jl — 6 testsets with 32 tests (including 20 no-false-negative checks)
- Learning:
  - Dict-based cells avoid pre-allocating huge 3D array for sparse grids
  - Cell neighborhood search radius = ceil(search_radius / cell_size) cells in each direction
  - Query returns all segments in neighboring cells; caller does exact distance filtering
- Next: VESSEL-1012 (Hemodynamics: Poiseuille resistance + pressure/flow)

### 2026-03-05: VESSEL-1012 [PASS]
- Attempted: Implement hemodynamic computations (Poiseuille resistance, flow, pressure, validation)
- Result: 607 new tests pass (3471 total), all acceptance criteria met
- Regression gate: 3471 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/hemodynamics.jl — compute_resistances! (AK.foreachindex), _topo_order (BFS), compute_flows!, compute_pressures!, validate_hemodynamics
  - test/test_hemodynamics.jl — 16 testsets with 607 tests (including 500-segment AK verification and CCO-grown tree integration)
  - src/VesselTree.jl — include + exports for hemodynamics
- Learning:
  - Generation-based traversal breaks after add_bifurcation! because child generations are not updated when reparented to continuation segment. BFS from root is the correct approach.
  - Flow splits inversely proportional to subtree resistance (conductance-weighted), not equally among terminals
  - Subtree resistance: R_sub = R_self + 1/(sum of child conductances) for bifurcation, R_self for terminal
  - Total flow = (P_root - P_term) / R_total_tree; pressures computed top-down as P_distal = P_proximal - Q*R
  - Poiseuille formula R = 8*mu*L/(pi*r^4) implemented as AK kernel; r^4 scaling means small radius changes have huge resistance effects
- Next: VESSEL-1013 (Flow-weighted radius assignment)

### 2026-03-05: VESSEL-1013 [PASS]
- Attempted: Implement flow-weighted radius assignment with territory-based flow distribution
- Result: 21 new tests pass (3492 total), P3 Hemodynamics milestone complete
- Regression gate: 3492 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/hemodynamics.jl — assign_terminal_flows! (AK nn-distance for territory), recompute_radii_from_flow!, _compute_total_resistance
  - test/test_hemodynamics.jl — 5 new testsets for flow-weighted assignment
  - src/VesselTree.jl — new exports
- Learning:
  - Territory weight via nn_dist^3 gives reasonable volume proxy for 3D Voronoi-like territory
  - Terminal radius proportional to Q^(1/gamma) from Murray's law generalization (Q ∝ r^gamma)
  - AK.foreachindex used for pairwise terminal distance computation; AK.minimum for nn lookup
  - Flow-weighted pipeline: assign_terminal_flows! → recompute_radii_from_flow! → compute_resistances! → compute_flows! → compute_pressures!
- Next: VESSEL-1014 (Kassab: Strahler ordering + asymmetry sampling)

### 2026-03-05: VESSEL-1014 [PASS]
- Attempted: Implement Kassab Strahler ordering, asymmetry sampling, daughter radii computation
- Result: 438 new tests pass (3930 total), all acceptance criteria met
- Regression gate: 3930 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/kassab.jl — assign_strahler_orders! (AK kernel), sample_asymmetry, compute_daughter_radii, apply_kassab_radii!
  - test/test_kassab.jl — 14 testsets with 438 tests (including 200-segment AK verification and KS distribution test)
  - src/VesselTree.jl — include + exports
- Learning:
  - AK.foreachindex works for Strahler ordering: classify_order loop inside kernel closure works correctly
  - Beta(2.5, 0.8) theoretical median is ~0.812, not 0.76 (Kassab's empirical median from the paper)
  - apply_kassab_radii! does top-down asymmetry assignment then bottom-up Murray's law propagation
  - compute_daughter_radii: r_large = r_parent / (1 + asymmetry^gamma)^(1/gamma)
- Next: VESSEL-1015 (Kassab: Connectivity matrix validation)

### 2026-03-05: VESSEL-1015 [PASS]
- Attempted: Implement connectivity matrix construction and validation against Kassab 1993
- Result: 24 new tests pass (3954 total), all acceptance criteria met
- Regression gate: 3954 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/kassab.jl — build_empirical_connectivity, validate_connectivity (chi-squared), constrain_connectivity!, _find_nearest_valid_order
  - test/test_kassab.jl — 4 new testsets for connectivity matrix
  - src/VesselTree.jl — new exports
- Learning:
  - CCO-generated trees have very small segments classified as order 0 that can be parents (unlike real vasculature)
  - Chi-squared test only meaningful for cells with expected count >= 0.5 (sparse cells skipped)
  - constrain_connectivity! adjusts daughter radius to mean diameter of nearest valid order
  - Chisq distribution from Distributions.jl used for p-value computation
- Next: VESSEL-1016 (Kassab: Integrated constraint enforcement in growth loop)

### 2026-03-05: VESSEL-1016 [PASS]
- Attempted: Integrate Kassab constraints into CCO growth loop
- Result: 40 new tests pass (3994 total), P4 Kassab Constraints milestone complete
- Regression gate: 3994 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/growth.jl — grow_tree! gains kassab=true keyword for asymmetric daughter radii
  - test/test_growth.jl — 5 new testsets for Kassab-integrated growth
- Learning:
  - Simple approach: sample asymmetry per bifurcation, assign r_large to continuation, r_small to new terminal
  - Murray's law update_radii! at end ensures consistency regardless of initial assignment
  - Backward compatible: kassab=false (default) preserves old behavior
- Next: VESSEL-1017 (Barabasi: Chi/rho parameters + junction geometry)

### 2026-03-05: VESSEL-1017 [PASS]
- Attempted: Implement Barabasi chi/rho parameters and junction geometry
- Result: 145 new tests pass (4139 total), all acceptance criteria met
- Regression gate: 4139 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/barabasi.jl — compute_chi, compute_rho, classify_junction, compute_junction_angles, apply_junction_geometry!, _find_perpendicular, _update_seg_length!
  - test/test_barabasi.jl — 15 testsets with 145 tests
  - src/VesselTree.jl — include + exports
- Learning:
  - Sprouting regime (rho < 0.83): large daughter continues straight, small at ~90 deg
  - Branching regime (rho >= 0.83): linear interpolation from threshold to symmetric angle
  - _find_perpendicular uses cross product with axis of smallest component
  - apply_junction_geometry! preserves segment lengths while rotating endpoints
  - compute_chi(1.0, π) fails because π is Irrational — must use Float64(π)
- Next: VESSEL-1018 (Barabasi: Trifurcation detection + handling)

### 2026-03-05: VESSEL-1018 [PASS]
- Attempted: Implement trifurcation detection and handling (Barabasi chi > 0.83)
- Result: 120 new tests pass (4259 total), all acceptance criteria met
- Regression gate: 4259 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/barabasi.jl — compute_trifurcation_angles, check_trifurcation_merge, merge_to_trifurcation!, apply_trifurcation_geometry!
  - src/growth.jl — grow_tree! gains trifurcation=true keyword
  - test/test_trifurcation.jl — 14 testsets with 120 tests
  - src/VesselTree.jl — new exports
  - test/runtests.jl — include test_trifurcation.jl
- Learning:
  - check_trifurcation_merge scans all bifurcations, computes chi = 2*pi*r/d to distal point
  - merge_to_trifurcation! delegates to add_segment! which handles the bifurcation→trifurcation topology upgrade
  - apply_trifurcation_geometry! places 3 daughters at azimuthal angles 0, 2pi/3, 4pi/3 with polar angle proportional to inverse radius
  - Trifurcation rate depends on tree density and parent radii; in dense trees more merges occur
  - Murray's law verified at all trifurcations: r_parent^gamma = r1^gamma + r2^gamma + r3^gamma
- Next: VESSEL-1019 (Barabasi: Surface cost function + integration into growth)

### 2026-03-05: VESSEL-1019 [PASS]
- Attempted: Implement surface cost with junction estimate and integrate Barabasi geometry into growth
- Result: 68 new tests pass (4327 total), P5 Barabasi Geometry milestone complete
- Regression gate: 4327 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/kamiya.jl — junction_surface_estimate, surface_cost_with_junction
  - src/growth.jl — grow_tree! gains barabasi=true keyword; applies apply_junction_geometry! after bifurcation
  - test/test_surface_cost.jl — 10 testsets with 68 tests
  - test/runtests.jl — include test_surface_cost.jl
- Learning:
  - Sprouting junctions have lower surface cost: base * r_small * 0.5 vs base * (r_large + r_small) * 0.5
  - Junction cost is additive to tube cost in surface_cost_with_junction
  - Barabasi geometry application rotates daughter endpoints after radii update
  - With barabasi=true, intersection checks may reject more candidates (rotated segments can intersect)
  - Angle distribution from grown trees is measurable; mean angle is reasonable
- Next: VESSEL-1020 (Validation framework + diameter/asymmetry tests)

### 2026-03-05: VESSEL-1020 [PASS]
- Attempted: Implement validation framework with diameter, asymmetry, L/D, and angle statistics
- Result: 130 new tests pass (4457 total), all acceptance criteria met
- Regression gate: 4457 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/validation.jl — ValidationReport, compute_diameters! (AK), compute_ld_ratios! (AK), compute_asymmetry_ratios, compute_branching_angles, validate_tree, print_report
  - test/test_validation.jl — 11 testsets with 130 tests
  - src/VesselTree.jl — include + exports
  - test/runtests.jl — include test_validation.jl
- Learning:
  - AK.foreachindex for bulk diameter and L/D computation works cleanly
  - Asymmetry and angle computation sequential (iterate bifurcations only)
  - Approximate z-test p-value for diameter validation (full KS would need HypothesisTests.jl)
  - ValidationReport uses @kwdef for convenient construction
- Next: VESSEL-1021 (Full validation suite: connectivity + trifurcation + P(lambda))

### 2026-03-05: VESSEL-1021 [PASS]
- Attempted: Complete validation suite with connectivity, trifurcation %, P(lambda), Murray deviation
- Result: 51 new tests pass (4508 total), P6 Validation milestone complete
- Regression gate: 4508 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/validation.jl — compute_murray_deviation, compute_trifurcation_pct, compute_lambda_distribution, extra dict in ValidationReport, updated validate_tree and print_report
  - test/test_validation_full.jl — 10 testsets with 51 tests
  - test/runtests.jl — include test_validation_full.jl
- Learning:
  - validate_connectivity takes (empirical_matrix, reference_matrix) not (matrix, params)
  - P(lambda) = min_daughter_length / (2*pi*r_parent) measures how close a junction is to trifurcation
  - Murray deviation at all junctions after update_radii! is < 1e-6 as expected
  - Extra metrics stored in Dict{Symbol, Any} for extensibility
- Next: VESSEL-1022 (Territory partitioning + collision avoidance)

### 2026-03-05: VESSEL-1022 [PASS]
- Attempted: Implement territory partitioning, collision avoidance, coronary configs
- Result: 24 new tests pass (4532 total), all acceptance criteria met
- Regression gate: 4532 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/forest.jl — TreeConfig, TerritoryMap, initialize_territories, query_territory, sample_in_territory, check_inter_tree_collision, coronary_tree_configs
  - test/test_forest.jl — 9 testsets with 24 tests
  - src/VesselTree.jl — include + exports
  - test/runtests.jl — include test_forest.jl
- Learning:
  - Territory map uses ~20 cells per dimension for coarse Voronoi-like partition
  - Nearest-root assignment gives reasonable territory shapes
  - check_inter_tree_collision uses AK.minimum on distance buffer
  - coronary_tree_configs: LAD 40%, LCX 25%, RCA 35% with approximate ostia positions
- Next: VESSEL-1023 (Multi-tree forest growth)

### 2026-03-05: VESSEL-1023 [PASS]
- Attempted: Implement multi-tree forest growth with round-robin and collision avoidance
- Result: 46 new tests pass (4578 total), P7 Multi-Tree Forest milestone complete
- Regression gate: 4578 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/forest.jl — CoronaryForest, generate_coronary_forest (round-robin + territory + collision), validate_forest
  - test/test_forest_growth.jl — 7 testsets with 46 tests
  - src/VesselTree.jl — new exports
  - test/runtests.jl — include test_forest_growth.jl
- Learning:
  - Round-robin growth: iterate over configs, try one terminal per tree per round
  - Inter-tree collision: check point distance to all other trees' segments via AK.minimum
  - Root segment direction/length from TreeConfig; capacity estimated at 3x target_terminals
  - validate_forest runs validate_tree on each sub-tree
- Next: VESSEL-1024 (JLD2 save/load + VTP export)

### 2026-03-05: VESSEL-1028 [PASS]
- Attempted: Integrate SpatialGrid into grow_tree! and generate_coronary_forest for grid-accelerated growth
- Result: 237 new tests pass (6394 total), all acceptance criteria met
- Regression gate: 6394 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/growth.jl — _domain_size, _grid_cell_size, _insert_segment_to_grid!, _find_nearest_via_grid, grow_tree! with grid
  - src/forest.jl — generate_coronary_forest with per-tree spatial grids
  - test/test_grid_growth.jl — 8 testsets with 237 tests (grid vs brute force, Murray's law, connectivity, kassab, forest)
  - test/runtests.jl — include test_grid_growth.jl
- Learning:
  - Grid activation threshold (n >= 200) preserves exact brute-force results for small trees, matching existing tests with same RNG seed
  - Grid rebuilt every 100 additions; incremental insert between rebuilds (stale entries are OK since distance is computed from current SoA arrays)
  - search_radius = d_thresh * 5.0 provides good coverage; falls back to brute force if too few candidates
  - _domain_size helper extracts domain sizing logic shared between grow_tree! and generate_coronary_forest
- Next: VESSEL-1029 (Statistical subdivision engine)

### 2026-03-05: VESSEL-1029 [PASS]
- Attempted: Implement statistical subdivision engine using Kassab connectivity matrix
- Result: 77 new tests pass (6471 total), all acceptance criteria met
- Regression gate: 6471 tests pass, 0 fail, 0 error
- Files created/modified:
  - src/subdivision.jl — estimate_total_segments, estimate_subdivision_capacity, subdivide_terminals!, _subdivide_recursive!, _add_subdivision_daughter!, _random_daughter_direction
  - src/VesselTree.jl — include + exports
  - test/test_subdivision.jl — 12 testsets with 77 tests
  - test/runtests.jl — include test_subdivision.jl
- Learning:
  - Poisson sampling for daughter counts preserves CM expected values; round() loses too many small entries
  - Cascading bifurcations needed when Poisson gives >3 daughters: use slot-aware approach checking actual child count before each add
  - E[S(5)] ≈ 78, E[S(7)] ≈ 259, E[S(11)] ≈ 1506 segments per terminal with current CM values
  - CM values are conservative; full 6M segments would need larger CM entries or more CCO terminals
  - _find_perpendicular from barabasi.jl reused for direction computation
  - add_segment! errors at 4th child — must cascade before 3rd child used up
- Next: VESSEL-1030 (Post-hoc Kassab radius refinement)
