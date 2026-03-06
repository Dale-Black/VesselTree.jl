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
