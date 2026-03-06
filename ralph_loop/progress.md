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
