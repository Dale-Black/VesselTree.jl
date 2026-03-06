# VesselTree.jl

Physiologically accurate vascular tree generation from large arteries (~3-5mm) down to capillaries (8um cutoff).

## Quick Reference

```bash
# Run tests
julia --project=. test/runtests.jl

# Load package
julia --project=. -e 'using VesselTree'
```

## Architecture

- **CCO** (Constrained Constructive Optimization): Space-filling tree growth algorithm
- **Kassab-Molloi Morphometry**: 12 Strahler orders, diameter 8um-4500um, ~6M segments
- **Barabasi Surface Optimization**: Junction geometry (sprouting vs branching)

## Key Conventions

- **AcceleratedKernels.jl (AK)**: ALL computational kernels must use AK. Import as `import AcceleratedKernels as AK`. Key functions: `foreachindex`, `reduce`, `mapreduce`, `sort`, `sortperm`, `any`, `all`, `count`, `accumulate`. No `argmin`/`argmax` — use `mapreduce` with `neutral` parameter instead.
- **SoA layout**: Structure-of-Arrays for GPU compatibility. Separate arrays per field (e.g., `proximal_x`, `proximal_y`), not arrays of structs.
- **Murray's law**: gamma = 7/3 (Huo-Kassab 2007), NOT 3.0
- **Vessel cutoff**: 8um (capillary diameter)
- **Pre-allocation + views**: Growing arrays use pre-allocated capacity with `@view field[1:n]` for AK operations.

## AK Gotchas

- AK has ~20-70us overhead per call; crossover vs Base at ~100K elements
- `BitVector` not supported — convert to `Vector{Bool}` first
- `merge_sort_by_key!` is GPU-only; use `sortperm` on CPU
- `mapreduce` requires `neutral` keyword (e.g., `neutral=(Inf, 0)`)
- `sqrt` of negative in AK kernel throws `DomainError` — guard inputs explicitly

## Scientific Constants (src/VesselTree.jl)

- `MURRAY_GAMMA = 7/3` — vascular branching exponent
- `VESSEL_CUTOFF_UM = 8.0` — capillary diameter in micrometers
- `BLOOD_VISCOSITY = 0.0035` — Pa*s (3.5 cP)
- `ROOT_PRESSURE = 13332.0` — Pa (100 mmHg)
- `TERMINAL_PRESSURE = 3999.6` — Pa (30 mmHg)
- `TRIFURCATION_CHI_TH = 0.83` — Barabasi trifurcation threshold
- `SPROUTING_RHO_TH = 0.83` — Barabasi sprouting threshold
- `ASYMMETRY_ALPHA = 2.5`, `ASYMMETRY_BETA = 0.8` — Beta distribution params

## Dependencies

- AcceleratedKernels.jl (0.4) — cross-architecture parallel algorithms
- Distributions.jl (0.25) — statistical distributions (Beta for asymmetry)
- StaticArrays.jl (1) — fixed-size arrays for 3D coordinates
- LinearAlgebra, Random (stdlib)
