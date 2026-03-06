# VesselTree.jl — Guardrails & Anti-Cheat Rules

## RULE #1: AcceleratedKernels.jl IS MANDATORY — NO EXCEPTIONS

Every function that operates on arrays of segment data MUST use AK.jl.

**FORBIDDEN in computational kernels:**
```julia
# NEVER do this for parallel-eligible work:
for i in 1:n_segments
    distances[i] = compute_distance(...)
end

# NEVER do this:
@threads for i in eachindex(segments)
    ...
end

# NEVER do this:
map(f, segments)  # when f is a computational kernel

# NEVER use broadcasting as a substitute for AK:
distances .= point_distance.(segments, Ref(candidate))  # NO
```

**REQUIRED — the AK pattern:**
```julia
import AcceleratedKernels as AK

# Element-wise parallel operations
AK.foreachindex(output_array) do i
    output_array[i] = compute_something(input_arrays..., i)
end

# Reduction operations
result = AK.reduce(+, array)
idx = AK.argmin(array)
```

**Allowed WITHOUT AK (must be commented `# sequential: <reason>`):**
- Tree topology traversal (parent/child walks)
- Building adjacency/children dictionaries
- Topological sort computation
- File I/O, export formatting, JSON parsing
- One-time setup (parameter construction, domain initialization)
- The outer CCO growth loop itself (inherently sequential: add one terminal at a time)

**The test:** If a function takes an array and produces an array, it needs AK. If it walks a tree structure following pointers, it can be sequential.

---

## RULE #2: SoA DATA LAYOUT FOR GEOMETRY

Segment geometry MUST use Structure-of-Arrays layout for AK/GPU compatibility.

**FORBIDDEN:**
```julia
# AoS — kills GPU performance, incompatible with AK
struct Segment
    proximal::SVector{3, Float64}
    distal::SVector{3, Float64}
    radius::Float64
end
segments = Vector{Segment}(undef, n)
```

**REQUIRED:**
```julia
# SoA — GPU-friendly, AK-compatible
mutable struct SegmentData{V <: AbstractVector}
    proximal_x::V; proximal_y::V; proximal_z::V
    distal_x::V;   distal_y::V;   distal_z::V
    radius::V
    # ...
    n::Int  # number of active segments
end
```

Topology (parent_id, child_ids, strahler_order, is_terminal) stays in `Vector{Int32}` / `Vector{Bool}` — these are sequential-access only.

---

## RULE #3: MURRAY EXPONENT IS 2.33, NOT 3.0

The bifurcation law is: `r_parent^gamma = r_daughter1^gamma + r_daughter2^gamma`

**gamma = 7/3 ~ 2.333** (Huo & Kassab 2007)

NOT Murray's classical 3.0. The meta-analysis pooled value is 2.39 (95% CI: 2.24-2.54). We use 7/3 as the theoretically derived value.

If you find yourself writing `gamma = 3.0`: STOP. You are wrong.

---

## RULE #4: VESSEL CUTOFF IS 8 MICROMETERS

The minimum vessel diameter is **8 um** (capillary level, Kassab Order 0).

NOT 50um. NOT 100um. The full tree goes all the way to capillaries.

This means ~2 million segments per major coronary artery. AK.jl is essential.

---

## RULE #5: REGRESSION GATE IS MANDATORY

Before EVERY commit that changes source code:

```bash
cd /Users/daleblack/Documents/dev/julia/VesselTree.jl
julia --project=. test/runtests.jl 2>&1 | tail -10
```

**Test count is a RATCHET**: it must never decrease.

| After Story   | Minimum Tests |
|--------------|---------------|
| VESSEL-1005  | 20            |
| VESSEL-1006  | 40            |
| VESSEL-1007  | 60            |
| VESSEL-1010  | 100           |
| VESSEL-1014  | 150           |
| VESSEL-1017  | 200           |
| VESSEL-1019  | 250           |
| VESSEL-1021  | 280           |
| VESSEL-1025  | 350           |
| VESSEL-1033  | 21000         |
| VESSEL-1035  | 21050         |
| VESSEL-1036  | 21100         |
| VESSEL-1037  | 21150         |
| VESSEL-1038  | 21200         |
| VESSEL-1040  | 21250         |
| VESSEL-1043  | 21300         |

---

## RULE #6: KASSAB DATA IS GROUND TRUTH — NO FABRICATION

**ZERO TOLERANCE for fabricated, approximate, or made-up data.**

Every number in `parameters.jl` MUST trace to a specific table and page number in a published paper:
- **Kassab et al. 1993** (Am J Physiol Heart 265:H350-H365) — Tables 1-9
- **Jiang et al. 1994** (J Appl Physiol 76:882-892) — diameter-defined Strahler method
- **Kassab & Fung 1995** (Ann Biomed Eng 23:13-20) — Murray's law validation
- **Huo & Kassab 2007/2009** — gamma=7/3, volume scaling

PDFs at: `/Users/daleblack/Documents/vessel tree growing papers/`

### Three SEPARATE connectivity matrices (NOT one generic "coronary" CM):
- RCA: Kassab 1993 Table 6 (11 orders)
- LAD: Kassab 1993 Table 7 (11 orders)
- LCX: Kassab 1993 Table 8 (10 orders)

### Real data (SEGMENTS) from Kassab 1993:
```
RCA Segments:
Order  D(um)±SD       L(um)±SD       N_seg
1     10.8±2.85      96±78          854,880
2     21.5±3.48      92±74          163,233
3     39.1±5.64      100±83          43,858
4     69.2±10.5      123±111         12,285
5     124±20.2       176±166          5,449
6     228±40.4       311±299          3,279
7     408±68.8       503±440          1,563
8     703±118        810±720            550
9    1211±210       1309±1260           135
10   2167±290       1997±1750            36
11   3480±440       2610±1870             6
```

### Key distinction: Segments vs Elements
- **Segment**: vessel between two bifurcation nodes
- **Element**: series of same-order segments (Kassab's fundamental unit)
- The **connectivity matrix operates on ELEMENTS**, not segments
- Element diameter = mean of segment diameters; element length = sum of segment lengths
- Validate against ELEMENT-level distributions

See `ralph_loop/research/kassab_1993_real_data.md` for all three arteries' complete tables.

Generated trees MUST be validated against these distributions. A tree that doesn't match Kassab is wrong, no matter how pretty it looks.

---

## RULE #7: BARABASI GEOMETRY IS NOT OPTIONAL

Junction geometry MUST follow the Barabasi 2026 surface optimization rules:

1. Compute `rho = r_small / r_large` at every bifurcation
2. If `rho < 0.83`: **SPROUTING** — small branch at ~90 deg, large continues straight
3. If `rho >= 0.83`: **BRANCHING** — both daughters deflect, angle depends on rho
4. Compute `chi = w/r` — if `chi > 0.83`, check for stable trifurcation

Do NOT default to symmetric 120-degree bifurcations everywhere. Real vascular trees show a bimodal angle distribution (sprouting peak + branching peak).

---

## RULE #8: ONE STORY AT A TIME

Do NOT start a story until ALL of its `blockedBy` stories have status `"done"`.

Do NOT work on multiple stories simultaneously.

Do NOT skip ahead to "more interesting" stories.

---

## RULE #9: CHECKPOINT EVERY 15 MINUTES

You WILL be killed by a timeout. Protect your work:

1. After completing ANY test file: commit
2. After implementing ANY feature: commit
3. After 15 minutes of work: commit (even WIP)
4. BEFORE running long Julia commands: commit

```bash
cd /Users/daleblack/Documents/dev/julia/VesselTree.jl
git add [changed files]
git commit -m "VESSEL-1XXX: WIP — [what you did so far]"
```

---

## RULE #10: DISCOVERY STORIES DO NOT WRITE SOURCE CODE

DISCOVERY stories (VESSEL-1000 through 1003) produce ONLY:
- Markdown files in `ralph_loop/research/`
- Updates to `prd.json` (refining stories based on findings)
- No changes to `src/` or `test/`
- No changes to `Project.toml`

---

## RULE #11: TESTS MUST VERIFY BEHAVIOR, NOT JUST CONSTRUCTION

**FORBIDDEN (trivial tests):**
```julia
tree = VascularTree()
@test tree isa VascularTree  # Tests nothing useful
```

**REQUIRED (behavioral tests):**
```julia
# Test Murray's law is actually enforced
seg = add_test_bifurcation!(tree, r_parent=100.0)
update_radii!(tree, params)
r1, r2 = get_daughter_radii(tree, seg)
@test r_parent^2.33 ≈ r1^2.33 + r2^2.33  rtol=1e-10

# Test AK kernel produces correct distances
distances = similar(seg_radius)
compute_distances!(distances, seg_data, candidate)
# Verify against scalar reference implementation
for i in 1:n
    expected = scalar_point_segment_distance(candidate, seg_data, i)
    @test distances[i] ≈ expected  rtol=1e-12
end
```

---

## RULE #12: GENERAL-PURPOSE DESIGN

VesselTree.jl is NOT coronary-only. The architecture must support:
- **Coronary**: LAD, LCX, RCA (initial target)
- **Cerebral**: Circle of Willis, MCA, ACA, PCA
- **Pulmonary**: Pulmonary arterial tree

This means:
- Morphometric parameters are a configurable struct, not hardcoded
- Domain shapes are abstract (ellipsoid, mesh, SDF)
- `kassab_coronary_params()` is one preset among potentially many
- Tree names are strings, not hardcoded `:LAD`/`:LCX`/`:RCA`

Do NOT hardcode coronary-specific logic into the core algorithm.

---

## Anti-Cheat Patterns

### Tests must exercise AK kernels on real data
```julia
# BAD: Only tests scalar helper
@test point_segment_distance(p, a, b) ≈ 5.0

# GOOD: Tests AK kernel on array data
n = 1000
seg_data = random_segment_data(n)
distances = similar(seg_data.radius)
compute_distances!(distances, seg_data, candidate)
@test length(distances) == n
@test all(distances .>= 0)
@test argmin(distances) == brute_force_nearest(seg_data, candidate)
```

### Validation must use statistical tests, not eyeballing
```julia
# BAD: "Looks about right"
@test mean(diameters) > 0

# GOOD: KS test against Kassab distribution
using HypothesisTests
ks = ApproximateTwoSampleKSTest(generated_diameters, kassab_diameters)
@test pvalue(ks) > 0.05
```

### No hardcoded test expectations
```julia
# BAD: Fragile — breaks if algorithm improves
@test tree.n_terminals == 1000

# GOOD: Tests invariant
@test tree.n_terminals >= target_terminals * 0.95  # At least 95% success rate
@test all(seg_data.radius[1:tree.n] .> 0)          # All radii positive
```

---

## RULE #13: DIAMETER-DEFINED STRAHLER ORDERING (JIANG 1994)

Do NOT use simple diameter binning for Strahler order assignment. Kassab uses an ITERATIVE
diameter-defined Strahler system (Jiang et al. 1994):

1. Initialize with topological Strahler ordering
2. Compute per-order mean D_n and SD_n
3. Redefine bounds: lower = [(D_{n-1}+SD_{n-1}) + (D_n-SD_n)]/2, upper = [(D_n+SD_n) + (D_{n+1}-SD_{n+1})]/2
4. Reassign orders by diameter range
5. Iterate until < 1% of vessels change order

See `ralph_loop/research/jiang_1994_methodology.md` for full algorithm.

---

## RULE #14: NO CASCADE STUBS IN SUBDIVISION

Subdivision must produce ONLY bifurcations (2 children). When the CM says a parent element has
N daughters, create a chain of N-1 bifurcations (each peeling off one daughter), NOT a single
node with N children connected by tiny stub segments.

This naturally produces the S/E ratio (multiple segments per element) and keeps trifurcation
rate below the Barabasi target of ~8.3%.

---

## RULE #15: ELEMENT-LEVEL VALIDATION

All Kassab validation must operate on ELEMENTS (series of same-order segments), not on
individual segments. The connectivity matrix, diameter distributions, length distributions,
and asymmetry ratios in Kassab 1993 are all element-level measurements.

Implement element grouping BEFORE running validation statistics.

---

## RULE #16: P10 (EXACT KASSAB PARITY) IS THE CURRENT PRIORITY

P10 stories (VESSEL-1034 through VESSEL-1043) have priority 9 and are the current focus.

Previously, P9 stories (VESSEL-1028 through VESSEL-1033) had priority 7 and ran before P8 stories (priority 8).

P9 uses a **hybrid CCO + statistical subdivision** approach:
1. **CCO phase**: Grow upper-order skeleton (orders 5-11) with spatial grid acceleration. ~1000-2000 terminals per artery. This gives realistic spatial layout with intersection checking.
2. **Subdivision phase**: Recursively split each terminal down to order 0 (8um capillaries) using Kassab's connectivity matrix. Each order-k terminal gets daughters per CM[m+1, k+1]. NO intersection checking at capillary scale — pure statistical tree generation.
3. **Refinement phase**: Apply Kassab asymmetry radii POST-HOC (not during growth). Apply Barabasi junction geometry.

This produces ~4-6M total segments in < 5 minutes. The connectivity matrix IS the fast path — it encodes the complete branching statistics directly.

**Key insight**: You do NOT need to grow millions of capillaries one-at-a-time with CCO. The CM tells you exactly how many daughters of each order to create. CCO handles the physics of large-vessel layout; the CM handles the statistics of microvasculature branching.

Do NOT attempt to grow millions of terminals with the CCO loop. That is O(n²) and will take weeks. Use subdivision.

---

## Quick Reference

```bash
# Working directory
cd /Users/daleblack/Documents/dev/julia/VesselTree.jl

# Run all tests
julia --project=. test/runtests.jl

# Run single test file
julia --project=. -e 'using Test; include("test/test_types.jl")'

# Package load test
julia --project=. -e 'using VesselTree; println("OK")'

# Instantiate after Project.toml change
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Commit pattern
git add [files] && git commit -m "VESSEL-1XXX: Description"
```
