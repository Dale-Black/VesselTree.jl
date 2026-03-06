# VesselTree.jl — Physiologically Accurate Vascular Tree Generator

You are an autonomous agent building **VesselTree.jl**, a Julia package that generates physiologically accurate vascular trees from large arteries (~3-5mm) down to capillaries (8um cutoff). The package is general-purpose (coronary, cerebral, pulmonary) with the initial implementation targeting coronary vasculature.

## What You Are Building

Three complementary algorithms fused into one generator:

1. **Constrained Constructive Optimization (CCO)** — Space-filling tree growth via iterative terminal addition
2. **Kassab-Molloi Morphometry** — Experimentally measured coronary statistics (diameter/length by Strahler order, connectivity matrices, asymmetry distributions)
3. **Barabasi Surface Optimization** — Physics of junction geometry (sprouting vs branching regimes, stable trifurcations when chi > 0.83)

The result: synthetic vascular trees that match both global statistical distributions AND local geometric principles.

## CRITICAL: AcceleratedKernels.jl Is MANDATORY

**Every computational kernel MUST use AcceleratedKernels.jl (AK).** This is non-negotiable.

The package must be CPU/GPU-native from day one. With ~2 million segments in a full coronary tree, plain loops are unacceptable for:
- Distance computations (finding nearest segment to candidate point)
- Intersection testing (checking new segments against existing tree)
- Cost evaluation (surface/volume cost across all candidate connections)
- Radius propagation (Murray's law updates across tree levels)
- Validation statistics (diameter histograms, angle distributions)
- Mesh generation (cylinder vertices for export)

### AK Pattern — The Only Way

```julia
import AcceleratedKernels as AK

# CORRECT: AK.foreachindex for parallel element-wise operations
function compute_distances!(distances, seg_px, seg_py, seg_pz, seg_dx, seg_dy, seg_dz, cx, cy, cz)
    AK.foreachindex(distances) do i
        distances[i] = point_segment_distance(
            cx, cy, cz,
            seg_px[i], seg_py[i], seg_pz[i],
            seg_dx[i], seg_dy[i], seg_dz[i]
        )
    end
end

# CORRECT: AK.reduce for parallel aggregation
best_idx = AK.argmin(distances)

# FORBIDDEN: Plain for loops for computational work
# for i in 1:n_segments  <-- NEVER DO THIS for parallel-eligible work
```

### What IS Allowed Without AK

Sequential tree topology operations that are inherently serial:
- Tree traversal (walking parent→child or child→parent)
- Building adjacency lists
- Topological sorting
- I/O and export formatting

These MUST be clearly marked with a `# sequential: <reason>` comment.

## Architecture

```
VesselTree.jl/
├── src/
│   ├── VesselTree.jl           # Main module
│   ├── types.jl                # SoA data structures (GPU-compatible)
│   ├── parameters.jl           # Kassab + Barabasi parameters
│   ├── domain.jl               # Perfusion domain (SDF, sampling)
│   ├── spatial.jl              # Spatial indexing (grid, AK-accelerated)
│   ├── distance.jl             # Point-segment distance kernels (AK)
│   ├── intersection.jl         # Segment intersection tests (AK)
│   ├── kamiya.jl               # Bifurcation point optimization
│   ├── growth.jl               # Main CCO growth loop
│   ├── murray.jl               # Murray's law radius propagation
│   ├── hemodynamics.jl         # Poiseuille resistance, pressure/flow
│   ├── kassab.jl               # Strahler ordering, asymmetry, connectivity
│   ├── barabasi.jl             # Chi/rho parameters, junction geometry
│   ├── forest.jl               # Multi-tree (LAD+LCX+RCA)
│   ├── validation.jl           # Statistical validation vs Kassab
│   ├── export.jl               # VTP, STL, JLD2 output
│   └── visualization.jl        # Makie plotting
├── test/
│   ├── runtests.jl
│   ├── test_types.jl
│   ├── test_distance.jl
│   ├── test_intersection.jl
│   ├── test_growth.jl
│   ├── test_kassab.jl
│   ├── test_barabasi.jl
│   ├── test_hemodynamics.jl
│   ├── test_validation.jl
│   └── test_forest.jl
├── examples/
├── data/kassab_tables/
└── ralph_loop/
```

## Data Layout: Structure of Arrays (SoA)

For AK.jl compatibility and GPU performance, geometry data uses SoA layout:

```julia
mutable struct SegmentData{V <: AbstractVector{Float64}}
    # Geometry (SoA — each field is a contiguous array)
    proximal_x::V;  proximal_y::V;  proximal_z::V
    distal_x::V;    distal_y::V;    distal_z::V
    radius::V
    seg_length::V

    # Hemodynamics
    flow::V
    pressure_proximal::V
    pressure_distal::V
    resistance::V

    # Count
    n::Int
end
```

Topology (parent/child IDs, Strahler orders, flags) stays in plain `Vector{Int32}` / `Vector{Bool}` on CPU — these are accessed sequentially during tree walks.

## Workflow — Every Iteration

### STEP 0: Environment Check

```bash
cd /Users/daleblack/Documents/dev/julia/VesselTree.jl
julia --project=. -e 'using VesselTree'  # Must load (after VESSEL-1004)
```

### STEP 1: Read State (DO THIS FIRST)

Read these files in order:
1. `ralph_loop/prd.json` — story list + status
2. `ralph_loop/guardrails.md` — rules you MUST follow
3. `ralph_loop/progress.md` — previous learnings

### STEP 2: Pick a Story

Find the highest-priority UNBLOCKED story:
- Status must be `"open"`
- All `blockedBy` stories must have status `"done"`
- Priority 0 > 1 > 2 (lower number = higher priority)
- Ties: pick lowest story ID

### STEP 3: Execute

#### DISCOVERY stories (VESSEL-1000 through 1003)
1. Use WebFetch / WebSearch to study the specified sources
2. Read reference codebases via `gh api` or WebFetch
3. Write findings to `ralph_loop/research/{topic}.md`
4. Do NOT write VesselTree.jl source code
5. Do NOT modify Project.toml
6. DO update prd.json if research reveals the planned approach needs adjustment

#### IMPL stories (most stories)
1. Write tests FIRST when feasible
2. Implement using AK.jl for all computational kernels
3. Run tests: `julia --project=. test/runtests.jl`
4. All tests must pass (including previously passing tests)
5. Run regression gate before committing

#### VERIFY stories
1. Run full test suite
2. Run validation against Kassab data
3. Document results
4. Fix any issues found

### STEP 4: Regression Gate (BEFORE every commit)

```bash
cd /Users/daleblack/Documents/dev/julia/VesselTree.jl
julia --project=. test/runtests.jl 2>&1 | tail -5
# Test count must be >= previous count (ratchet — never goes down)
```

### STEP 5: Update PRD + Progress

After completing a story:

1. Update `prd.json`:
   - Change story status: `"open"` -> `"done"`
   - Add `completedNote` field

2. Append to `progress.md`:
```markdown
### YYYY-MM-DD: VESSEL-1XXX [PASS/FAIL]
- Attempted: what you tried
- Result: what happened
- Regression gate: X tests pass (was Y)
- Learning: key insight for future agents
- Next: what comes next
```

### STEP 6: Commit

```bash
cd /Users/daleblack/Documents/dev/julia/VesselTree.jl
git add [specific changed files]
git commit -m "VESSEL-1XXX: Description"
```

### STEP 7: Exit

- If all stories done: output `RALPH_COMPLETE`
- If blocked by external issue: output `RALPH_BLOCKED` with explanation
- Otherwise: exit normally (loop.sh will restart you)

## Key Scientific Constants

```
Murray exponent (gamma):     7/3 = 2.33  (NOT 3.0! Huo-Kassab 2007)
Vessel cutoff diameter:      8 um (capillary level)
Root pressure:               100 mmHg (aortic)
Terminal pressure:            30 mmHg (arteriolar)
Blood viscosity:              3.5 cP (0.0035 Pa*s)
Trifurcation threshold (chi): 0.83 (Barabasi 2026)
Sprouting threshold (rho):    0.83 (blood vessels)
Asymmetry ratio median:       0.76 (Kassab 1993)
Length/Diameter ratio:         16.8 (Kassab 1993)
```

## Key References

1. **Schreiner & Buxbaum 1993** — CCO algorithm (IEEE Trans Biomed Eng)
2. **Kassab et al. 1993** — Coronary morphometry (Am J Physiol Heart 265:H350-H365)
3. **Jiang et al. 1994** — Diameter-defined Strahler system (J Appl Physiol 76:882-892)
4. **Kassab & Fung 1995** — Arteriolar bifurcation + Murray's law (Ann Biomed Eng 23:13-20)
5. **Huo & Kassab 2007** — gamma = 7/3 derivation (Biophys J)
6. **Huo & Kassab 2009** — Vascular volume scaling law (Biophys J 96:347-353)
7. **Barabasi et al. 2026** — Surface optimization (Nature 649:315-322)
8. **OpenCCO** — https://github.com/OpenCCO-team/OpenCCO (C++ reference)
9. **min-surf-netw** — https://github.com/Barabasi-Lab/min-surf-netw (Python reference)

### Paper PDFs (local)
All papers are available as PDFs for direct reading:
```
/Users/daleblack/Documents/vessel tree growing papers/
  kassab-et-al-1993-morphometry-of-pig-coronary-arterial-trees.pdf
  jiang-et-al-1994-diameter-defined-strahler-system-and-connectivity-matrix-of-the-pulmonary-arterial-tree.pdf
  BF02368296.pdf  (Kassab & Fung 1995)
  main.pdf         (Huo & Kassab 2009)
```

### Research Notes (extracted data)
```
ralph_loop/research/
  kassab_1993_real_data.md    — All tables from Kassab 1993 (needs PDF verification)
  jiang_1994_methodology.md   — Diameter-defined Strahler algorithm
  kassab_parity_gaps.md        — Gap analysis: what's wrong and how to fix it
```

## Checkpoint Discipline

Commit after:
- Every test file that passes
- Every feature implementation that works
- Every 15 minutes of work (even if incomplete — WIP commits are fine)
- BEFORE running long commands

Never lose more than 15 minutes of work to a timeout.

## P9: Kassab-Scale Realism — The Critical Milestone

P9 stories (VESSEL-1028 through VESSEL-1033, priority 7) run BEFORE P8 (Export, priority 8).

**THE CORE ALGORITHM — Hybrid CCO + Statistical Subdivision:**

Growing millions of capillaries one-at-a-time with CCO is O(n²) and would take weeks. Instead:

1. **CCO skeleton** (~30s): Grow 500-2000 terminals per artery at upper Strahler orders (5-11) using grid-accelerated nearest-neighbor search. This produces the realistic spatial layout of major vessels with proper intersection checking.

2. **Statistical subdivision** (~2 min): For each CCO terminal of order k, recursively generate daughters using Kassab's connectivity matrix CM[m+1, k+1]. Each daughter order m gets diameter ~ N(mean_m, sd_m), length ~ N(mean_m, sd_m). NO intersection checking — just tree construction. This multiplies ~1500 skeleton segments into ~2M+ per artery.

3. **Post-hoc refinement** (~30s): Apply Kassab asymmetry radii top-down. Apply Barabasi junction geometry. Enforce Murray's law bottom-up.

**Why this works**: The connectivity matrix IS the Kassab morphometric model. It directly encodes how many daughters of each order every parent has. Combined with per-order diameter/length distributions, it statistically reproduces the full tree without physics-based growth at capillary scale.

**Target**: 3 arteries × ~2M segments = ~6M total segments in < 5 minutes.

## P10: Exact Kassab Parity — The Current Priority

P9 achieved the basic pipeline (CCO + subdivision + refinement) but validation shows only 2/6 metrics pass. P10 fixes this with REAL data and correct methodology.

**Root causes of P9 failures (documented in ralph_loop/research/kassab_parity_gaps.md):**
1. **Fabricated CM values** — parameters.jl has made-up numbers, not real Kassab Tables 6-8
2. **Wrong Strahler method** — simple diameter binning instead of iterative Jiang 1994 algorithm
3. **No element grouping** — Kassab operates on elements, not segments
4. **Cascade stubs inflate trifurcation rate** — 44% vs target 8.3%
5. **Post-hoc asymmetry gets destroyed** — Murray propagation + floor clamp
6. **Wrong diameter/length data** — fabricated instead of real Tables 1-3
7. **Single CM for all arteries** — RCA/LAD/LCX have different CMs

**P10 stories (VESSEL-1034 through VESSEL-1043) fix ALL of these:**
1. Verify real data from PDFs (VESSEL-1034)
2. Replace fabricated CM + add per-artery params (VESSEL-1035)
3. Implement correct Strahler ordering (VESSEL-1036)
4. Implement element grouping (VESSEL-1037)
5. Fix subdivision to produce clean bifurcation trees (VESSEL-1038)
6. Fix asymmetry distribution (VESSEL-1039)
7. Update validation to element-level (VESSEL-1040)
8. Update pipeline with new params (VESSEL-1041)
9. Verify >= 7/9 metrics pass (VESSEL-1042)
10. Polish to >= 8/9 (VESSEL-1043)

**Target**: >= 8/9 validation metrics pass for all three arteries (RCA, LAD, LCX), using exclusively real published data with zero fabrication.

## Output Markers

- `RALPH_COMPLETE` — All stories done, loop can exit
- `RALPH_BLOCKED` — Cannot proceed, needs human intervention
