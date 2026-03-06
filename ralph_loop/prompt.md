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
2. **Kassab et al. 1993** — Coronary morphometry (Am J Physiol Heart)
3. **Huo & Kassab 2007** — gamma = 7/3 derivation (Biophys J)
4. **Barabasi et al. 2026** — Surface optimization (Nature 649:315-322)
5. **OpenCCO** — https://github.com/OpenCCO-team/OpenCCO (C++ reference)
6. **min-surf-netw** — https://github.com/Barabasi-Lab/min-surf-netw (Python reference)

## Checkpoint Discipline

Commit after:
- Every test file that passes
- Every feature implementation that works
- Every 15 minutes of work (even if incomplete — WIP commits are fine)
- BEFORE running long commands

Never lose more than 15 minutes of work to a timeout.

## Output Markers

- `RALPH_COMPLETE` — All stories done, loop can exit
- `RALPH_BLOCKED` — Cannot proceed, needs human intervention
