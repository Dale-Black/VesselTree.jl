# Kassab Parity — Gap Analysis

Current validation: 2/6 metrics pass (33%). This document identifies every gap
between the current implementation and exact Kassab parity, with root causes and fixes.

## Gap 1: FABRICATED Connectivity Matrix

**Problem**: `parameters.jl` contains hand-fabricated CM values labeled "approximate values".
These do not match ANY table in Kassab 1993.

**Evidence**: Current CM has simple patterns like CM[1,2]=2.3, CM[2,3]=2.1 — real Kassab
data has CM[1,2]=2.75 for RCA, 3.18 for LAD/LCX. Off-diagonal entries are completely wrong.

**Fix**: Replace with real Tables 6 (RCA), 7 (LAD), 8 (LCX) from Kassab 1993.
Need three separate parameter presets: `kassab_rca_params()`, `kassab_lad_params()`, `kassab_lcx_params()`.

**Impact**: Connectivity chi-squared test, segment count ratios, overall tree structure.

## Gap 2: Wrong Strahler Ordering Method

**Problem**: Current `assign_strahler_orders!` in `kassab.jl` uses simple diameter binning with
fixed thresholds from `params.diameter_bounds`. This is NOT Kassab's method.

**Kassab's method** (Jiang 1994): Iterative diameter-defined Strahler ordering:
1. Start from topological Strahler
2. Compute per-order mean/SD
3. Redefine bounds as midpoints between (mean+SD) and (mean-SD) of adjacent orders
4. Reassign, iterate until <1% change

**Fix**: Implement `assign_diameter_defined_strahler!` with convergence loop.

**Impact**: Every metric that depends on order classification.

## Gap 3: Asymmetry Distribution Wrong

**Problem**: Current median asymmetry is 0.94 (should be ~0.76). Post-hoc asymmetry from
`apply_full_kassab_radii!` gets washed out by Murray's law bottom-up propagation + floor clamping.

**Root causes**:
1. Floor clamping at 8um pushes small daughters up, increasing asymmetry ratio
2. Bottom-up Murray propagation overwrites the top-down asymmetry assignments
3. Beta(2.5, 0.8) theoretical median is ~0.81, not 0.76

**Fix**:
- Apply asymmetry DURING subdivision (not post-hoc) so it's baked into the tree structure
- Remove or loosen floor clamping (let Murray's law set radii naturally)
- Verify Beta parameters against Kassab's actual asymmetry data
- Consider using element-level asymmetry (not segment-level)

**Impact**: Asymmetry KS test, diameter distributions, overall tree realism.

## Gap 4: Trifurcation Rate Way Too High (44% vs ~8.3%)

**Problem**: 44% of junctions are trifurcations. Kassab/Barabasi data shows ~8.3%.

**Root cause**: Subdivision cascade stubs. When Poisson sampling gives >3 daughters for a
parent, the code creates "cascade stubs" — tiny segments that bridge to additional child slots.
These stubs get classified as trifurcations because they have 3 children.

**Fix**:
- Eliminate cascade stubs: limit daughters to 2 per parent (pure bifurcation tree)
- If CM says 5 daughters, create a chain of bifurcations (each splits off one daughter)
- Only allow real trifurcations when Barabasi chi > 0.83 criterion is met

**Impact**: Trifurcation percentage, junction angle distributions.

## Gap 5: Length Distributions Completely Wrong

**Problem**: 0/10 orders pass KS test for length.

**Root causes**:
1. CCO-generated segments have lengths determined by spatial layout (not Kassab distributions)
2. Subdivision uses per-order Normal(mean, sd) but the mean/sd values in parameters.jl
   are fabricated (not from real Kassab Table 1/2/3)
3. No distinction between segment length and element length
4. CCO stub segments (very short) pollute length distributions

**Fix**:
- Use real Kassab segment/element length data from Tables 1-3
- During subdivision, sample from correct per-order distributions
- For CCO skeleton segments, either: (a) adjust lengths post-hoc, or (b) use them as-is
  since they represent spatial layout and validate only subdivision segments
- Implement element grouping and validate element lengths (not segment lengths)

**Impact**: Length KS tests, L/D ratios.

## Gap 6: Diameter Distributions Mostly Wrong

**Problem**: Only 1-3/10 orders pass KS test for diameter.

**Root causes**:
1. `apply_full_kassab_radii!` shifts diameters away from per-order N(mean, sd)
2. Floor clamping at 8um creates a delta function at order 0 instead of N(8, 2)
3. Diameter mean/sd values in parameters.jl may not match real Kassab data exactly
4. Murray's law bottom-up propagation changes all radii

**Fix**:
- Use real Kassab segment diameter data from Tables 1-3
- Assign diameters during subdivision from correct distributions
- Do NOT overwrite with post-hoc radius refinement (or make it much gentler)
- Remove hard floor clamp; use soft distribution-based minimum

**Impact**: Diameter KS tests, order classification.

## Gap 7: No Per-Artery Parameters

**Problem**: Single `kassab_coronary_params()` used for all three arteries.
Kassab 1993 shows RCA, LAD, and LCX have different CMs, different N_elements,
different diameter/length distributions, and different number of orders.

**Fix**: Create `kassab_rca_params()`, `kassab_lad_params()`, `kassab_lcx_params()`.
LCX has only 10 orders (not 12). Update `coronary_tree_configs` to use artery-specific params.

**Impact**: All per-artery validation metrics.

## Gap 8: No Element Grouping

**Problem**: No concept of "element" in the code. All validation operates on segments.
Kassab's CM and statistics are element-based.

**Fix**: Implement `group_into_elements(tree, params)` that:
1. Assigns diameter-defined Strahler orders
2. Walks the tree and groups consecutive same-order segments into elements
3. Returns element-level statistics (diameter = mean of segments, length = sum of segments)
4. CM is built from elements

**Impact**: All Kassab validation metrics.

## Priority Order for Fixes

1. **Replace CM with real data** (Gap 1) — highest impact, easiest fix
2. **Per-artery parameters** (Gap 7) — needed for accurate CM
3. **Implement element grouping** (Gap 8) — needed for correct validation
4. **Fix Strahler ordering** (Gap 2) — needed for correct order assignment
5. **Fix trifurcation rate** (Gap 4) — structural fix to subdivision
6. **Fix diameter distributions** (Gap 6) — bake into subdivision
7. **Fix length distributions** (Gap 5) — bake into subdivision
8. **Fix asymmetry** (Gap 3) — bake into subdivision
