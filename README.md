# VesselTree v4.0

## Overview

This repository contains the current Julia implementation of a vascular tree growth pipeline for:

- importing an existing proximal vasculature from an ASCII `.nrb` file,
- building an organ domain from NURBS surfaces,
- continuing growth from the imported vessels into the distal microvascular bed,
- optionally applying Kassab-style statistical subdivision down to cutoff scale,
- exporting trees for downstream flow simulation.

The current production path is the XCAT heart workflow, but the code now includes a general `NRB` interface so the same pattern can be reused for other organs such as:

- lower leg,
- thigh,
- brain white matter,
- any other organ with an `.nrb` surface model and existing proximal vessels.

The intended use case is:

1. start from an organ `.nrb` file,
2. import the existing vessel trunks,
3. continue distal growth inside the organ domain,
4. optionally resolve all remaining orders down to `8 um` terminal diameter.

## Repository Layout

- `src/`
  Core library code.
- `examples/`
  Runnable scripts for inspection, generation, continuation, hemodynamics, and full-cutoff runs.
- `output/`
  Generated figures, CSVs, and exported Wenbo-format text files.
- `model_CSVs/`
  Legacy CSV domain inputs used by the earlier coronary volume workflow.
- `notebooks/`
  Existing notebook-based workflows.

## Main Entry Points

### Generic `.nrb` interfaces

These are the interfaces for future organs:

- `NRBTreeSpec`
- `NRBOrganSpec`
- `nrb_shell_domain(...)`
- `nrb_import_fixed_trees(...)`
- `generate_nrb_continuation_forest(...)`
- `generate_nrb_kassab_forest(...)`

These live in:

- [nrb_growth.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/src/nrb_growth.jl)

### XCAT heart helpers

The heart-specific wrappers remain available for convenience:

- `xcat_heart_organ_spec(...)`
- `xcat_import_coronary_trees(...)`
- `generate_xcat_coronary_forest(...)`
- `generate_xcat_kassab_coronary(...)`

These are thin wrappers over the generic `NRB` workflow.

## Recommended Workflows

### 1. Inspect a new `.nrb` file

Use:

- [xcat_nrb_inspect.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_nrb_inspect.jl)

This confirms:

- object names,
- relative placement,
- control-net sanity.

### 2. Inspect the organ shell domain

Use:

- [xcat_shell_domain_inspect.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_shell_domain_inspect.jl)

This verifies:

- outer surface,
- cavity exclusions,
- interior point sampling.

### 3. Inspect imported proximal trees

Use:

- [xcat_coronary_centerlines_inspect.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_coronary_centerlines_inspect.jl)
- [xcat_import_fixed_trees.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_import_fixed_trees.jl)

These verify:

- centerlines,
- branch connectivity,
- imported `VascularTree` topology.

### 4. Continue growth from imported trunks

Use:

- [xcat_continue_forest_demo.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_continue_forest_demo.jl)

This is the lightweight continuation-only demo.

### 5. Generate a deeper XCAT coronary forest

Use:

- [xcat_kassab_coronary.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_kassab_coronary.jl)

This runs:

- imported fixed trunks,
- continuation growth,
- Kassab subdivision,
- optional Barabasi geometry.

### 6. Full-cutoff generation to `8 um`

Use:

- [xcat_kassab_full_cutoff.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_kassab_full_cutoff.jl)

This is the current reference path for forcing all terminals to:

- `Strahler order 0`
- `8 um` terminal diameter

Important:

- this mode is computationally heavy,
- geometry application is usually disabled during smoke tests,
- segment counts can quickly reach hundreds of thousands to millions per tree.

### 7. Full-cutoff calibration against coronary targets

Use:

- [xcat_calibrate_full_cutoff.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_calibrate_full_cutoff.jl)

This script is intended for parameter sweeps and prints:

- configured continuation targets,
- actual in-domain territory fractions,
- proximal imported-root radius diagnostics,
- baseline and hyperemic flow,
- a simple absolute-error score against the published hyperemic targets.

It reads its main knobs from environment variables, including:

- `LAD_TARGET`
- `LCX_TARGET`
- `RCA_TARGET`
- `XCAT_MAX_TREE_CAPACITY`
- `XCAT_APPLY_GEOMETRY`

### 8. Hemodynamics and Wenbo export

Use:

- [xcat_kassab_hemodynamics.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_kassab_hemodynamics.jl)
- [xcat_kassab_full_cutoff_hemodynamics.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_kassab_full_cutoff_hemodynamics.jl)

These scripts:

- compute baseline flow,
- compute vasodilated flow,
- export Wenbo-format `.txt` trees when appropriate.

## Important Modeling Notes

### Pressure conditions

For explicit arterial trees, the flow comparison uses:

- inlet pressure: `100 mmHg`
- outlet pressure: `15 mmHg`

The `15 mmHg` outlet corresponds to pre-capillary pressure, not venous pressure.

### Vasodilation vs baseline

For the coronary paper comparisons:

- the published `242 / 116 / 214 mL/min` numbers correspond to the **vasodilated / hyperemic** state,
- not the baseline state.

This means:

- baseline flow should be lower,
- hyperemic flow is the correct quantity to compare against the paper target values.

### Full-cutoff mode

`enforce_full_cutoff=true` is the current mechanism for forcing the subdivision process to continue until all terminals reach cutoff scale.

Without it, some branches can terminate prematurely at higher order because the raw connectivity-matrix draw may produce too few daughters for a full terminal bifurcation chain.

## Using the Generic `.nrb` Interface for Other Organs

The simplest starting point is:

- [nrb_generic_template.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/nrb_generic_template.jl)

To adapt it:

1. Identify the outer organ surface name.
2. Identify any cavity or exclusion surfaces.
3. Identify the proximal vessel surfaces that should be imported.
4. Group those vessel surfaces into one or more `NRBTreeSpec`s.
5. Provide a reference surface name if root orientation depends on a parent conduit.

Minimal example:

```julia
spec = NRBOrganSpec(
    "Example Organ",
    "example_outer_surface",
    ["example_cavity_1", "example_cavity_2"],
    [
        NRBTreeSpec(
            "TREE_A",
            ["tree_a_root", "tree_a_mid", "tree_a_distal"];
            target_terminals=40,
            territory_fraction=0.6,
        ),
        NRBTreeSpec(
            "TREE_B",
            ["tree_b_root", "tree_b_distal"];
            target_terminals=25,
            territory_fraction=0.4,
        ),
    ],
    "example_reference_surface",
)
```

Then run:

```julia
result = generate_nrb_kassab_forest(
    nrb_path,
    spec,
    params;
    handoff_order=6,
    subdivision_max_order=6,
)
```

## Capacity and Performance

Deep subdivision can be expensive.

The main controls are:

- `target_terminals`
  Controls the continuation skeleton size.
- `handoff_order`
  Controls how deep the usual subdivision goes.
- `subdivision_max_order`
  Allows explicit control of subdivision depth independent of the continuation skeleton.
- `enforce_full_cutoff`
  Forces cutoff-scale resolution.
- `apply_geometry`
  Can be disabled during large smoke tests to avoid spending most runtime on final geometry updates.
- `max_tree_capacity`
  Hard cap on per-tree segment allocation.

If a run fails with `SegmentData at capacity`, increase:

- `max_tree_capacity`

or reduce:

- `target_terminals`
- `subdivision_max_order`

## Current Known Limitations

- The generic `.nrb` interface assumes vessel surface groups can be represented as ordered surface chains plus explicit parent-child connections.
- Full-cutoff runs can become memory-heavy very quickly.
- The built-in hemodynamics solver is useful for trend-checking, but the Wenbo export path remains important for comparison against the existing lab-specific flow workflow.
- Large full-cutoff runs may need `apply_geometry=false` during exploration, then a separate geometry pass later if required.
- Imported proximal coronary roots from XCAT still show evidence of cap-artifact bias at the aorta-facing end; this can depress inlet radius and remains an active source of flow mismatch during calibration.
- Current XCAT calibration is closest for LAD/RCA; LCX remains the hardest branch to match to the published hyperemic target without pushing the other trees off target.

## Quick Start Commands

From the `v4` project directory:

```bash
julia --project=. examples/xcat_kassab_coronary.jl
```

```bash
julia --project=. examples/xcat_kassab_hemodynamics.jl
```

```bash
julia --project=. examples/xcat_kassab_full_cutoff.jl
```

## Recommended Onboarding Order for New Lab Members

1. Run the `.nrb` inspection example.
2. Run the shell-domain inspection example.
3. Run the fixed-tree import example.
4. Run the continuation demo.
5. Run the Kassab coronary example.
6. Run the hemodynamics example.
7. Use the calibration script when tuning full-cutoff targets.
8. Only then move to the heaviest full-cutoff runs.

This sequence is much faster for debugging than jumping directly to the heaviest workflow.
