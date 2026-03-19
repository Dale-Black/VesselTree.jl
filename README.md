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

### 9. Segment-level contrast transport

Contrast transport is now available as a separate library layer on top of the
tree hemodynamics.

Main entry points:

- `gamma_variate_input(...)`
- `root_pulse_input(...)`
- `simulate_contrast_transport(...)`
- `simulate_forest_contrast(...)`
- `segment_mass_mg(...)`

Reference examples:

- [xcat_contrast_transport.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_contrast_transport.jl)
- [xcat_contrast_viewer.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_contrast_viewer.jl)

Browser viewer notes:

- the viewer exports a self-contained `index.html` into `v4/output/contrast_viewer/`,
- it can also start a local HTTP server on `127.0.0.1:8008` by default,
- the browser UI includes a draggable time slider plus Play/Pause controls,
- the default viewer now shows every `0.05 s` time step,
- the default sampling keeps more distal/small vessels visible,
- gray lines show the vascular skeleton and colored markers show the segment-average iodine concentration over time,
- the default color map is tuned for low-concentration visibility: near `0` is gray, around `3 mg/mL` is blue, and near the peak is red.

Launch command:

```bash
JULIA_DEPOT_PATH=/tmp/julia_depot:/home/molloi-lab/.julia \
/home/molloi-lab/.julia/juliaup/julia-1.10.10+0.x64.linux.gnu/bin/julia \
  --project="v4" v4/examples/xcat_contrast_viewer.jl
```

Then open `http://127.0.0.1:8008/` in a browser.

Current model assumptions:

- each vascular segment is treated as a well-mixed control volume,
- segment transit time is computed from `tau = V / Q`,
- outlet concentration is the inlet concentration delayed by the segment transit time,
- segment-average concentration is updated by a mass-balance ODE.

This is a reduced transport model driven by the tree flow solution; it is not a
full 3D Navier-Stokes simulation.

### 10. Unified browser visualization

The recommended browser viewer is now:

- [xcat_unified_viewer.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_unified_viewer.jl)

This exports a single interactive `index.html` into `v4/output/unified_viewer/`
and can also launch a local HTTP server on `127.0.0.1:8008`.

The unified viewer combines these layers in one scene, each with a toggle:

- domain surface
- cardiac cavity points
- sampled original XCAT vessel surfaces
- imported fixed trunks
- grown vascular tree
- time-varying iodine concentration

Viewer notes:

- the iodine layer uses the same low-concentration-enhanced color map as the standalone contrast viewer,
- the default time step is `0.05 s`,
- hover on iodine points shows concentration, segment length, and segment diameter,
- geometry layers are intentionally decimated to keep the browser responsive,
- the page is intended for SSH workflows, so port-forwarding `8008` is usually required.

Launch command:

```bash
JULIA_DEPOT_PATH=/tmp/julia_depot:/home/molloi-lab/.julia /home/molloi-lab/.julia/juliaup/julia-1.10.10+0.x64.linux.gnu/bin/julia   --project="v4" v4/examples/xcat_unified_viewer.jl
```

Then open `http://127.0.0.1:8008/` in a browser.

The current XCAT example is wired to the latest visually validated coronary setup, including:

- full-cutoff subdivision to `8 um`,
- per-tree proximal fixed-prefix blending,
- a single browser scene for geometry plus iodine flow inspection.

### 11. Save the current run to disk

Use:

- [xcat_save_artifacts.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/examples/xcat_save_artifacts.jl)

This script saves three groups of artifacts for the current XCAT coronary run:

- timestamped tree snapshots for each branch, using `BRANCH-TIMESTAMP` filenames,
- dynamic contrast data, including per-segment time courses plus a voxelized iodine concentration volume,
- a fused `.nrb` model that appends the grown vascular tree back into the original XCAT file.

Default output location:

- `v4/output/saved_runs/<TIMESTAMP>/`

Directory layout:

- `trees/`
  Contains branch-specific snapshots such as `LAD-<TIMESTAMP>.jld2`, `LCX-<TIMESTAMP>.csv`, and Wenbo-format `.txt` exports.
- `contrast/`
  Contains `contrast-segments-<TIMESTAMP>.jld2` and `contrast-volume-<TIMESTAMP>.jld2`.
  The volume export always stores a sparse voxel representation and will also save a dense `4D` array automatically when the estimated array size is reasonable.
- `model/`
  Contains `xcat-grown-<TIMESTAMP>.nrb`, which preserves the original XCAT surfaces and appends the grown tree as additional vessel surfaces.

Useful environment variables:

- `XCAT_ARTIFACT_OUTPUT_DIR`
- `XCAT_CONTRAST_DENSE_MODE=auto|always|never`
- `XCAT_CONTRAST_VOXEL_SPACING_MM`
- `XCAT_CONTRAST_VOXEL_SUPERSAMPLE`

Launch command:

```bash
JULIA_DEPOT_PATH=/tmp/julia_depot:/home/molloi-lab/.julia \
/home/molloi-lab/.julia/juliaup/julia-1.10.10+0.x64.linux.gnu/bin/julia \
  --project="v4" v4/examples/xcat_save_artifacts.jl
```

### 12. Compile `basis_simulator` phantom frames

Use the new notebook:

- [xcat_basis_frames.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/notebooks/xcat_basis_frames.jl)

This notebook does **not** modify `basis_simulator`. Instead, it compiles one full phantom frame per requested time point by combining:

- the baseline XCAT raw label volume,
- the non-contrast XCAT material spreadsheet,
- the newly grown distal coronary tree,
- the time-varying segment iodine concentrations.

Current behavior:

- writes one frame per requested time point (for example `t = 0:9 s`),
- keeps the original XCAT raw grid as the simulator grid,
- overlays only the newly grown distal tree rather than rewriting the original proximal XCAT vessels,
- emits a full label mask plus a `label => XA.Material` dictionary for each frame, ready for downstream `basis_simulator` use.

The notebook defaults to `run_pipeline = false` so it is safe to open without immediately launching a heavy XCAT regeneration run.

Core bridge functions now live in `VesselTree`:

- `load_xcat_raw_labels(...)`
- `load_xcat_materials_from_xlsx(...)`
- `estimate_xcat_raw_alignment(...)`
- `compile_basis_material_frame(...)`
- `save_basis_frame(...)`
- `generate_basis_frames(...)`

### 13. Run GE scanner smoke scans on compiled frames

Use the new notebook:

- [xcat_basis_ge_scan.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/notebooks/xcat_basis_ge_scan.jl)

This notebook stays outside the simulator core and simply **applies** the existing `BasisSimulator` GE workflow to one or more precompiled basis frames. It is intended as the bridge from:

- `xcat_basis_frames.jl` output
- to GE-style forward projection + FDK reconstruction

Current behavior:

- loads one or more frame manifests from `xcat_basis_frames.json`,
- crops around the dynamic coronary overlay,
- optionally downsamples the cropped mask for faster smoke testing,
- builds a `BS.Phantom(mask, materials_dict, voxel_size_cm)`,
- runs the existing GE scanner simulation and FDK reconstruction,
- saves one scan artifact set per frame.

Recommended smoke-scan controls are exposed through environment variables:

- `XCAT_RUN_GE_SCAN=true`
- `XCAT_SCAN_ALL_FRAMES=true|false`
- `XCAT_FRAME_INDEX`
- `XCAT_GE_DOWNSAMPLE_FACTOR`
- `XCAT_GE_RECON_XY_CAP`
- `XCAT_GE_RECON_SLICES_CAP`
- `XCAT_GE_CALIBRATE_WATER=true|false`
- `XCAT_GE_PROTOCOL_VIEWS`
- `XCAT_GE_PROTOCOL_MA`
- `XCAT_GE_SIM_FIDELITY=high|medium|low`
- `XCAT_GE_SCAN_OUTPUT_DIR`

This notebook requires **Julia 1.11+** because `basis_simulator` is pinned to that range. The default `VesselTree` development workflow can stay on Julia 1.10; only the GE scan notebook needs the newer runtime.

### 14. View reconstructed CT time series in the browser

Use the new notebook:

- [xcat_ct_viewer.jl](/media/molloi-lab/2TB3/wenbo%20playground/flow%20simulation%20tree%20generation/v4/notebooks/xcat_ct_viewer.jl)

This notebook loads one `ge_scan_run.json` manifest and exports a lightweight browser viewer for the reconstructed CT volumes.

Current viewer features:

- a **time slider** across the reconstructed frames (for example `t = 0..9 s`),
- a **slice slider**,
- axial / coronal / sagittal plane switching,
- Play / Pause controls for time playback,
- simple display-window presets based on the reconstructed intensity range.

By default it points at the latest smoke-scan run:

- `v4/output/ge_scans_hu_all_full/20260318T160257/ge_scan_run.json`

Useful environment variables:

- `CT_VIEWER_RUN_MANIFEST`
- `CT_VIEWER_OUTPUT_DIR`
- `CT_VIEWER_HOST`
- `CT_VIEWER_PORT`
- `CT_VIEWER_PREFER_HU=true|false`
- `CT_VIEWER_Q_LOW`
- `CT_VIEWER_Q_HIGH`
- `CT_VIEWER_SERVE=true|false`

The notebook reuses the existing static HTTP helper and does not rerun any simulation; it only visualizes previously saved scan outputs.

The current recommended default is the full HU-calibrated 10-frame run (`t = 0..9 s`) at `424 x 424 x 128`.

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
- Contrast transport currently assumes a directed arterial tree without local recirculation or axial dispersion inside each segment; the first implementation tracks segment-average concentration rather than within-segment concentration profiles.

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
7. Run the contrast-transport example if you need dynamic iodine concentration.
8. Run the unified viewer example if you need one browser scene for domain, vessels, and iodine flow.
9. Use the calibration script when tuning full-cutoff targets.
10. Only then move to the heaviest full-cutoff runs.

This sequence is much faster for debugging than jumping directly to the heaviest workflow.


### 15. Export reconstructed CT frames to DICOM

You can convert any completed `ge_scan_run.json` plus its `*_scan.jld2` files into DICOM series **without rerunning the scanner**.

Notebook entry point:

- `v4/notebooks/xcat_export_dicom.jl`

Python exporter used by the notebook:

- `v4/scripts/export_ct_dicoms.py`

Example:

```bash
XCAT_RUN_DICOM_EXPORT=true \
XCAT_DICOM_RUN_MANIFEST="/path/to/ge_scan_run.json" \
XCAT_DICOM_OUTPUT_DIR="/path/to/dicom_export" \
julia --project="/path/to/v4" \
  -e 'include("/path/to/v4/notebooks/xcat_export_dicom.jl")'
```

The exporter writes:

- one DICOM series per time point (`t000`, `t001`, ...)
- one axial slice per `.dcm` file
- a `dicom_export_manifest.json` summary in the export directory

