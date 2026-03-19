### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 9ad4ac20-b5ad-4b0a-9d44-000000000001
md"""
# XCAT GE Recon -> DICOM Export

This notebook exports already reconstructed GE scan volumes to DICOM series **without rerunning the simulator**.

It wraps `scripts/export_ct_dicoms.py`, which reads the existing `ge_scan_run.json` plus the corresponding `*_scan.jld2` files and emits one DICOM series per time point.
"""

# ╔═╡ 9ad4ac20-b5ad-4b0a-9d44-000000000002
begin
    const DEFAULT_RUN_MANIFEST = joinpath(dirname(@__DIR__), "output", "ge_scans_hu_all_full", "20260318T160257", "ge_scan_run.json")
    const DEFAULT_OUTPUT_DIR = joinpath(dirname(@__DIR__), "output", "dicom_exports")
    const DEFAULT_PYTHON = normpath(joinpath(dirname(@__DIR__), "..", ".conda", "bin", "python"))
    const DEFAULT_SCRIPT = joinpath(dirname(@__DIR__), "scripts", "export_ct_dicoms.py")
    const ENV_TRUE = ("1", "true", "yes")
end

# ╔═╡ 9ad4ac20-b5ad-4b0a-9d44-000000000003
begin
    run_export = lowercase(get(ENV, "XCAT_RUN_DICOM_EXPORT", "false")) in ENV_TRUE
    run_manifest_path = get(ENV, "XCAT_DICOM_RUN_MANIFEST", DEFAULT_RUN_MANIFEST)
    output_dir = get(ENV, "XCAT_DICOM_OUTPUT_DIR", DEFAULT_OUTPUT_DIR)
    python_exe = get(ENV, "XCAT_DICOM_PYTHON", DEFAULT_PYTHON)
    script_path = get(ENV, "XCAT_DICOM_SCRIPT", DEFAULT_SCRIPT)
    overwrite = lowercase(get(ENV, "XCAT_DICOM_OVERWRITE", "true")) in ENV_TRUE
end

# ╔═╡ 9ad4ac20-b5ad-4b0a-9d44-000000000004
export_cmd = `$(python_exe) $(script_path) --run-manifest $(run_manifest_path) --output-dir $(output_dir) $(overwrite ? "--overwrite" : "")`

# ╔═╡ 9ad4ac20-b5ad-4b0a-9d44-000000000005
run_export ? run(export_cmd) : nothing

# ╔═╡ 9ad4ac20-b5ad-4b0a-9d44-000000000006
md"""
## Run from the terminal

```bash
XCAT_RUN_DICOM_EXPORT=true \
XCAT_DICOM_RUN_MANIFEST="$(run_manifest_path)" \
XCAT_DICOM_OUTPUT_DIR="$(output_dir)" \
$(Base.julia_cmd()) --project="$(dirname(@__DIR__))" \
  -e 'include("$(joinpath(dirname(@__DIR__), "notebooks", "xcat_export_dicom.jl"))")'
```
"""
