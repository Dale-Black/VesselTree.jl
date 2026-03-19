### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000001
begin
    import Pkg
    Pkg.activate(dirname(@__DIR__))
    Pkg.instantiate()
end

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000002
using VesselTree

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000003
md"""
# XCAT CT Time-Series Viewer

This notebook loads precomputed GE scan outputs and exports a browser viewer with:

- a **time slider** for `t = 0..9 s`,
- a **slice slider**,
- axial / coronal / sagittal plane switching,
- a lightweight local HTTP server for browser inspection.
"""

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000004
_env_bool(name::AbstractString, default::Bool) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000005
begin
    default_run_manifest = joinpath(dirname(@__DIR__), "output", "ge_scans_hu_all_full", "20260318T160257", "ge_scan_run.json")
    run_manifest_path = get(ENV, "CT_VIEWER_RUN_MANIFEST", default_run_manifest)
    output_dir = get(ENV, "CT_VIEWER_OUTPUT_DIR", joinpath(dirname(@__DIR__), "output", "ct_viewer"))
    host = get(ENV, "CT_VIEWER_HOST", "127.0.0.1")
    port = parse(Int, get(ENV, "CT_VIEWER_PORT", "8008"))
    prefer_hu = _env_bool("CT_VIEWER_PREFER_HU", true)
    q_low = parse(Float64, get(ENV, "CT_VIEWER_Q_LOW", "0.01"))
    q_high = parse(Float64, get(ENV, "CT_VIEWER_Q_HIGH", "0.995"))
    serve_viewer = _env_bool("CT_VIEWER_SERVE", true)
    block_server = _env_bool("CT_VIEWER_BLOCK", true)
end

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000006
viewer = load_ct_viewer_data(
    run_manifest_path;
    prefer_hu=prefer_hu,
    q_low=q_low,
    q_high=q_high,
    title="XCAT CT Time-Series Viewer",
)

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000007
html_path = export_ct_viewer_html(joinpath(output_dir, "index.html"), viewer; title=viewer.title)

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000008
server = if serve_viewer
    serve_static_directory(output_dir; host=host, port=port)
else
    nothing
end

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000009
md"""
## Status

- run manifest: `$(run_manifest_path)`
- HTML: `$(html_path)`
- frames: `$(length(viewer.times))`
- unit: `$(viewer.value_unit)`
- display range: `$(viewer.display_min)` to `$(viewer.display_max)`
- browser URL: `http://$(host):$(port)/`
"""

# ╔═╡ c1a90036-3c5d-4d0a-a001-100000000010
if serve_viewer && block_server
    println("CT viewer HTML written to: ", html_path)
    println("Open this in your browser: http://", host, ":", port, "/")
    println("Server is running. Press Ctrl+C to stop.")
    while true
        sleep(3600.0)
    end
elseif serve_viewer && server !== nothing
    sleep(parse(Float64, get(ENV, "CT_VIEWER_SERVE_SECONDS", "2.0")))
    close(server)
end

# ╔═╡ Cell order:
# ╠═c1a90036-3c5d-4d0a-a001-100000000001
# ╠═c1a90036-3c5d-4d0a-a001-100000000002
# ╠═c1a90036-3c5d-4d0a-a001-100000000003
# ╠═c1a90036-3c5d-4d0a-a001-100000000004
# ╠═c1a90036-3c5d-4d0a-a001-100000000005
# ╠═c1a90036-3c5d-4d0a-a001-100000000006
# ╠═c1a90036-3c5d-4d0a-a001-100000000007
# ╠═c1a90036-3c5d-4d0a-a001-100000000008
# ╠═c1a90036-3c5d-4d0a-a001-100000000009
# ╠═c1a90036-3c5d-4d0a-a001-100000000010
