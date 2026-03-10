### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ d64a065b-258a-47df-86e4-eba5c8966829
begin
	import Pkg
	Pkg.activate(dirname(@__DIR__))
end

# ╔═╡ 38d9074f-a4ad-4ae0-9c22-4a57a9ecfaf3
# ╠═╡ show_logs = false
using VesselTree

# ╔═╡ 6f7c71c6-a2c4-4e69-af63-3cc382df760d
using CairoMakie

# ╔═╡ 86de57c5-bb12-425f-a2d3-218f805f6761
using Random

# ╔═╡ e5c54f95-2461-46e7-9870-d411df930300
using Statistics

# ╔═╡ e4965a7a-0800-49fb-93b9-be63caf2111f
md"""
# Coronary Arterial Tree Generation — VesselTree.jl

Pure CCO (Constrained Constructive Optimization) growth matching the svVascularize approach. Murray's law radii (γ = 7/3). Domain is an ellipsoid shell modeling the myocardial wall.
"""

# ╔═╡ 87657f8c-3c2b-4f06-b89d-b77de7bdae8b
md"""
## Step 1: Define the Heart Domain

The domain is an **ellipsoid shell** modeling the myocardial wall:
- **a = 50mm** (left-right), **b = 35mm** (anterior-posterior), **c = 45mm** (base-apex)
- **thickness = 30%** → wall ~10-15mm, matching real myocardium
"""

# ╔═╡ bd7d2e5c-f5ea-4dad-a131-0b61d4b99de6
domain = default_coronary_domain()  # EllipsoidShellDomain

# ╔═╡ 7f5bdb8c-d128-4324-9b57-db15adac6de3
md"""
**Domain:** `EllipsoidShellDomain`
- Center: $(domain.center)
- Semi-axes: $(domain.semi_axes) mm
- Shell thickness: $(domain.thickness) ($(round(domain.thickness * 100))%)
- Wall thickness: ~$(round(minimum(domain.semi_axes) * domain.thickness, digits=1)) – $(round(maximum(domain.semi_axes) * domain.thickness, digits=1)) mm
"""

# ╔═╡ 28d33479-2192-4cc2-8e2f-d1dd2772fb69
md"""
## Step 2: Define Root Positions (Seed Points)

| Artery | Territory | Anatomical Position |
|:-------|:----------|:-------------------|
| **LAD** | Anterior wall (40%) | Anterior, slightly left |
| **LCX** | Left/posterior wall (25%) | Left-posterior |
| **RCA** | Right wall (35%) | Right-anterior |
"""

# ╔═╡ b26965de-6069-4c3d-a9c9-8be15d74bf10
begin
	configs = VesselTree.coronary_tree_configs(domain)

	config_rows = join([
		"| $(cfg.name) | ($(round.(cfg.root_position, digits=1))) | $(cfg.root_radius) mm | $(cfg.target_terminals) | $(round(cfg.territory_fraction * 100))% |"
		for cfg in configs
	], "\n")
end;

# ╔═╡ 52055f05-592d-473a-88a8-eb4a5034065d
Markdown.parse("""
| Artery | Root Position (mm) | Root Radius | Target Terminals | Territory |
|:-------|:-------------------|:------------|:-----------------|:----------|
$config_rows
""")

# ╔═╡ 6df19e82-fa9f-4e26-bc54-7199d330abb8
md"""
## Step 3: Growth Parameters

Murray's law exponent γ = 7/3 (Huo-Kassab 2007). Per-artery Kassab morphometric params provide terminal radii for CCO handoff.
"""

# ╔═╡ 33573674-ea88-406b-b3c6-9013b976ee37
params_rca = kassab_rca_params();

# ╔═╡ 0a987a32-e3d5-44c5-b6a9-d2cf0b29bde1
md"""
γ = $(params_rca.gamma), $(params_rca.n_orders) Strahler orders
"""

# ╔═╡ acb6d732-8524-4e60-85d3-a5bc4fa56d86
md"""
## Step 4: CCO Progressive Growth (LAD only)

Watch the LAD grow progressively to verify it wraps around the shell correctly.
"""

# ╔═╡ 7fcc9a94-7cf4-4fdb-9e22-6e3bacd643a2
begin
	"""Plot vessels with single color, linewidth ∝ radius."""
	function plot_vessels_3d!(ax, seg, n; color=:red, nbins=15)
		n == 0 && return
		raw_w = [clamp(Float32(seg.radius[i] * 2000 / 50), 0.3f0, 8.0f0) for i in 1:n]
		wmin, wmax = extrema(raw_w)
		wmax == wmin && (wmax = wmin + 1f0)

		for b in 1:nbins
			lo = wmin + (b - 1) * (wmax - wmin) / nbins
			hi = wmin + b * (wmax - wmin) / nbins
			bw = (lo + hi) / 2
			pts = Point3f[]
			for i in 1:n
				w = raw_w[i]
				if (b == nbins ? w >= lo : lo <= w < hi)
					push!(pts, Point3f(seg.proximal_x[i], seg.proximal_y[i], seg.proximal_z[i]))
					push!(pts, Point3f(seg.distal_x[i], seg.distal_y[i], seg.distal_z[i]))
				end
			end
			isempty(pts) && continue
			alpha = clamp(0.3 + 0.6 * (bw - wmin) / (wmax - wmin), 0.3, 0.95)
			linesegments!(ax, pts; linewidth=bw, color=(color, alpha))
		end
	end

	"""
	Depth-based visualization matching Wenbo's `plot_tree_3d_depth_with_level_stats`.
	Color = depth from root. Linewidth = radius (percentile-scaled).
	`max_depth` filters how many levels from root to show.
	"""
	function plot_tree_3d_depth!(ax, tree;
			max_depth=typemax(Int), max_segments=typemax(Int),
			lw_min=0.3, lw_max=5.0, alpha=0.9)
		seg = tree.segments
		topo = tree.topology
		n = seg.n
		n == 0 && return 0

		# BFS depth from root
		depth = fill(-1, n)
		queue = Int[]
		for i in 1:n
			if Int(topo.parent_id[i]) <= 0
				depth[i] = 0
				push!(queue, i)
			end
		end
		local h = 1
		while h <= length(queue)
			i = queue[h]; h += 1
			for cid_raw in (topo.child1_id[i], topo.child2_id[i], topo.child3_id[i])
				cid = Int(cid_raw)
				(cid <= 0 || cid > n) && continue
				depth[cid] = depth[i] + 1
				push!(queue, cid)
			end
		end

		# Collect edges filtered by depth
		edges = Tuple{Int,Int,Int}[]
		for i in 1:n
			(depth[i] < 0 || depth[i] >= max_depth) && continue
			for cid_raw in (topo.child1_id[i], topo.child2_id[i], topo.child3_id[i])
				cid = Int(cid_raw)
				(cid <= 0 || cid > n) && continue
				push!(edges, (i, cid, depth[cid]))
			end
		end
		isempty(edges) && return 0

		if length(edges) > max_segments
			idx = sort(Random.shuffle(collect(1:length(edges)))[1:max_segments])
			edges = edges[idx]
		end

		# Color palette per depth
		depth_colors = [
			RGBf(0.8, 0.1, 0.1),    # 0 red
			RGBf(0.0, 0.45, 0.85),   # 1 blue
			RGBf(0.2, 0.7, 0.3),     # 2 green
			RGBf(0.85, 0.55, 0.0),   # 3 orange
			RGBf(0.6, 0.2, 0.8),     # 4 purple
			RGBf(0.0, 0.7, 0.7),     # 5 cyan
			RGBf(0.9, 0.4, 0.6),     # 6 pink
			RGBf(0.5, 0.5, 0.0),     # 7 olive
			RGBf(0.4, 0.4, 0.7),     # 8 slate
			RGBf(0.7, 0.3, 0.3),     # 9 brown
			RGBf(0.5, 0.6, 0.5),     # 10
			RGBf(0.6, 0.4, 0.2),     # 11
		]

		# Linewidth from radius
		radii = [Float64(seg.radius[c]) for (_, c, _) in edges]
		r0, r1 = length(radii) >= 10 ?
			(quantile(radii, 0.05), quantile(radii, 0.95)) : extrema(radii)
		denom = max(r1 - r0, 1e-12)

		edge_lw = [lw_min + clamp((Float64(seg.radius[c]) - r0) / denom, 0.0, 1.0) * (lw_max - lw_min) for (_, c, _) in edges]

		# Draw deepest first (behind), shallowest last (on top)
		max_d = maximum(e[3] for e in edges)
		nbins = 15
		for d in max_d:-1:0
			d_idx = [k for (k, (_, _, ed)) in enumerate(edges) if ed == d]
			isempty(d_idx) && continue

			col = depth_colors[mod1(d + 1, length(depth_colors))]
			a = clamp(alpha - d * 0.04, 0.3, alpha)

			lws_d = [edge_lw[k] for k in d_idx]
			lw_lo_d, lw_hi_d = extrema(lws_d)
			lw_range_d = max(lw_hi_d - lw_lo_d, 1e-12)

			for b in 1:nbins
				lo = lw_lo_d + (b - 1) * lw_range_d / nbins
				hi = lw_lo_d + b * lw_range_d / nbins
				bw = (lo + hi) / 2
				pts = Point3f[]
				for k in d_idx
					w = edge_lw[k]
					if (b == nbins ? w >= lo : lo <= w < hi)
						p, c, _ = edges[k]
						push!(pts, Point3f(seg.proximal_x[p], seg.proximal_y[p], seg.proximal_z[p]))
						push!(pts, Point3f(seg.distal_x[c], seg.distal_y[c], seg.distal_z[c]))
					end
				end
				isempty(pts) && continue
				linesegments!(ax, pts; linewidth=bw, color=(col, a))
			end
		end
		return length(edges)
	end

	"""Plot tree with a single base color, linewidth from radius."""
	function plot_tree_single_color!(ax, tree; color=:red,
			max_depth=typemax(Int), max_segments=typemax(Int),
			lw_min=0.3, lw_max=5.0, alpha=0.9)
		seg = tree.segments
		topo = tree.topology
		n = seg.n
		n == 0 && return 0

		# BFS depth
		depth = fill(-1, n)
		queue = Int[]
		for i in 1:n
			Int(topo.parent_id[i]) <= 0 && (depth[i] = 0; push!(queue, i))
		end
		local h = 1
		while h <= length(queue)
			i = queue[h]; h += 1
			for cid_raw in (topo.child1_id[i], topo.child2_id[i], topo.child3_id[i])
				cid = Int(cid_raw)
				(cid <= 0 || cid > n) && continue
				depth[cid] = depth[i] + 1
				push!(queue, cid)
			end
		end

		edges = Tuple{Int,Int}[]
		for i in 1:n
			(depth[i] < 0 || depth[i] >= max_depth) && continue
			for cid_raw in (topo.child1_id[i], topo.child2_id[i], topo.child3_id[i])
				cid = Int(cid_raw)
				(cid <= 0 || cid > n) && continue
				push!(edges, (i, cid))
			end
		end
		isempty(edges) && return 0

		if length(edges) > max_segments
			idx = sort(Random.shuffle(collect(1:length(edges)))[1:max_segments])
			edges = edges[idx]
		end

		radii = [Float64(seg.radius[c]) for (_, c) in edges]
		r0, r1 = length(radii) >= 10 ?
			(quantile(radii, 0.05), quantile(radii, 0.95)) : extrema(radii)
		denom = max(r1 - r0, 1e-12)

		nbins = 15
		lws = [lw_min + clamp((r - r0) / denom, 0.0, 1.0) * (lw_max - lw_min) for r in radii]
		lw_lo, lw_hi = extrema(lws)
		lw_range = max(lw_hi - lw_lo, 1e-12)

		for b in 1:nbins
			lo = lw_lo + (b - 1) * lw_range / nbins
			hi = lw_lo + b * lw_range / nbins
			bw = (lo + hi) / 2
			pts = Point3f[]
			for (k, (p, c)) in enumerate(edges)
				w = lws[k]
				if (b == nbins ? w >= lo : lo <= w < hi)
					push!(pts, Point3f(seg.proximal_x[p], seg.proximal_y[p], seg.proximal_z[p]))
					push!(pts, Point3f(seg.distal_x[c], seg.distal_y[c], seg.distal_z[c]))
				end
			end
			isempty(pts) && continue
			a = clamp(0.4 + 0.5 * (bw - lw_lo) / lw_range, 0.3, alpha)
			linesegments!(ax, pts; linewidth=bw, color=(color, a))
		end
		return length(edges)
	end

	# Helper: draw shell wireframe
	function draw_shell_wireframe!(ax, dom; nθ=30, nφ=15, color=(:grey70, 0.5))
		θr = range(0, 2π, length=nθ)
		φr = range(0, π, length=nφ)
		a, b, c = dom.semi_axes
		ell_x = [a * sin(p) * cos(t) for t in θr, p in φr]
		ell_y = [b * sin(p) * sin(t) for t in θr, p in φr]
		ell_z = [c * cos(p) for t in θr, p in φr]
		wireframe!(ax, ell_x, ell_y, ell_z; color=color, linewidth=0.3)
	end
end

# ╔═╡ ccce5307-01aa-4f91-a343-1c4d6eaaedf9
begin
	# Grow LAD progressively: 10 → 50 → 150 → 500 terminals
	cfg_lad = configs[1]  # LAD
	stages = [10, 50, 150, 500]
	stage_trees = []

	for nt in stages
		tr = VascularTree("LAD", 5000)
		dx, dy, dz = cfg_lad.root_direction
		dlen = sqrt(dx^2 + dy^2 + dz^2)
		dx /= dlen; dy /= dlen; dz /= dlen
		root_len = cfg_lad.root_radius * 5.0
		distal = (cfg_lad.root_position[1] + dx * root_len,
		          cfg_lad.root_position[2] + dy * root_len,
		          cfg_lad.root_position[3] + dz * root_len)
		add_segment!(tr, cfg_lad.root_position, distal, cfg_lad.root_radius, Int32(-1))
		grow_tree!(tr, domain, nt, params_rca; rng=MersenneTwister(42), kassab=true)
		push!(stage_trees, tr)
	end

	fig_prog = Figure(size=(600, 1200))

	for (col, tr, nt) in zip(1:4, stage_trees, stages)
		ax = Axis3(fig_prog[col, 1],
			title="$nt terminals ($(tr.segments.n) seg)",
			azimuth=1.3, elevation=0.3)
		draw_shell_wireframe!(ax, domain)
		plot_vessels_3d!(ax, tr.segments, tr.segments.n; color=:crimson)
	end

	fig_prog
end

# ╔═╡ 8cd162a5-f0c2-4cbd-8452-34dfc9553904
md"""
## Step 5: Terminal Diameter Monitoring
"""

# ╔═╡ fcb7cba1-34d8-48a6-81cf-0974349cea13
begin
	function terminal_diameter_stats(tree)
		seg = tree.segments
		topo = tree.topology
		n = seg.n
		term_diams = Float64[]
		for i in 1:n
			if topo.is_terminal[i]
				push!(term_diams, seg.radius[i] * 2000.0)  # diameter in μm
			end
		end
		return term_diams
	end

	function pct_below(diams, threshold_um)
		count(d -> d < threshold_um, diams) / length(diams) * 100
	end

	rows = String[]
	for (nt, tr) in zip(stages, stage_trees)
		diams = terminal_diameter_stats(tr)
		n_term = length(diams)
		min_d = round(minimum(diams), digits=1)
		max_d = round(maximum(diams), digits=1)
		mean_d = round(mean(diams), digits=1)
		p90 = round(pct_below(diams, 500.0), digits=1)
		push!(rows, "| $nt | $n_term | $min_d | $mean_d | $max_d | $(p90)% |")
	end
end

# ╔═╡ 5200ec9a-bde4-4e31-b969-46081b35c3f3
Markdown.parse("""
**Terminal diameter monitoring (CCO):**

| Target | Terminals | Min D (μm) | Mean D (μm) | Max D (μm) | % < 500μm |
|:-------|:----------|:-----------|:------------|:-----------|:-----------|
$(join(rows, "\n"))
""")

# ╔═╡ 4718d9a5-2dc9-408d-9afb-4f039765cace
md"""
## Step 6: Full CCO Growth — All 3 Arteries

Simultaneous round-robin CCO growth with inter-tree collision avoidance and territory partitioning. No subdivision — pure CCO matching Wenbo's svVascularize approach.
"""

# ╔═╡ a0109729-5424-438c-8adb-8def433a9cb6
begin
	# CCO-only growth WITH inter-tree collision avoidance
	# handoff_order=0 skips subdivision entirely
	forest = generate_kassab_coronary(
		domain, params_rca;
		rng=MersenneTwister(42),
		verbose=false,
		handoff_order=0,
	)
	cco_trees = Dict(name => tree for (name, tree) in forest.trees)

	cco_total = sum(t.segments.n for (_, t) in cco_trees)
end

# ╔═╡ adc5e9ae-fa33-48f1-98a1-396502a6d4c1
md"""
**CCO result:** $(cco_total) total segments (with inter-tree collision avoidance)

| Artery | Segments | Terminals |
|:-------|:---------|:----------|
| LAD | $(cco_trees["LAD"].segments.n) | $(cco_trees["LAD"].n_terminals) |
| LCX | $(cco_trees["LCX"].segments.n) | $(cco_trees["LCX"].n_terminals) |
| RCA | $(cco_trees["RCA"].segments.n) | $(cco_trees["RCA"].n_terminals) |
"""

# ╔═╡ 28f107f7-1f6f-4329-9f5b-a125f19675f3
md"""
### Per-Artery Views (depth-colored)
"""

# ╔═╡ 941e4826-8025-450b-bd3c-ea7367f423fe
begin
	artery_colors = Dict("LAD" => :crimson, "LCX" => :dodgerblue, "RCA" => :limegreen)
	draw_max_segments = 50000

	fig_per = Figure(size=(1000, 2400))

	for (row, name) in enumerate(["LAD", "LCX", "RCA"])
		tree = cco_trees[name]
		n = tree.segments.n

		ax = Axis3(fig_per[row, 1],
			title="$name — CCO ($n seg)",
			xlabel="X", ylabel="Y", zlabel="Z",
			titlesize=14,
			azimuth=1.3, elevation=0.3)
		plot_tree_3d_depth!(ax, tree;
			max_depth=typemax(Int),
			max_segments=draw_max_segments,
			lw_min=0.3, lw_max=5.0, alpha=0.9)
	end

	fig_per
end

# ╔═╡ d2000001-0000-0000-0000-000000000001
begin
	# Combined 3-artery view — per-artery color to distinguish trees
	fig_3d = Figure(size=(1000, 800))

	ax_combined = Axis3(fig_3d[1, 1],
		title="All Arteries — CCO",
		xlabel="X", ylabel="Y", zlabel="Z",
		titlesize=14,
		azimuth=1.3, elevation=0.3)

	for name in ["RCA", "LCX", "LAD"]
		tree = cco_trees[name]
		plot_tree_single_color!(ax_combined, tree;
			color=artery_colors[name],
			# max_depth=typemax(Int),
			max_depth=20,
			max_segments=draw_max_segments, lw_min=0.3, lw_max=5.0, alpha=0.9)
	end

	fig_3d
end

# ╔═╡ 69e142f1-228e-409d-a096-09b391d3a026
md"""
## Step 7: Terminal Diameter Check
"""

# ╔═╡ 091f9ae7-29f9-4514-b449-91fdc8c31e7e
begin
	thresholds_um = [50.0, 100.0, 200.0, 500.0, 1000.0]

	all_rows = String[]
	for name in ["LAD", "LCX", "RCA"]
		tree = cco_trees[name]
		diams = terminal_diameter_stats(tree)
		n_term = length(diams)
		min_d = round(minimum(diams), digits=1)
		max_d = round(maximum(diams), digits=1)

		pcts = [string(round(pct_below(diams, th), digits=1), "%") for th in thresholds_um]
		push!(all_rows, "| $name | $n_term | $min_d | $max_d | $(join(pcts, " | ")) |")
	end
end

# ╔═╡ 4630961f-8aed-4dbb-b3bc-0c97984d56b1
Markdown.parse("""
**Terminal diameter analysis (CCO):**

| Artery | Terminals | Min D (μm) | Max D (μm) | <50μm | <100μm | <200μm | <500μm | <1000μm |
|:-------|:----------|:-----------|:-----------|:------|:-------|:-------|:-------|:--------|
$(join(all_rows, "\n"))
""")

# ╔═╡ fbb64594-1b4e-4131-bcfb-de324dcc6b65
md"""
## Step 8: Domain Constraint Verification
"""

# ╔═╡ f78616ac-404f-4f80-a495-20fad2b4ea38
begin
	domain_rows = String[]
	for name in ["LAD", "LCX", "RCA"]
		tree = cco_trees[name]
		seg = tree.segments
		n = seg.n
		n_prox_out = 0
		n_dist_out = 0
		for i in 1:n
			pp = (seg.proximal_x[i], seg.proximal_y[i], seg.proximal_z[i])
			dp = (seg.distal_x[i], seg.distal_y[i], seg.distal_z[i])
			!VesselTree.in_domain(domain, pp) && (n_prox_out += 1)
			!VesselTree.in_domain(domain, dp) && (n_dist_out += 1)
		end
		pct_in = round((1 - (n_prox_out + n_dist_out) / (2 * n)) * 100, digits=2)
		push!(domain_rows, "| $name | $n | $n_prox_out | $n_dist_out | $(pct_in)% |")
	end
end

# ╔═╡ a35d2510-dbf9-414f-a36d-e8ea9ac20944
Markdown.parse("""
**Domain constraint check:**

| Artery | Segments | Proximal Outside | Distal Outside | % Inside |
|:-------|:---------|:-----------------|:---------------|:---------|
$(join(domain_rows, "\n"))
""")

# ╔═╡ b8b78b5b-48f4-4507-8a02-d777081cc7de
md"""
## Step 9: Diameter Distribution
"""

# ╔═╡ cb61549f-7fe2-47d5-91fa-ca6b5d7272bd
begin
	fig_hist = Figure(size=(1200, 400))

	for (idx, name) in enumerate(["LAD", "LCX", "RCA"])
		tree = cco_trees[name]
		seg = tree.segments
		n = seg.n
		all_diameters = [seg.radius[i] * 2000 for i in 1:n if seg.radius[i] > 0]

		ax_h = Axis(fig_hist[1, idx],
			title="$name ($(n) seg)",
			xlabel="Diameter (μm)", ylabel="Count",
			xscale=log10,
			xticks=[1, 10, 100, 1000])

		hist!(ax_h, filter(d -> d > 0, all_diameters);
			bins=10.0 .^ range(-0.5, 3.8, length=40),
			color=(artery_colors[name], 0.7))
	end

	fig_hist
end

# ╔═╡ d3000001-0000-0000-0000-000000000001
md"""
## Step 10: Export for Wenbo's Flow Simulation

Export each tree as a text file matching Wenbo's format (11 columns):

```
node_id  parent_id  direction  diameter(μm)  length(μm)  0  0  0  x  y  z
```

This is directly compatible with Wenbo's `flow_simulation_2025_Nov.ipynb` loader.
"""

# ╔═╡ d3000002-0000-0000-0000-000000000001
begin
	export_dir = joinpath(dirname(@__DIR__), "output")
	export_paths = export_forest_wenbo_txt(forest, export_dir)
end

# ╔═╡ 188533a9-2b9b-4e53-89c8-ba8f82ce236e
export_dir

# ╔═╡ d3000003-0000-0000-0000-000000000001
begin
	# Preview first 10 lines of LAD export
	lad_path = joinpath(export_dir, "LAD.txt")
	preview_lines = readlines(lad_path)
	preview = join(preview_lines[1:min(10, length(preview_lines))], "\n")

	md"""
	**Preview (LAD, first 10 lines):**
	```
	$(preview)
	```

	**Column format:** `node_id  parent_id  dir  diameter(μm)  length(μm)  0  0  0  x(mm)  y(mm)  z(mm)`
	"""
end

# ╔═╡ 5f955173-ad5c-48ef-b4aa-dfbbf371713c
md"""
## Summary

Pure CCO growth with Murray's law radii (γ = 7/3) and inter-tree collision avoidance. Domain is an ellipsoid shell modeling the myocardial wall. Text files exported in Wenbo's format for flow simulation.
"""

# ╔═╡ Cell order:
# ╠═d64a065b-258a-47df-86e4-eba5c8966829
# ╠═38d9074f-a4ad-4ae0-9c22-4a57a9ecfaf3
# ╠═6f7c71c6-a2c4-4e69-af63-3cc382df760d
# ╠═86de57c5-bb12-425f-a2d3-218f805f6761
# ╠═e5c54f95-2461-46e7-9870-d411df930300
# ╟─e4965a7a-0800-49fb-93b9-be63caf2111f
# ╟─87657f8c-3c2b-4f06-b89d-b77de7bdae8b
# ╠═bd7d2e5c-f5ea-4dad-a131-0b61d4b99de6
# ╟─7f5bdb8c-d128-4324-9b57-db15adac6de3
# ╟─28d33479-2192-4cc2-8e2f-d1dd2772fb69
# ╠═b26965de-6069-4c3d-a9c9-8be15d74bf10
# ╟─52055f05-592d-473a-88a8-eb4a5034065d
# ╟─6df19e82-fa9f-4e26-bc54-7199d330abb8
# ╠═33573674-ea88-406b-b3c6-9013b976ee37
# ╟─0a987a32-e3d5-44c5-b6a9-d2cf0b29bde1
# ╟─acb6d732-8524-4e60-85d3-a5bc4fa56d86
# ╟─7fcc9a94-7cf4-4fdb-9e22-6e3bacd643a2
# ╟─ccce5307-01aa-4f91-a343-1c4d6eaaedf9
# ╟─8cd162a5-f0c2-4cbd-8452-34dfc9553904
# ╠═fcb7cba1-34d8-48a6-81cf-0974349cea13
# ╟─5200ec9a-bde4-4e31-b969-46081b35c3f3
# ╟─4718d9a5-2dc9-408d-9afb-4f039765cace
# ╠═a0109729-5424-438c-8adb-8def433a9cb6
# ╟─adc5e9ae-fa33-48f1-98a1-396502a6d4c1
# ╟─28f107f7-1f6f-4329-9f5b-a125f19675f3
# ╟─941e4826-8025-450b-bd3c-ea7367f423fe
# ╟─d2000001-0000-0000-0000-000000000001
# ╟─69e142f1-228e-409d-a096-09b391d3a026
# ╠═091f9ae7-29f9-4514-b449-91fdc8c31e7e
# ╟─4630961f-8aed-4dbb-b3bc-0c97984d56b1
# ╟─fbb64594-1b4e-4131-bcfb-de324dcc6b65
# ╠═f78616ac-404f-4f80-a495-20fad2b4ea38
# ╟─a35d2510-dbf9-414f-a36d-e8ea9ac20944
# ╟─b8b78b5b-48f4-4507-8a02-d777081cc7de
# ╟─cb61549f-7fe2-47d5-91fa-ca6b5d7272bd
# ╟─d3000001-0000-0000-0000-000000000001
# ╠═d3000002-0000-0000-0000-000000000001
# ╠═188533a9-2b9b-4e53-89c8-ba8f82ce236e
# ╟─d3000003-0000-0000-0000-000000000001
# ╟─5f955173-ad5c-48ef-b4aa-dfbbf371713c
