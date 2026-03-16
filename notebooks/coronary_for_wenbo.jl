### A Pluto.jl notebook ###
# v0.20.21

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

Pure CCO (Constrained Constructive Optimization) growth matching the svVascularize approach. Murray's law radii (γ = 7/3). Domain is reconstructed from Wenbo's CSV heart volume.
"""

# ╔═╡ 87657f8c-3c2b-4f06-b89d-b77de7bdae8b
md"""
## Step 1: Define the Heart Domain

The domain is built from Wenbo's heart surface CSVs:
- `v4/model_CSVs/heart_points_unique.csv`
- `v4/model_CSVs/heart_normals_unique.csv`

This keeps the existing `VesselTree.jl` growth pipeline intact while replacing
the original analytic ellipsoid shell with the v3 heart volume.

Internally the CSV coordinates are rescaled from Wenbo's centimeter-like units
to VesselTree's millimeter convention so the growth geometry and Murray-law
radii stay in a physically consistent length scale.
"""

# ╔═╡ bd7d2e5c-f5ea-4dad-a131-0b61d4b99de6
domain = default_coronary_volume_domain()

# ╔═╡ e1d44ebb-e7e7-45e5-b649-6d166a0721d6


# ╔═╡ 7f5bdb8c-d128-4324-9b57-db15adac6de3
md"""
**Domain:** `$(nameof(typeof(domain)))`
- Center: $(domain.center)
- Bounding box min: $(domain.min_corner)
- Bounding box max: $(domain.max_corner)
- Monte Carlo volume estimate: $(round(domain.volume, digits=3))
- CSV → internal length scale: $(domain.length_scale)x
- Cached interior sample count: $(size(domain.interior_points, 1))
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

# ╔═╡ 6df19e82-fa9f-4e26-bc54-7199d330abb8
md"""
## Step 3: Growth Parameters

Murray's law exponent γ = 7/3 (Huo-Kassab 2007). Per-artery Kassab morphometric params provide the CCO handoff order and the statistically reconstructed microvascular tree below that handoff.
"""

# ╔═╡ 33573674-ea88-406b-b3c6-9013b976ee37
params_rca = kassab_rca_params();

# ╔═╡ d9f3a362-4a35-41b6-9316-e8efc6753f0a
begin
	handoff_order = 5
	target_terminals = Dict(
		"LAD" => 2250,
		"LCX" => 1050,
		"RCA" => 2500,
	)
	territory_fractions = Dict(
		"LAD" => 0.40,
		"LCX" => 0.25,
		"RCA" => 0.35,
	)
	root_radii_mm = Dict(
		"LAD" => 1.75,
		"LCX" => 2.25,
		"RCA" => 2.25,
	)

	base_configs = VesselTree.coronary_tree_configs(domain)
	tree_configs = [
		VesselTree.TreeConfig(
			cfg.name,
			cfg.root_position,
			get(root_radii_mm, cfg.name, cfg.root_radius),
			cfg.root_direction,
			get(target_terminals, cfg.name, cfg.target_terminals),
			get(territory_fractions, cfg.name, cfg.territory_fraction),
		)
		for cfg in base_configs
	]
end

# ╔═╡ b26965de-6069-4c3d-a9c9-8be15d74bf10
begin
	configs = tree_configs

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

# ╔═╡ 0a987a32-e3d5-44c5-b6a9-d2cf0b29bde1
md"""
γ = $(params_rca.gamma), $(params_rca.n_orders) Strahler orders

Paper-consistent default for flow comparison:
- `handoff_order = $(handoff_order)` so each CCO terminal is statistically subdivided into the downstream microvascular tree.
- `handoff_order = 0` leaves only the CCO skeleton, which gives far too few parallel terminal pathways to compare against the 2007/2008 porcine reconstructions.
- `target_terminals = $(target_terminals)` sets the CCO skeleton resolution before statistical subdivision.
- `territory_fractions = $(territory_fractions)` sets the intended myocardial territory split.
- `root_radii_mm = $(root_radii_mm)` applies the paper-derived order-11 root-radius scale.

Edit only this one parameter block when calibrating toward the paper:
- `handoff_order`
- `target_terminals["LAD"]`
- `target_terminals["LCX"]`
- `target_terminals["RCA"]`
- `territory_fractions["LAD"]`
- `territory_fractions["LCX"]`
- `territory_fractions["RCA"]`
- `root_radii_mm["LAD"]`
- `root_radii_mm["LCX"]`
- `root_radii_mm["RCA"]`
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
			lw_min=0.3, lw_max=5.0, alpha=0.9,
			min_diameter_um=0.0)
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
				(2.0 * Float64(seg.radius[cid]) * 1000.0) < min_diameter_um && continue
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
			lw_min=0.3, lw_max=5.0, alpha=0.9,
			min_diameter_um=0.0)
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
				(2.0 * Float64(seg.radius[cid]) * 1000.0) < min_diameter_um && continue
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

	# Helper: draw analytic shell or point-cloud surface context
	function draw_shell_wireframe!(ax, dom; nθ=30, nφ=15, color=(:grey70, 0.5))
		if dom isa CSVVolumeDomain
			n = size(dom.surface_points, 1)
			step = max(1, cld(n, 2500))
			pts = dom.surface_points[1:step:end, :]
			scatter!(ax, pts[:, 1], pts[:, 2], pts[:, 3];
				color=color, markersize=2)
			return
		end

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
		if !in_domain(domain, distal)
			distal = project_to_domain(domain, distal)
		end
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
## Step 6: Full Coronary Growth — All 3 Arteries

Simultaneous round-robin CCO growth with inter-tree collision avoidance and territory partitioning, followed by Kassab connectivity-matrix subdivision below the handoff order.

This keeps the existing v4 growth model intact:
- CCO builds the resolving skeleton in the v3 heart volume
- Murray's law sets radii throughout
- Kassab morphometry fills in the distal microvascular tree without hand-tuning diameters for flow targets
"""

# ╔═╡ a0109729-5424-438c-8adb-8def433a9cb6
begin
	forest = generate_kassab_coronary(
		domain, params_rca;
		rng=MersenneTwister(42),
		verbose=false,
		handoff_order=handoff_order,
		tree_configs=tree_configs,
	)
	cco_trees = Dict(name => tree for (name, tree) in forest.trees)

	cco_total = sum(t.segments.n for (_, t) in cco_trees)
end

# ╔═╡ adc5e9ae-fa33-48f1-98a1-396502a6d4c1
md"""
**Tree generation result:** $(cco_total) total segments (`handoff_order = $(handoff_order)`)

**Configured CCO target terminals:** `$(target_terminals)`
**Configured territory fractions:** `$(territory_fractions)`

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

# ╔═╡ d2100001-0000-0000-0000-000000000001
begin
	macro_max_depth = 14
	macro_min_diameter_um = 180.0
end

# ╔═╡ 941e4826-8025-450b-bd3c-ea7367f423fe
begin
	artery_colors = Dict("LAD" => :crimson, "LCX" => :dodgerblue, "RCA" => :limegreen)
	draw_max_segments = 50000

	fig_per = Figure(size=(1000, 2400))

	for (row, name) in enumerate(["LAD", "LCX", "RCA"])
		tree = cco_trees[name]
		n = tree.segments.n

		ax = Axis3(fig_per[row, 1],
			title="$name — full tree ($n seg)",
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

# ╔═╡ d2100002-0000-0000-0000-000000000001
md"""
### Macro Trunk Views

For million-segment trees, the full-depth plot becomes a microvascular cloud.
This view shows only the larger, proximal branches:
- `max_depth = $(macro_max_depth)`
- `min_diameter = $(macro_min_diameter_um) μm`
"""

# ╔═╡ d2100003-0000-0000-0000-000000000001
begin
	fig_macro = Figure(size=(1000, 2400))

	for (row, name) in enumerate(["LAD", "LCX", "RCA"])
		tree = cco_trees[name]
		n = tree.segments.n

		ax = Axis3(fig_macro[row, 1],
			title="$name — macro trunk view ($n seg total)",
			xlabel="X", ylabel="Y", zlabel="Z",
			titlesize=14,
			azimuth=1.3, elevation=0.3)
		draw_shell_wireframe!(ax, domain)
		plot_tree_3d_depth!(ax, tree;
			max_depth=macro_max_depth,
			min_diameter_um=macro_min_diameter_um,
			max_segments=draw_max_segments,
			lw_min=0.8, lw_max=8.0, alpha=0.95)
	end

	fig_macro
end

# ╔═╡ d2000001-0000-0000-0000-000000000001
begin
	# Combined 3-artery view — per-artery color to distinguish trees
	fig_3d = Figure(size=(1000, 800))

	ax_combined = Axis3(fig_3d[1, 1],
		title="All Arteries — macro trunk view",
		xlabel="X", ylabel="Y", zlabel="Z",
		titlesize=14,
		azimuth=1.3, elevation=0.3)

	for name in ["RCA", "LCX", "LAD"]
		tree = cco_trees[name]
		plot_tree_single_color!(ax_combined, tree;
			color=artery_colors[name],
			max_depth=macro_max_depth,
			min_diameter_um=macro_min_diameter_um,
			max_segments=draw_max_segments, lw_min=0.8, lw_max=8.0, alpha=0.95)
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
node_id  parent_id  direction  diameter(μm)  length(cm)  0  0  0  x(cm)  y(cm)  z(cm)
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

	**Column format:** `node_id  parent_id  dir  diameter(μm)  length(cm)  0  0  0  x(cm)  y(cm)  z(cm)`
	"""
end

# ╔═╡ 5f955173-ad5c-48ef-b4aa-dfbbf371713c
md"""
## Summary

CCO skeleton growth with Murray's law radii (γ = 7/3), Kassab statistical subdivision, and inter-tree collision avoidance inside Wenbo's reconstructed heart volume. Text files are exported in Wenbo's format for flow simulation.
"""

# ╔═╡ f1000001-0000-0000-0000-000000000001
md"""
## Step 11: Flow Simulation

**1:1 port of Wenbo's Python `flow_simulation_2025_Nov.ipynb`.**

Reads the exported `.txt` files, builds a binary tree, computes Poiseuille resistance with Pries (1994) in-vivo viscosity, prunes oversized leaf nodes ("leakages"), and solves for flow/pressure.

This port keeps Wenbo's file convention of `diameter(μm)` and `length(cm)`.
For comparison with the 2008 Molloi/Wong normal-tree model, the outlet
pressure is set to pre-capillary pressure `Pc = Pv + 10 mmHg = 15 mmHg`
when venous pressure is `Pv = 5 mmHg`.

Vasodilation option:
- the flow solver exposes a generic `vasodilation(...)` function
- default hyperemic parameters matching the 2008 paper are:
  `factor = 1.6`, `min_diameter = 0 μm`, `max_diameter = 400 μm`
- this changes only the hemodynamic evaluation, not the generated tree geometry
  or Murray-law radii stored in the tree/export

**Changes from Python original:**
- Preserves Wenbo's `length(cm) → meters` conversion via `length × 0.01`
- Uses pre-capillary outlet pressure `Pc = 15 mmHg` for arterial-tree flow
- Supports optional vasodilation via a configurable diameter-band scaling function
- Added pruning of oversized leaf nodes (leaked large vessels that were never subdivided)
- Handles single-child nodes (series connection) after pruning

| Python constant | Value | Julia name |
|:----------------|:------|:-----------|
| `Pin` | 100 mmHg × 133.32 = 13332.24 Pa | `FLOW_Pin` |
| `Pout` | 15 mmHg × 133.32 = 1999.84 Pa | `FLOW_Pout` |
| `TOLERANCE` | 0.001 | `FLOW_TOLERANCE` |
| `micron2m` | 1e-6 | `FLOW_micron2m` |
| `m3pers2mLpermin` | 6e7 | `FLOW_m3s_to_mLmin` |
| `vasodilation` | `factor=1.6, min=0 μm, max=400 μm` | `FLOW_ENABLE_VASODILATION` |
"""

# ╔═╡ f1000002-0000-0000-0000-000000000001
begin
	"""Vessel segment node for flow simulation (matches Python `Node` class)."""
	mutable struct FlowNode
		label::Int
		parent_label::Int
		diameter::Float64        # μm
		len::Float64             # cm (Wenbo export convention)
		direction::String
		Rs::Float64              # Segment resistance (N·s/m⁵)
		Rc::Float64              # Crown resistance (subtree equivalent)
		Rs_decoupled::Float64    # Resistance used in iteration
		pressure::Float64        # N/m²
		flow::Float64            # m³/s
		flow_old::Float64        # Previous iteration flow
		parent::Union{FlowNode, Nothing}
		left_branch::Union{FlowNode, Nothing}
		right_branch::Union{FlowNode, Nothing}

		function FlowNode(label, parent_label, diameter, len, direction)
			new(label, parent_label, diameter, len, direction,
				0.0, 0.0, 0.0,    # Rs, Rc, Rs_decoupled
				0.0, 0.0, 0.0,    # pressure, flow, flow_old
				nothing, nothing, nothing)
		end
	end

	# Constants (matching Wenbo's Python exactly)
	FLOW_mmHg2Pa      = 133.32239
	FLOW_micron2m     = 1e-6
	FLOW_m3s_to_mLmin = 1e6 * 60.0
	FLOW_TOLERANCE    = 0.001
	FLOW_Pin          = 100 * FLOW_mmHg2Pa   # 100 mmHg inlet
	FLOW_Pout         = 15 * FLOW_mmHg2Pa    # pre-capillary outlet (Pv=5 mmHg, Pc=15 mmHg)
	FLOW_ENABLE_VASODILATION = false
	FLOW_VASODILATION_MIN_DIAMETER_UM = 0.0
	FLOW_VASODILATION_MAX_DIAMETER_UM = 400.0
	FLOW_VASODILATION_FACTOR = 1.6
end

# ╔═╡ f1000003-0000-0000-0000-000000000001
begin
	# ── Pries (1994) in-vivo viscosity ──────────────────────────────────
	# 1:1 port of Python calculate_viscosity_pries(diameter_um)
	function flow_viscosity_pries(diameter_um::Float64)
		D = max(diameter_um, 2.0)
		Hd = 0.45
		mu_plasma = 0.012   # Poise

		term1 = 6.0 * exp(-0.085 * D)
		term2 = 3.2 - 2.44 * exp(-0.06 * D^0.645)
		mu_045_star = term1 + term2

		term_c = (0.8 + exp(-0.075 * D)) *
				 (-1.0 + 1.0 / (1.0 + 1e-11 * D^12)) +
				 1.0 / (1.0 + 1e-11 * D^12)

		viscosity_rel = (1.0 + (mu_045_star - 1.0) *
						 ((1.0 - Hd)^term_c - 1.0) /
						 ((1.0 - 0.45)^term_c - 1.0) *
						 (D / (D - 1.1))^2) *
						(D / (D - 1.1))^2

		return viscosity_rel * mu_plasma * 0.1   # Poise → Pa·s
	end

	function vasodilation(
		diameter_um::Float64;
		factor::Float64=FLOW_VASODILATION_FACTOR,
		min_diameter_um::Float64=FLOW_VASODILATION_MIN_DIAMETER_UM,
		max_diameter_um::Float64=FLOW_VASODILATION_MAX_DIAMETER_UM,
	)
		if min_diameter_um <= diameter_um <= max_diameter_um
			return diameter_um * factor
		end
		return diameter_um
	end

	function flow_effective_diameter_um(diameter_um::Float64)
		if FLOW_ENABLE_VASODILATION
			return vasodilation(
				diameter_um;
				factor=FLOW_VASODILATION_FACTOR,
				min_diameter_um=FLOW_VASODILATION_MIN_DIAMETER_UM,
				max_diameter_um=FLOW_VASODILATION_MAX_DIAMETER_UM,
			)
		end
		return diameter_um
	end

	# ── Two-pass tree loader ───────────────────────────────────────────
	# 1:1 port of Python load_unordered_tree(file_path)
	function flow_load_tree(file_path::String)
		println("Scanning file: $file_path...")
		isfile(file_path) || (println("Error: File not found."); return nothing)

		node_map = Dict{Int, FlowNode}()

		# PASS 1: Create all nodes
		for (line_num, line) in enumerate(eachline(file_path))
			words = split(strip(line))
			length(words) < 5 && continue
			startswith(words[1], "#") && continue

			label = parse(Int, words[1])
			parent_label = parse(Int, words[2])
			direction = lowercase(String(words[3]))
			diameter = parse(Float64, words[4])
			seg_len = parse(Float64, words[5])

			if haskey(node_map, label)
				println("Warning: Duplicate node label $label at line $line_num")
				continue
			end
			node_map[label] = FlowNode(label, parent_label, diameter, seg_len, direction)
		end

		println("File loaded. Found $(length(node_map)) nodes. Linking tree structure...")

		# PASS 2: Link parent ↔ child
		tree_root = nothing
		linked_count = 0

		for (label, node) in node_map
			if node.parent_label == -1
				tree_root = node
				println("Root found: Node $(node.label)")
				continue
			end
			parent = get(node_map, node.parent_label, nothing)
			if parent !== nothing
				node.parent = parent
				if node.direction == "l"
					parent.left_branch = node
				elseif node.direction == "r"
					parent.right_branch = node
				end
				linked_count += 1
			else
				println("Warning: Node $label has missing parent $(node.parent_label)")
			end
		end

		if tree_root === nothing
			println("Critical Error: No root node (parent -1) found!")
			return nothing
		end
		println("Tree linked successfully. $linked_count connections made.")
		return tree_root
	end

	# ── Prune oversized leaf nodes ─────────────────────────────────────
	# Removes "leaked" large vessels that were never subdivided
	function flow_prune_leaves!(root::FlowNode; max_leaf_diameter_um::Float64=50.0)
		pruned = 0
		stack = FlowNode[root]
		while !isempty(stack)
			node = pop!(stack)
			# Check left child
			lc = node.left_branch
			if lc !== nothing
				is_leaf = lc.left_branch === nothing && lc.right_branch === nothing
				if is_leaf && lc.diameter > max_leaf_diameter_um
					node.left_branch = nothing
					pruned += 1
				else
					push!(stack, lc)
				end
			end
			# Check right child
			rc = node.right_branch
			if rc !== nothing
				is_leaf = rc.left_branch === nothing && rc.right_branch === nothing
				if is_leaf && rc.diameter > max_leaf_diameter_um
					node.right_branch = nothing
					pruned += 1
				else
					push!(stack, rc)
				end
			end
		end
		println("Pruned $pruned oversized leaf nodes (diameter > $max_leaf_diameter_um μm)")
		return pruned
	end

	# ── Resistance calculation ─────────────────────────────────────────
	# 1:1 port of Python update_tree_physics(current_node)
	# Rs = 128·μ·L / (π·D⁴)  (Poiseuille's law)
	# Rc = Rs + parallel(Rc_left, Rc_right) for bifurcations
	#    = Rs + Rc_child for single-child nodes (after pruning)
	#    = Rs for leaf nodes
	function flow_update_physics!(node::Union{FlowNode, Nothing})
		node === nothing && return 1e-20

		eff_diameter = flow_effective_diameter_um(node.diameter)
		mu = flow_viscosity_pries(eff_diameter)
		L_m = node.len * 0.01                   # cm → m (Wenbo file convention)
		D_m = eff_diameter * FLOW_micron2m      # μm → m
		L_m <= 0.0 && (L_m = 1e-9)

		Rs_SI = (128.0 * mu * L_m) / (π * D_m^4)
		Rs_SI < 1e-20 && (Rs_SI = 1e-20)
		node.Rs = Rs_SI
		node.Rs_decoupled = Rs_SI

		Rc_left  = flow_update_physics!(node.left_branch)
		Rc_right = flow_update_physics!(node.right_branch)

		has_left  = node.left_branch !== nothing
		has_right = node.right_branch !== nothing

		if !has_left && !has_right
			node.Rc = node.Rs
		elseif has_left && has_right
			R_par = (Rc_left * Rc_right) / (Rc_left + Rc_right)
			node.Rc = node.Rs + R_par
		elseif has_left
			node.Rc = node.Rs + Rc_left
		else
			node.Rc = node.Rs + Rc_right
		end
		return node.Rc
	end

	# ── Flow/pressure propagation ─────────────────────────────────────
	# 1:1 port of Python calculate_flow_pressure(upstream_node, error_sq)
	function flow_calc_pressure!(node::FlowNode, error_sq::Float64)
		has_left  = node.left_branch !== nothing
		has_right = node.right_branch !== nothing

		if !has_left && !has_right
			error_sq += (node.flow - node.flow_old)^2
			node.flow_old = node.flow
			return error_sq
		end

		pressure_drop = node.flow * node.Rs_decoupled
		P_bif = node.pressure - pressure_drop
		ΔP = P_bif - FLOW_Pout

		if has_left && has_right
			node.left_branch.pressure  = P_bif
			node.right_branch.pressure = P_bif
			node.left_branch.flow  = ΔP / max(node.left_branch.Rc, 1e-20)
			node.right_branch.flow = ΔP / max(node.right_branch.Rc, 1e-20)
		elseif has_left
			node.left_branch.pressure = P_bif
			node.left_branch.flow = ΔP / max(node.left_branch.Rc, 1e-20)
		else
			node.right_branch.pressure = P_bif
			node.right_branch.flow = ΔP / max(node.right_branch.Rc, 1e-20)
		end

		error_sq += (node.flow - node.flow_old)^2
		node.flow_old = node.flow

		has_left  && (error_sq = flow_calc_pressure!(node.left_branch, error_sq))
		has_right && (error_sq = flow_calc_pressure!(node.right_branch, error_sq))
		return error_sq
	end

	# ── Leaf node analysis ─────────────────────────────────────────────
	# 1:1 port of Python analyze_leaf_nodes(root_node)
	function flow_analyze_leaves(root::FlowNode)
		leaves = FlowNode[]
		stack = FlowNode[root]
		while !isempty(stack)
			n = pop!(stack)
			if n.left_branch === nothing && n.right_branch === nothing
				push!(leaves, n)
			else
				n.left_branch  !== nothing && push!(stack, n.left_branch)
				n.right_branch !== nothing && push!(stack, n.right_branch)
			end
		end

		diams = [n.diameter for n in leaves]
		d_min, d_max = extrema(diams)
		d_avg = sum(diams) / length(diams)
		n_large = count(d -> d > 9.0, diams)
		pct = round(n_large / length(diams) * 100; digits=2)

		println("Total Leaf Nodes: $(length(leaves))")
		println("Leaf Diameter Stats (microns):")
		println("  Min: $(round(d_min; digits=2))")
		println("  Max: $(round(d_max; digits=2))  <-- This number should not be really big")
		println("  Avg: $(round(d_avg; digits=2))")
		println("Count of 'Large Leaves' (>9.00 μm): $n_large ($pct%)")
		return leaves
	end

	# ── Main driver ────────────────────────────────────────────────────
	# 1:1 port of Python run_dcad_simulation(input_file_path)
	function flow_run_simulation(file_path::String; max_leaf_diameter_um::Float64=50.0)
		# 1. Load tree
		root = flow_load_tree(file_path)
		root === nothing && (println("Simulation aborted."); return nothing, nothing)

		# 2. Prune oversized leaves
		flow_prune_leaves!(root; max_leaf_diameter_um)

		# 3. Calculate physics (Rs & Rc)
		println("Calculating physics (Rs & Rc)...")
		total_resistance = flow_update_physics!(root)
		println("Total Tree Resistance: $(round(total_resistance; sigdigits=5)) N*s/m^5")

		# 4. Initialize flow
		root.pressure = FLOW_Pin
		root.flow = (FLOW_Pin - FLOW_Pout) / root.Rc
		root.flow_old = root.flow
		flow_mLmin = root.flow * FLOW_m3s_to_mLmin
		println("Initial Inlet Flow: $(round(flow_mLmin; digits=4)) mL/min")

		# 5. Iterative solver (matches Python loop exactly)
		println("Starting simulation loop...")
		max_iterations = 50
		converged = false

		for iteration in 1:max_iterations
			error_sq = flow_calc_pressure!(root, 0.0)
			std_error = sqrt(error_sq)

			iteration % 10 == 0 && println("Iteration $iteration: Error = $std_error")

			if std_error < FLOW_TOLERANCE
				converged = true
				println("Simulation finished in $iteration iterations.")
				break
			end
		end
		converged || println("Simulation finished at max iterations ($max_iterations).")

		final_flow = root.flow * FLOW_m3s_to_mLmin
		println("Final Inlet Flow: $(round(final_flow; digits=4)) mL/min")

		return root, final_flow
	end

	md"Flow simulation functions loaded."
end

# ╔═╡ f1000004-0000-0000-0000-000000000001
md"""
### Run Flow Simulation on Exported Trees

Reads the `.txt` files exported in Step 10 and runs Wenbo's flow simulation on each artery.
"""

# ╔═╡ f1000005-0000-0000-0000-000000000001
begin
	flow_results = Dict{String, Tuple{FlowNode, Float64}}()

	for name in ["LAD", "LCX", "RCA"]
		txt_path = joinpath(export_dir, "$name.txt")
		if isfile(txt_path)
			println("\n" * "="^60)
			println("  $name")
			println("="^60)
			root_node, final_flow = flow_run_simulation(txt_path; max_leaf_diameter_um=50.0)
			if root_node !== nothing
				flow_results[name] = (root_node, final_flow)
				println()
				flow_analyze_leaves(root_node)
			end
		else
			println("Skipping $name — file not found: $txt_path")
		end
	end

	println("\n" * "="^60)
	println("  SUMMARY")
	println("="^60)
	for name in ["LAD", "LCX", "RCA"]
		if haskey(flow_results, name)
			_, flow_val = flow_results[name]
			println("  $name: $(round(flow_val; digits=4)) mL/min")
		end
	end
end

# ╔═╡ f1000006-0000-0000-0000-000000000001
begin
	flow_summary_rows = String[]
	for name in ["LAD", "LCX", "RCA"]
		if haskey(flow_results, name)
			root_n, flow_val = flow_results[name]
			leaves = flow_analyze_leaves(root_n)
			diams_f = [n.diameter for n in leaves]
			d_min_f, d_max_f = extrema(diams_f)
			d_avg_f = round(sum(diams_f) / length(diams_f); digits=1)
			push!(flow_summary_rows, "| $name | $(round(flow_val; digits=4)) | $(length(leaves)) | $(round(d_min_f; digits=1)) | $(round(d_max_f; digits=1)) | $d_avg_f |")
		end
	end
end

# ╔═╡ 6277c28f-6ebe-4209-81d6-19e75cadc3e2
Markdown.parse("""
**Flow Simulation Results (Pries viscosity, 100 mmHg inlet, 15 mmHg pre-capillary outlet, vasodilation enabled = $(FLOW_ENABLE_VASODILATION)):**

| Artery | Flow (mL/min) | Leaf Nodes | Min D (μm) | Max D (μm) | Avg D (μm) |
|:-------|:--------------|:-----------|:-----------|:-----------|:-----------|
$(join(flow_summary_rows, "\n"))

*Oversized leaf nodes (>50 μm diameter) pruned before simulation.*

$(FLOW_ENABLE_VASODILATION ?
"*Vasodilation enabled: segments with diameters in [$(Int(FLOW_VASODILATION_MIN_DIAMETER_UM)), $(Int(FLOW_VASODILATION_MAX_DIAMETER_UM))] μm use a $(FLOW_VASODILATION_FACTOR)× dilated diameter during resistance calculation.*" :
"*Vasodilation disabled: flow uses the exported diameters directly.*")
""")

# ╔═╡ Cell order:
# ╠═d64a065b-258a-47df-86e4-eba5c8966829
# ╠═38d9074f-a4ad-4ae0-9c22-4a57a9ecfaf3
# ╠═6f7c71c6-a2c4-4e69-af63-3cc382df760d
# ╠═86de57c5-bb12-425f-a2d3-218f805f6761
# ╠═e5c54f95-2461-46e7-9870-d411df930300
# ╟─e4965a7a-0800-49fb-93b9-be63caf2111f
# ╟─87657f8c-3c2b-4f06-b89d-b77de7bdae8b
# ╠═bd7d2e5c-f5ea-4dad-a131-0b61d4b99de6
# ╠═e1d44ebb-e7e7-45e5-b649-6d166a0721d6
# ╟─7f5bdb8c-d128-4324-9b57-db15adac6de3
# ╟─28d33479-2192-4cc2-8e2f-d1dd2772fb69
# ╠═b26965de-6069-4c3d-a9c9-8be15d74bf10
# ╟─52055f05-592d-473a-88a8-eb4a5034065d
# ╟─6df19e82-fa9f-4e26-bc54-7199d330abb8
# ╠═33573674-ea88-406b-b3c6-9013b976ee37
# ╠═d9f3a362-4a35-41b6-9316-e8efc6753f0a
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
# ╠═d2100001-0000-0000-0000-000000000001
# ╟─941e4826-8025-450b-bd3c-ea7367f423fe
# ╟─d2100002-0000-0000-0000-000000000001
# ╟─d2100003-0000-0000-0000-000000000001
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
# ╟─f1000001-0000-0000-0000-000000000001
# ╠═f1000002-0000-0000-0000-000000000001
# ╠═f1000003-0000-0000-0000-000000000001
# ╟─f1000004-0000-0000-0000-000000000001
# ╠═f1000005-0000-0000-0000-000000000001
# ╠═f1000006-0000-0000-0000-000000000001
# ╟─6277c28f-6ebe-4209-81d6-19e75cadc3e2
