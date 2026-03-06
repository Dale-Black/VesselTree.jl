### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ a1000001-b000-4c00-8d00-e00000000002
begin
	import Pkg
	Pkg.activate(dirname(@__DIR__))
end

# ╔═╡ f4c0b7fb-637d-4bd7-adff-4a06fe91c081
using VesselTree

# ╔═╡ 4378fe0a-1dd5-4b66-b6ae-b0668e12b944
using CairoMakie

# ╔═╡ fb9ffd67-0ab9-4974-ba52-d64742a58a08
using Random, Statistics

# ╔═╡ 80b575e3-d4f9-4250-904f-19e88f059029
md"""
# Coronary Vasculature Generator

Generate a complete coronary arterial tree from ~3mm stems down to **8 micrometer capillaries** using the hybrid CCO + Kassab subdivision pipeline.

- **CCO** (Constrained Constructive Optimization) -- space-filling skeleton growth
- **Kassab 1993 morphometry** -- connectivity matrix, per-order diameter/length distributions
- **Barabasi 2026 surface optimization** -- junction geometry (sprouting vs branching)
- **Murray's law** -- gamma = 7/3 (Huo-Kassab 2007)
- **AcceleratedKernels.jl** -- parallel distance/intersection kernels
"""

# ╔═╡ dcbc6e0f-c3c5-4ef1-86c4-64ca5a5bc850
md"""
## Step 1: CCO Skeleton Growth

First we grow a space-filling skeleton using CCO. Each iteration adds one terminal with intersection checking and Murray's law radius updates. We show progressive growth at 10, 30, and 100 terminals.
"""

# ╔═╡ 31d6ee6d-896a-4dab-a0fe-98e551461314
begin
	params_rca = kassab_rca_params()
	domain = SphereDomain((0.0, 0.0, 0.0), 50.0)
	rng = MersenneTwister(42)
end;

# ╔═╡ 9201e42e-a20a-494b-a34b-316f0c536098
md"""
### Progressive CCO Growth
"""

# ╔═╡ 3c536f7a-9713-40f5-bac9-645cb3998d1a
"""
Plot vessel segments with linewidth proportional to radius. Batches by discretized
linewidth since CairoMakie needs uniform width per linesegments! call.
"""
function plot_vessels!(ax, seg, n; color=:red, nbins=15, alpha=0.85)
	n == 0 && return

	# Pick the two axes with the most spatial extent for the projection
	xr = extrema(seg.proximal_x[i] for i in 1:n)
	yr = extrema(seg.proximal_y[i] for i in 1:n)
	zr = extrema(seg.proximal_z[i] for i in 1:n)
	spans = [(xr[2]-xr[1], :x), (yr[2]-yr[1], :y), (zr[2]-zr[1], :z)]
	sort!(spans, by=first, rev=true)
	ax1, ax2 = spans[1][2], spans[2][2]
	getcoord(field, i) = field == :x ? (seg.proximal_x[i], seg.distal_x[i]) :
	                     field == :y ? (seg.proximal_y[i], seg.distal_y[i]) :
	                                   (seg.proximal_z[i], seg.distal_z[i])

	raw_widths = [clamp(Float32(sqrt(seg.radius[i]) * 8), 0.3f0, 8.0f0) for i in 1:n]
	wmin, wmax = extrema(raw_widths)
	wmax == wmin && (wmax = wmin + 1f0)

	for b in 1:nbins
		lo = wmin + (b - 1) * (wmax - wmin) / nbins
		hi = wmin + b * (wmax - wmin) / nbins
		bin_w = (lo + hi) / 2

		pts = Point2f[]
		for i in 1:n
			w = raw_widths[i]
			if (b == nbins ? w >= lo : lo <= w < hi)
				p1, d1 = getcoord(ax1, i)
				p2, d2 = getcoord(ax2, i)
				push!(pts, Point2f(p1, p2))
				push!(pts, Point2f(d1, d2))
			end
		end
		isempty(pts) && continue
		linesegments!(ax, pts; linewidth=bin_w, color=(color, alpha))
	end
end

# ╔═╡ a2000001-b000-4c00-8d00-e00000000003
begin
	# Root: start at top of sphere, grow downward
	root_pos = (0.0, 0.0, 50.0)
	root_end = (0.0, 0.0, 42.5)  # 7.5mm root segment along -z

	# Grow CCO skeleton at 3 stages
	tree_10 = VascularTree("RCA", 5000)
	add_segment!(tree_10, root_pos, root_end, 1.5, Int32(-1))
	grow_tree!(tree_10, domain, 10, params_rca; rng=MersenneTwister(42))

	tree_30 = VascularTree("RCA", 5000)
	add_segment!(tree_30, root_pos, root_end, 1.5, Int32(-1))
	grow_tree!(tree_30, domain, 30, params_rca; rng=MersenneTwister(42))

	tree_100 = VascularTree("RCA", 5000)
	add_segment!(tree_100, root_pos, root_end, 1.5, Int32(-1))
	grow_tree!(tree_100, domain, 100, params_rca; rng=MersenneTwister(42))

	md"CCO skeleton: **$(tree_10.segments.n)** -> **$(tree_30.segments.n)** -> **$(tree_100.segments.n)** segments"
end

# ╔═╡ a2000002-b000-4c00-8d00-e00000000004
begin
	fig_cco = Figure(size=(1200, 400), backgroundcolor=:grey10)

	for (idx, (tree_s, label)) in enumerate([
		(tree_10, "10 terminals"),
		(tree_30, "30 terminals"),
		(tree_100, "100 terminals"),
	])
		ax = Axis(fig_cco[1, idx],
			title=label,
			xlabel="mm", ylabel="mm",
			backgroundcolor=:grey10,
			titlecolor=:white, xlabelcolor=:white, ylabelcolor=:white,
			xticklabelcolor=:grey70, yticklabelcolor=:grey70,
			aspect=DataAspect())
		plot_vessels!(ax, tree_s.segments, tree_s.segments.n; color=:limegreen)
	end

	fig_cco
end

# ╔═╡ b3000001-c000-4d00-9e00-f00000000001
md"""
## 3D Coronary Forest

Three arteries (LAD, LCX, RCA) grown with different root directions. Each artery fills its own plane around the heart surface -- together they form a true 3D vascular network. Subdivision extends each to capillary scale (~8 um).
"""

# ╔═╡ b3000002-c000-4d00-9e00-f00000000002
begin
	domain_3d = EllipsoidDomain((0.0, 0.0, 0.0), (50.0, 40.0, 35.0))
	configs_3d = [
		TreeConfig("LAD", (-5.0, 5.0, 30.0), 1.5, (1.0, -1.0, -1.0), 50, 0.40),
		TreeConfig("LCX", (-5.0, -5.0, 30.0), 1.2, (1.0, 1.0, -1.0), 30, 0.25),
		TreeConfig("RCA", (5.0, 0.0, 30.0), 1.3, (-1.0, 0.0, -1.0), 40, 0.35),
	]
	t_3d = @elapsed forest_3d = generate_kassab_coronary(
		domain_3d, params_rca;
		tree_configs=configs_3d,
		rng=MersenneTwister(42),
		verbose=false,
	)
	total_3d = sum(t.segments.n for (_, t) in forest_3d.trees)

	md"""
	**3-artery forest:** $(total_3d) segments in $(round(t_3d, digits=1))s
	"""
end

# ╔═╡ b3000003-c000-4d00-9e00-f00000000003
begin
	fig_3d = Figure(size=(1000, 800), backgroundcolor=:grey10)
	ax3d = Axis3(fig_3d[1, 1],
		title="Coronary Arterial Forest (3D)",
		xlabel="x (mm)", ylabel="y (mm)", zlabel="z (mm)",
		backgroundcolor=:grey10,
		titlecolor=:white,
		xlabelcolor=:white, ylabelcolor=:white, zlabelcolor=:white,
		xticklabelcolor=:grey70, yticklabelcolor=:grey70, zticklabelcolor=:grey70,
		azimuth=1.2, elevation=0.35,
	)

	artery_colors_3d = Dict("LAD" => :crimson, "LCX" => :dodgerblue, "RCA" => :limegreen)

	for (name, tree) in forest_3d.trees
		seg = tree.segments
		n = tree.segments.n
		col = artery_colors_3d[name]

		raw_w = [clamp(Float32(sqrt(seg.radius[i]) * 8), 0.3f0, 6.0f0) for i in 1:n]
		wmin, wmax = extrema(raw_w)
		wmax == wmin && (wmax = wmin + 1f0)

		for b in 1:12
			lo = wmin + (b - 1) * (wmax - wmin) / 12
			hi = wmin + b * (wmax - wmin) / 12
			bw = (lo + hi) / 2

			pts = Point3f[]
			for i in 1:n
				w = raw_w[i]
				if (b == 12 ? w >= lo : lo <= w < hi)
					push!(pts, Point3f(seg.proximal_x[i], seg.proximal_y[i], seg.proximal_z[i]))
					push!(pts, Point3f(seg.distal_x[i], seg.distal_y[i], seg.distal_z[i]))
				end
			end
			isempty(pts) && continue
			linesegments!(ax3d, pts; linewidth=bw, color=(col, 0.8))
		end
	end

	Legend(fig_3d[1, 2],
		[LineElement(color=c, linewidth=3) for c in [:crimson, :dodgerblue, :limegreen]],
		["LAD ($(forest_3d.trees["LAD"].segments.n))",
		 "LCX ($(forest_3d.trees["LCX"].segments.n))",
		 "RCA ($(forest_3d.trees["RCA"].segments.n))"],
		labelcolor=:white, framecolor=:grey40, backgroundcolor=:grey20)

	fig_3d
end

# ╔═╡ cbca48e2-3fbb-4b6e-9f87-2998578795f1
md"""
## Step 2: Kassab Subdivision to Capillaries

The CCO skeleton captures spatial layout of major vessels. Now we subdivide each terminal down to **8 um capillaries** using Kassab's connectivity matrix -- this encodes how many daughters of each Strahler order every parent has. No intersection checking needed; pure statistical branching.
"""

# ╔═╡ a2000003-b000-4c00-8d00-e00000000005
begin
	configs_full = [TreeConfig("RCA", root_pos, 1.5, (0.0, 0.0, -1.0), 100, 1.0)]
	t_gen = @elapsed forest_full = generate_kassab_coronary(
		domain, params_rca;
		tree_configs=configs_full,
		rng=MersenneTwister(42),
		verbose=false,
	)
	tree_full = forest_full.trees["RCA"]
	n_full = tree_full.segments.n
	radii_full = tree_full.segments.radius[1:n_full]
	dmin_um = round(minimum(radii_full) * 2 * 1000, digits=2)
	dmax_um = round(maximum(radii_full) * 2 * 1000, digits=0)

	md"""
	**Full pipeline result:** $(n_full) segments in $(round(t_gen, digits=1))s

	| Metric | Value |
	|:-------|:------|
	| Segments | **$(n_full)** |
	| Min diameter | **$(dmin_um) um** |
	| Max diameter | **$(dmax_um) um** |
	| Diameter range | **$(round(dmax_um / dmin_um, digits=0))x** |
	| Murray exponent | $(params_rca.gamma) |
	"""
end

# ╔═╡ 4200ca7d-d92a-44e3-a145-e4689c1d0435
md"""
## Before vs After Subdivision

Left: CCO skeleton only (100 terminals, ~200 segments). Right: after Kassab subdivision (~30K segments down to capillaries).
"""

# ╔═╡ a2000004-b000-4c00-8d00-e00000000006
begin
	fig_compare = Figure(size=(1200, 500), backgroundcolor=:grey10)

	ax_before = Axis(fig_compare[1, 1],
		title="CCO Skeleton ($(tree_100.segments.n) segments)",
		xlabel="mm", ylabel="mm",
		backgroundcolor=:grey10,
		titlecolor=:white, xlabelcolor=:white, ylabelcolor=:white,
		xticklabelcolor=:grey70, yticklabelcolor=:grey70,
		aspect=DataAspect())
	plot_vessels!(ax_before, tree_100.segments, tree_100.segments.n; color=:crimson)

	ax_after = Axis(fig_compare[1, 2],
		title="After Subdivision ($(n_full) segments)",
		xlabel="mm", ylabel="mm",
		backgroundcolor=:grey10,
		titlecolor=:white, xlabelcolor=:white, ylabelcolor=:white,
		xticklabelcolor=:grey70, yticklabelcolor=:grey70,
		aspect=DataAspect())
	plot_vessels!(ax_after, tree_full.segments, n_full; color=:crimson)

	fig_compare
end

# ╔═╡ 13855190-c1a3-4c1e-803b-e341b627d6cb
md"""
## Color by Vessel Diameter

Visualize the full range from ~3mm stems (bright) down to capillaries (dark). The color scale spans 3 orders of magnitude.
"""

# ╔═╡ a2000005-b000-4c00-8d00-e00000000007
begin
	seg_full = tree_full.segments

	# Auto-select projection axes (two with most extent)
	_xr = extrema(seg_full.proximal_x[i] for i in 1:n_full)
	_yr = extrema(seg_full.proximal_y[i] for i in 1:n_full)
	_zr = extrema(seg_full.proximal_z[i] for i in 1:n_full)
	_spans = [(_xr[2]-_xr[1], :x), (_yr[2]-_yr[1], :y), (_zr[2]-_zr[1], :z)]
	sort!(_spans, by=first, rev=true)
	_ax1, _ax2 = _spans[1][2], _spans[2][2]
	_get(field, i, prox) = field == :x ? (prox ? seg_full.proximal_x[i] : seg_full.distal_x[i]) :
	                       field == :y ? (prox ? seg_full.proximal_y[i] : seg_full.distal_y[i]) :
	                                     (prox ? seg_full.proximal_z[i] : seg_full.distal_z[i])

	widths_full = [clamp(Float32(sqrt(seg_full.radius[i]) * 8), 0.2f0, 8.0f0) for i in 1:n_full]
	colors_full = [log10(max(seg_full.radius[i] * 2000, 0.1)) for i in 1:n_full]
	wmin_f, wmax_f = extrema(widths_full)
	wmax_f == wmin_f && (wmax_f = wmin_f + 1f0)

	fig_diam = Figure(size=(950, 700), backgroundcolor=:grey10)
	ax_diam = Axis(fig_diam[1, 1],
		title="RCA -- Vessel Diameter (log scale)",
		xlabel="mm", ylabel="mm",
		backgroundcolor=:grey10,
		titlecolor=:white, xlabelcolor=:white, ylabelcolor=:white,
		xticklabelcolor=:grey70, yticklabelcolor=:grey70,
		aspect=DataAspect())

	cmap = Makie.to_colormap(:inferno)
	clo, chi = log10(1), log10(4000)

	nbins_d = 20
	for b in 1:nbins_d
		lo = wmin_f + (b - 1) * (wmax_f - wmin_f) / nbins_d
		hi = wmin_f + b * (wmax_f - wmin_f) / nbins_d
		bin_w = (lo + hi) / 2

		pts = Point2f[]
		cols = RGBAf[]
		for i in 1:n_full
			w = widths_full[i]
			if (b == nbins_d ? w >= lo : lo <= w < hi)
				push!(pts, Point2f(_get(_ax1, i, true), _get(_ax2, i, true)))
				push!(pts, Point2f(_get(_ax1, i, false), _get(_ax2, i, false)))
				t = clamp((colors_full[i] - clo) / (chi - clo), 0, 1)
				idx = clamp(round(Int, t * (length(cmap) - 1)) + 1, 1, length(cmap))
				c = cmap[idx]
				push!(cols, c)
				push!(cols, c)
			end
		end
		isempty(pts) && continue
		linesegments!(ax_diam, pts; linewidth=bin_w, color=cols)
	end

	Colorbar(fig_diam[1, 2], limits=(clo, chi),
		colormap=:inferno, label="log10(diameter um)",
		labelcolor=:white, ticklabelcolor=:grey70)

	fig_diam
end

# ╔═╡ a9b70262-9607-4f8e-9e80-d64b4b8eb580
md"""
## Diameter Distribution

The histogram spans from capillaries (~8 um) to the main stem (~3000 um). Dashed lines show Kassab's Strahler order boundaries.
"""

# ╔═╡ 120fa235-0791-469d-ae41-ae47c64b5f41
begin
	all_diameters = [seg_full.radius[i] * 2000 for i in 1:n_full]

	fig_hist = Figure(size=(800, 400))
	ax_hist = Axis(fig_hist[1, 1],
		title="Vessel Diameter Distribution ($(n_full) segments)",
		xlabel="Diameter (um)", ylabel="Count",
		xscale=log10,
		xticks=[1, 3, 10, 30, 100, 300, 1000, 3000])

	hist!(ax_hist, filter(d -> d > 0, all_diameters);
		bins=10.0 .^ range(-0.5, 3.8, length=50),
		color=(:crimson, 0.7))

	for bound in params_rca.diameter_bounds[2:end-1]
		vlines!(ax_hist, [bound]; color=:grey60, linestyle=:dash, linewidth=0.7)
	end

	fig_hist
end

# ╔═╡ a2000006-b000-4c00-8d00-e00000000008
md"""
## Per-Order Statistics

Segment counts and mean diameters by Strahler order, compared to Kassab 1993 reference values.
"""

# ╔═╡ a2000007-b000-4c00-8d00-e00000000009
begin
	orders = tree_full.topology.strahler_order[1:n_full]
	order_stats = []
	for o in sort(unique(orders))
		o < 0 && continue
		idxs = findall(==(o), orders)
		ds = [seg_full.radius[i] * 2000 for i in idxs]
		ref_d = o + 1 <= length(params_rca.diameter_mean) ? params_rca.diameter_mean[o+1] : NaN
		push!(order_stats, (
			order=o,
			n_seg=length(idxs),
			mean_d=round(mean(ds), digits=1),
			ref_d=round(ref_d, digits=1),
		))
	end

	header = "| Order | Segments | Mean D (um) | Kassab D (um) |\n|:------|:---------|:------------|:--------------|\n"
	rows = join([
		"| $(s.order) | $(s.n_seg) | $(s.mean_d) | $(s.ref_d) |"
		for s in order_stats
	], "\n")

	Markdown.parse(header * rows)
end

# ╔═╡ ce6d7a19-f1a9-41cd-8f1c-c510d192bf44
md"""
## Kassab Validation Report Card

9-metric element-level validation against Kassab 1993 morphometric data.
"""

# ╔═╡ 478158cf-3a47-4ae1-934c-fdc7c5972f5a
let
	card = generate_report_card(tree_full, params_rca)
	buf = IOBuffer()
	print_report_card(buf, card)
	@info String(take!(buf))
end

# ╔═╡ Cell order:
# ╟─80b575e3-d4f9-4250-904f-19e88f059029
# ╠═a1000001-b000-4c00-8d00-e00000000002
# ╠═f4c0b7fb-637d-4bd7-adff-4a06fe91c081
# ╠═4378fe0a-1dd5-4b66-b6ae-b0668e12b944
# ╠═fb9ffd67-0ab9-4974-ba52-d64742a58a08
# ╟─dcbc6e0f-c3c5-4ef1-86c4-64ca5a5bc850
# ╠═31d6ee6d-896a-4dab-a0fe-98e551461314
# ╟─9201e42e-a20a-494b-a34b-316f0c536098
# ╠═3c536f7a-9713-40f5-bac9-645cb3998d1a
# ╠═a2000001-b000-4c00-8d00-e00000000003
# ╟─a2000002-b000-4c00-8d00-e00000000004
# ╟─b3000001-c000-4d00-9e00-f00000000001
# ╠═b3000002-c000-4d00-9e00-f00000000002
# ╟─b3000003-c000-4d00-9e00-f00000000003
# ╟─cbca48e2-3fbb-4b6e-9f87-2998578795f1
# ╠═a2000003-b000-4c00-8d00-e00000000005
# ╟─4200ca7d-d92a-44e3-a145-e4689c1d0435
# ╟─a2000004-b000-4c00-8d00-e00000000006
# ╟─13855190-c1a3-4c1e-803b-e341b627d6cb
# ╟─a2000005-b000-4c00-8d00-e00000000007
# ╟─a9b70262-9607-4f8e-9e80-d64b4b8eb580
# ╟─120fa235-0791-469d-ae41-ae47c64b5f41
# ╟─a2000006-b000-4c00-8d00-e00000000008
# ╟─a2000007-b000-4c00-8d00-e00000000009
# ╟─ce6d7a19-f1a9-41cd-8f1c-c510d192bf44
# ╠═478158cf-3a47-4ae1-934c-fdc7c5972f5a
