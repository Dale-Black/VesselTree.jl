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
				push!(pts, Point2f(seg.proximal_x[i], seg.proximal_y[i]))
				push!(pts, Point2f(seg.distal_x[i], seg.distal_y[i]))
			end
		end
		isempty(pts) && continue
		linesegments!(ax, pts; linewidth=bin_w, color=(color, alpha))
	end
end

# ╔═╡ a2000001-b000-4c00-8d00-e00000000003
begin
	# Grow CCO skeleton at 3 stages
	tree_10 = VascularTree("RCA", 5000)
	VesselTree.add_segment!(tree_10, (0.0, 0.0, 50.0), (0.0, 0.0, 45.0), 1.5, Int32(-1))
	grow_tree!(tree_10, domain, 10, params_rca; rng=MersenneTwister(42))

	tree_30 = VascularTree("RCA", 5000)
	VesselTree.add_segment!(tree_30, (0.0, 0.0, 50.0), (0.0, 0.0, 45.0), 1.5, Int32(-1))
	grow_tree!(tree_30, domain, 30, params_rca; rng=MersenneTwister(42))

	tree_100 = VascularTree("RCA", 5000)
	VesselTree.add_segment!(tree_100, (0.0, 0.0, 50.0), (0.0, 0.0, 45.0), 1.5, Int32(-1))
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
			xlabel="x (mm)", ylabel="y (mm)",
			backgroundcolor=:grey10,
			titlecolor=:white, xlabelcolor=:white, ylabelcolor=:white,
			xticklabelcolor=:grey70, yticklabelcolor=:grey70,
			aspect=DataAspect())
		plot_vessels!(ax, tree_s.segments, tree_s.segments.n; color=:limegreen)
	end

	fig_cco
end

# ╔═╡ cbca48e2-3fbb-4b6e-9f87-2998578795f1
md"""
## Step 2: Kassab Subdivision to Capillaries

The CCO skeleton captures spatial layout of major vessels. Now we subdivide each terminal down to **8 um capillaries** using Kassab's connectivity matrix -- this encodes how many daughters of each Strahler order every parent has. No intersection checking needed; pure statistical branching.
"""

# ╔═╡ a2000003-b000-4c00-8d00-e00000000005
begin
	configs_full = [TreeConfig("RCA", (0.0, 0.0, 50.0), 1.5, (0.0, 0.0, -1.0), 100, 1.0)]
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
		xlabel="x (mm)", ylabel="y (mm)",
		backgroundcolor=:grey10,
		titlecolor=:white, xlabelcolor=:white, ylabelcolor=:white,
		xticklabelcolor=:grey70, yticklabelcolor=:grey70,
		aspect=DataAspect())
	plot_vessels!(ax_before, tree_100.segments, tree_100.segments.n; color=:crimson)

	ax_after = Axis(fig_compare[1, 2],
		title="After Subdivision ($(n_full) segments)",
		xlabel="x (mm)", ylabel="y (mm)",
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
	widths_full = [clamp(Float32(sqrt(seg_full.radius[i]) * 8), 0.2f0, 8.0f0) for i in 1:n_full]
	colors_full = [log10(max(seg_full.radius[i] * 2000, 0.1)) for i in 1:n_full]
	wmin_f, wmax_f = extrema(widths_full)
	wmax_f == wmin_f && (wmax_f = wmin_f + 1f0)

	fig_diam = Figure(size=(950, 700), backgroundcolor=:grey10)
	ax_diam = Axis(fig_diam[1, 1],
		title="RCA -- Vessel Diameter (log scale)",
		xlabel="x (mm)", ylabel="y (mm)",
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
				push!(pts, Point2f(seg_full.proximal_x[i], seg_full.proximal_y[i]))
				push!(pts, Point2f(seg_full.distal_x[i], seg_full.distal_y[i]))
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
	Text(String(take!(buf)))
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
