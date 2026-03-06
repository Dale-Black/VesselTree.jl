module VesselTreeMakieExt

using VesselTree
using CairoMakie
using Distributions

"""
    plot_tree(tree; color_by=:order, colormap=:viridis)

3D visualization of tree segments as colored lines.
Line width proportional to log(radius). Color by :order, :radius, :flow, or :generation.
"""
function VesselTree.plot_tree(tree::VascularTree; color_by::Symbol=:order, colormap=:viridis)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    colors = _get_color_values(seg, topo, n, color_by)
    widths = [max(0.5, log(1 + seg.radius[i] * 1000) * 2) for i in 1:n]

    fig = Figure(size=(800, 600))
    ax = Axis3(fig[1, 1], title="VascularTree: $(tree.name)", xlabel="x", ylabel="y", zlabel="z")

    for i in 1:n
        lines!(ax,
            [seg.proximal_x[i], seg.distal_x[i]],
            [seg.proximal_y[i], seg.distal_y[i]],
            [seg.proximal_z[i], seg.distal_z[i]];
            color=colors[i],
            linewidth=widths[i],
            colormap=colormap,
            colorrange=(minimum(colors), maximum(colors) + eps()),
        )
    end

    Colorbar(fig[1, 2], limits=(minimum(colors), maximum(colors) + eps()),
        colormap=colormap, label=String(color_by))

    return fig
end

"""
    plot_tree_2d(tree; projection=:xy, color_by=:order, colormap=:viridis)

2D projection of tree segments.
"""
function VesselTree.plot_tree_2d(tree::VascularTree; projection::Symbol=:xy, color_by::Symbol=:order, colormap=:viridis)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    colors = _get_color_values(seg, topo, n, color_by)
    widths = [max(0.5, log(1 + seg.radius[i] * 1000) * 2) for i in 1:n]

    px1, py1, px2, py2 = _project_coords(seg, n, projection)

    fig = Figure(size=(800, 600))
    xlabel, ylabel = _projection_labels(projection)
    ax = Axis(fig[1, 1], title="VascularTree: $(tree.name) ($(projection))", xlabel=xlabel, ylabel=ylabel, aspect=DataAspect())

    for i in 1:n
        lines!(ax,
            [px1[i], px2[i]],
            [py1[i], py2[i]];
            color=colors[i],
            linewidth=widths[i],
            colormap=colormap,
            colorrange=(minimum(colors), maximum(colors) + eps()),
        )
    end

    Colorbar(fig[1, 2], limits=(minimum(colors), maximum(colors) + eps()),
        colormap=colormap, label=String(color_by))

    return fig
end

"""
    plot_validation_report(tree, params)

Multi-panel validation figure: diameter histogram, asymmetry distribution,
angle histogram, connectivity heatmap.
"""
function VesselTree.plot_validation_report(tree::VascularTree, params::MorphometricParams)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    fig = Figure(size=(1200, 900))

    # Panel 1: Diameter histogram by order
    ax1 = Axis(fig[1, 1], title="Diameter by Strahler Order", xlabel="Diameter (um)", ylabel="Count")
    diameters_um = [seg.radius[i] * 2 * 1000 for i in 1:n]
    orders = [Int(topo.strahler_order[i]) for i in 1:n]
    hist!(ax1, diameters_um; bins=50, color=(:blue, 0.5))

    # Panel 2: Asymmetry ratio distribution
    ax2 = Axis(fig[1, 2], title="Asymmetry Ratio", xlabel="S = r_small/r_large", ylabel="Count")
    asymmetries = Float64[]
    for i in 1:n
        c1 = topo.child1_id[i]; c2 = topo.child2_id[i]
        (c1 <= 0 || c2 <= 0) && continue
        r1, r2 = seg.radius[c1], seg.radius[c2]
        r1 > 0 && r2 > 0 && push!(asymmetries, min(r1, r2) / max(r1, r2))
    end
    if !isempty(asymmetries)
        hist!(ax2, asymmetries; bins=30, color=(:green, 0.5), normalization=:pdf)
        # Overlay Beta(2.5, 0.8) PDF
        xs = range(0.01, 0.99, length=100)
        beta_dist = Beta(VesselTree.ASYMMETRY_ALPHA, VesselTree.ASYMMETRY_BETA)
        lines!(ax2, xs, pdf.(beta_dist, xs), color=:red, linewidth=2, label="Beta(2.5,0.8)")
        axislegend(ax2; position=:lt)
    end

    # Panel 3: Branching angle histogram
    ax3 = Axis(fig[2, 1], title="Branching Angles", xlabel="Angle (degrees)", ylabel="Count")
    angles = Float64[]
    for i in 1:n
        c1 = topo.child1_id[i]; c2 = topo.child2_id[i]
        (c1 <= 0 || c2 <= 0) && continue
        # Parent direction
        pdx = seg.distal_x[i] - seg.proximal_x[i]
        pdy = seg.distal_y[i] - seg.proximal_y[i]
        pdz = seg.distal_z[i] - seg.proximal_z[i]
        plen = sqrt(pdx^2 + pdy^2 + pdz^2)
        plen < 1e-15 && continue
        for ci in (c1, c2)
            ci <= 0 && continue
            cdx = seg.distal_x[ci] - seg.proximal_x[ci]
            cdy = seg.distal_y[ci] - seg.proximal_y[ci]
            cdz = seg.distal_z[ci] - seg.proximal_z[ci]
            clen = sqrt(cdx^2 + cdy^2 + cdz^2)
            clen < 1e-15 && continue
            cosang = (pdx*cdx + pdy*cdy + pdz*cdz) / (plen * clen)
            cosang = clamp(cosang, -1.0, 1.0)
            push!(angles, acos(cosang) * 180.0 / Float64(π))
        end
    end
    if !isempty(angles)
        hist!(ax3, angles; bins=36, color=(:orange, 0.5))
    end

    # Panel 4: Order distribution
    ax4 = Axis(fig[2, 2], title="Segments per Strahler Order", xlabel="Order", ylabel="Count")
    if !isempty(orders)
        max_ord = maximum(orders)
        order_counts = zeros(Int, max_ord + 1)
        for o in orders
            order_counts[o + 1] += 1
        end
        barplot!(ax4, 0:max_ord, order_counts, color=:steelblue)
    end

    return fig
end

"""
    plot_forest(forest; color_by=:tree, colormap=:Set1_4)

Color segments by tree membership or other attribute.
"""
function VesselTree.plot_forest(forest::VesselTree.CoronaryForest; color_by::Symbol=:tree, colormap=:Set1_4)
    fig = Figure(size=(900, 700))
    ax = Axis3(fig[1, 1], title="Coronary Forest", xlabel="x", ylabel="y", zlabel="z")

    tree_names = sort(collect(keys(forest.trees)))
    for (tidx, name) in enumerate(tree_names)
        tree = forest.trees[name]
        n = tree.segments.n
        seg = tree.segments

        for i in 1:n
            lines!(ax,
                [seg.proximal_x[i], seg.distal_x[i]],
                [seg.proximal_y[i], seg.distal_y[i]],
                [seg.proximal_z[i], seg.distal_z[i]];
                color=Makie.wong_colors()[mod1(tidx, 7)],
                linewidth=max(0.5, log(1 + seg.radius[i] * 1000)),
                label=(i == 1 ? name : nothing),
            )
        end
    end

    return fig
end

# --- Internal helpers ---

function _get_color_values(seg, topo, n, color_by)
    if color_by == :order
        return Float64[topo.strahler_order[i] for i in 1:n]
    elseif color_by == :radius
        return Float64[seg.radius[i] for i in 1:n]
    elseif color_by == :flow
        return Float64[seg.flow[i] for i in 1:n]
    elseif color_by == :generation
        return Float64[topo.generation[i] for i in 1:n]
    else
        return Float64[topo.strahler_order[i] for i in 1:n]
    end
end

function _project_coords(seg, n, projection)
    if projection == :xy
        return (seg.proximal_x[1:n], seg.proximal_y[1:n], seg.distal_x[1:n], seg.distal_y[1:n])
    elseif projection == :xz
        return (seg.proximal_x[1:n], seg.proximal_z[1:n], seg.distal_x[1:n], seg.distal_z[1:n])
    elseif projection == :yz
        return (seg.proximal_y[1:n], seg.proximal_z[1:n], seg.distal_y[1:n], seg.distal_z[1:n])
    else
        return (seg.proximal_x[1:n], seg.proximal_y[1:n], seg.distal_x[1:n], seg.distal_y[1:n])
    end
end

function _projection_labels(projection)
    if projection == :xy
        return ("x", "y")
    elseif projection == :xz
        return ("x", "z")
    elseif projection == :yz
        return ("y", "z")
    else
        return ("x", "y")
    end
end

end # module
