using CairoMakie
using StaticArrays
using LinearAlgebra
using Statistics

include(joinpath(dirname(@__DIR__), "src", "xcat_nrb.jl"))

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function plot_sampled_surface_grid!(
    ax,
    surface::XCATNurbsSurface;
    color=:steelblue,
    linewidth=1.0,
    alpha=0.35,
    n_u::Int=48,
    n_v::Int=24,
)
    points, _, _, _ = xcat_sample_surface(surface; n_u=n_u, n_v=n_v, orient_outward=false)
    for j in 1:size(points, 1)
        row = points[j, :]
        lines!(
            ax,
            [p[1] for p in row],
            [p[2] for p in row],
            [p[3] for p in row];
            color=(color, alpha),
            linewidth=linewidth,
        )
    end
    for i in 1:size(points, 2)
        col = points[:, i]
        lines!(
            ax,
            [p[1] for p in col],
            [p[2] for p in col],
            [p[3] for p in col];
            color=(color, alpha),
            linewidth=linewidth,
        )
    end
    return ax
end

function plot_centerline!(ax, centerline::XCATCenterline; color=:red, linewidth=3.0)
    pts = centerline.centers
    lines!(
        ax,
        [p[1] for p in pts],
        [p[2] for p in pts],
        [p[3] for p in pts];
        color=color,
        linewidth=linewidth,
    )
    start_p = first(pts)
    text!(ax, start_p[1], start_p[2], start_p[3]; text=centerline.name, color=color, fontsize=13)
    return ax
end

function make_centerline_figure(
    surfaces::AbstractVector{XCATNurbsSurface},
    centerlines::AbstractVector{XCATCenterline},
)
    fig = Figure(size=(1800, 1300))
    ax = Axis3(
        fig[1, 1];
        title="XCAT coronary surfaces and extracted centerlines",
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        aspect=:data,
    )

    palette = [
        :crimson,
        :darkorange,
        :goldenrod,
        :forestgreen,
        :dodgerblue,
        :mediumpurple,
        :deeppink,
    ]

    for (idx, surface) in enumerate(surfaces)
        color = palette[mod1(idx, length(palette))]
        plot_sampled_surface_grid!(ax, surface; color=color, linewidth=0.9, alpha=0.3)
    end

    for (idx, centerline) in enumerate(centerlines)
        color = palette[mod1(idx, length(palette))]
        plot_centerline!(ax, centerline; color=color, linewidth=3.0)
    end

    return fig
end

function make_tree_figure(trees::Dict{String, XCATCenterlineTree})
    fig = Figure(size=(1600, 1200))
    ax = Axis3(
        fig[1, 1];
        title="XCAT proximal coronary trees",
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        aspect=:data,
    )

    palette = Dict(
        "AORTA" => :deeppink,
        "LAD" => :forestgreen,
        "LCX" => :mediumpurple,
        "RCA" => :crimson,
    )

    for name in ["AORTA", "LAD", "LCX", "RCA"]
        haskey(trees, name) || continue
        tree = trees[name]
        color = palette[name]
        for seg in values(tree.segments)
            pts = seg.centers
            lines!(
                ax,
                [p[1] for p in pts],
                [p[2] for p in pts],
                [p[3] for p in pts];
                color=color,
                linewidth=name == "AORTA" ? 5.0 : 4.0,
            )
        end
        root_pts = tree.segments[tree.root_segment].centers
        start_p = first(root_pts)
        text!(ax, start_p[1], start_p[2], start_p[3]; text=name, color=color, fontsize=15)

        for conn in tree.connections
            parent_seg = tree.segments[conn.parent_segment]
            child_seg = tree.segments[conn.child_segment]
            p = parent_seg.centers[conn.parent_index]
            q = child_seg.centers[conn.child_index]
            scatter!(ax, [p[1]], [p[2]], [p[3]]; color=:black, markersize=10)
            lines!(ax, [p[1], q[1]], [p[2], q[2]], [p[3], q[3]]; color=(color, 0.45), linewidth=1.5, linestyle=:dash)
        end
    end

    return fig
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(@__DIR__), "output")
    mkpath(output_dir)

    surfaces = parse_xcat_nrb(nrb_path)
    selected = xcat_select_objects(
        surfaces;
        names=[
            "dias_aorta",
            "dias_lad1",
            "dias_lad2",
            "dias_lad3",
            "dias_lcx",
            "dias_rca1",
            "dias_rca2",
        ],
    )

    centerlines = [xcat_centerline_from_surface(surface) for surface in selected]
    trees = xcat_build_coronary_trees(centerlines)

    println("Extracted $(length(centerlines)) centerlines")
    println()
    for row in xcat_centerline_summary_rows(centerlines)
        println(
            rpad(row.name, 16),
            " n=", lpad(row.n, 3),
            " axial=", row.axial,
            " radius_mean=", round(row.radius_mean; digits=2),
            " mm range=[", round(row.radius_min; digits=2), ", ", round(row.radius_max; digits=2), "]",
        )
    end

    for line in centerlines
        out_path = joinpath(output_dir, "$(line.name)_centerline.csv")
        xcat_export_centerline_csv(line, out_path)
    end

    println()
    println("Coronary trees")
    for row in xcat_tree_summary_rows([trees["LAD"], trees["LCX"], trees["RCA"]])
        println(
            rpad(row.name, 6),
            " segments=", row.n_segments,
            " points=", row.n_points,
            " mean_radius=", round(row.mean_radius; digits=2),
            " mm max_gap=", round(row.max_gap; digits=2),
            " mm",
        )
    end

    for (tree_name, tree) in pairs(trees)
        tree_name == "AORTA" && continue
        for (seg_name, seg) in pairs(tree.segments)
            out_path = joinpath(output_dir, "$(tree_name)_$(seg_name)_segment.csv")
            xcat_export_centerline_csv(seg, out_path)
        end
    end

    fig = make_centerline_figure(selected, centerlines)
    fig_path = joinpath(output_dir, "xcat_coronary_centerlines.png")
    save(fig_path, fig)
    tree_fig = make_tree_figure(trees)
    tree_fig_path = joinpath(output_dir, "xcat_coronary_trees.png")
    save(tree_fig_path, tree_fig)
    println()
    println("Saved coronary centerline figure to: $fig_path")
    println("Saved coronary tree figure to: $tree_fig_path")
end

main()
