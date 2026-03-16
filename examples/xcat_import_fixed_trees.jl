using CairoMakie
using StaticArrays
using LinearAlgebra

include(joinpath(dirname(@__DIR__), "src", "types.jl"))
include(joinpath(dirname(@__DIR__), "src", "xcat_nrb.jl"))
include(joinpath(dirname(@__DIR__), "src", "xcat_integration.jl"))

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function plot_vascular_tree!(ax, tree::VascularTree; color=:steelblue, linewidth=2.5)
    seg = tree.segments
    for i in 1:seg.n
        lines!(
            ax,
            [seg.proximal_x[i], seg.distal_x[i]],
            [seg.proximal_y[i], seg.distal_y[i]],
            [seg.proximal_z[i], seg.distal_z[i]];
            color=color,
            linewidth=linewidth,
        )
    end
    root_id = tree.root_segment_id
    if root_id > 0
        text!(
            ax,
            seg.proximal_x[root_id],
            seg.proximal_y[root_id],
            seg.proximal_z[root_id];
            text=tree.name,
            color=color,
            fontsize=15,
        )
    end
    return ax
end

function make_imported_tree_figure(trees::Dict{String, VascularTree})
    fig = Figure(size=(1600, 1200))
    ax = Axis3(
        fig[1, 1];
        title="XCAT proximal trees imported into VascularTree",
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        aspect=:data,
    )
    palette = Dict(
        "LAD" => :forestgreen,
        "LCX" => :mediumpurple,
        "RCA" => :crimson,
    )
    for name in ["LAD", "LCX", "RCA"]
        haskey(trees, name) || continue
        plot_vascular_tree!(ax, trees[name]; color=palette[name], linewidth=3.0)
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
    imported = xcat_trees_to_vascular_trees(trees)

    println("Imported XCAT proximal trees into VascularTree")
    println()
    for name in ["LAD", "LCX", "RCA"]
        tree = imported[name]
        println(
            rpad(name, 6),
            " segments=", lpad(tree.segments.n, 4),
            " terminals=", lpad(tree.n_terminals, 3),
            " bifurcations=", lpad(tree.n_bifurcations, 3),
            " root_id=", tree.root_segment_id,
        )
    end

    fig = make_imported_tree_figure(imported)
    fig_path = joinpath(output_dir, "xcat_imported_fixed_trees.png")
    save(fig_path, fig)
    println()
    println("Saved imported-tree figure to: $fig_path")
end

main()
