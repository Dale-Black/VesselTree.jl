using CairoMakie
using StaticArrays

include(joinpath(dirname(@__DIR__), "src", "xcat_nrb.jl"))

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function plot_control_net!(ax, surface::XCATNurbsSurface; color=:steelblue, linewidth=1.0, alpha=0.7)
    pts = surface.control_points
    nu, nv = size(pts)

    for u in 1:nu
        row = pts[u, :]
        lines!(
            ax,
            [p[1] for p in row],
            [p[2] for p in row],
            [p[3] for p in row];
            color=(color, alpha),
            linewidth=linewidth,
        )
    end

    for v in 1:nv
        col = pts[:, v]
        lines!(
            ax,
            [p[1] for p in col],
            [p[2] for p in col],
            [p[3] for p in col];
            color=(color, alpha),
            linewidth=linewidth,
        )
    end

    c = xcat_center(surface)
    text!(ax, c[1], c[2], c[3]; text=surface.name, color=color, fontsize=12)
    return ax
end

function make_xcat_inspection_figure(
    surfaces::AbstractVector{XCATNurbsSurface};
    title::AbstractString="XCAT heart/coronary control-net inspection",
)
    fig = Figure(size=(1600, 1200))
    ax = Axis3(
        fig[1, 1];
        title=title,
        xlabel="X",
        ylabel="Y",
        zlabel="Z",
        aspect=:data,
    )

    palette = [
        :black,
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
        linewidth = occursin("pericardium", lowercase(surface.name)) ? 0.8 : 1.5
        alpha = occursin("pericardium", lowercase(surface.name)) ? 0.28 : 0.8
        plot_control_net!(ax, surface; color=color, linewidth=linewidth, alpha=alpha)
    end

    fig
end

function print_surface_summary(surfaces::AbstractVector{XCATNurbsSurface})
    println("Parsed $(length(surfaces)) XCAT objects")
    println()
    for row in xcat_summary_rows(surfaces)
        println(
            rpad(row.name, 20),
            " size=(", row.nu, ", ", row.nv, ")",
            " bounds=[(",
            round(row.min_x; digits=2), ", ",
            round(row.min_y; digits=2), ", ",
            round(row.min_z; digits=2), ") -> (",
            round(row.max_x; digits=2), ", ",
            round(row.max_y; digits=2), ", ",
            round(row.max_z; digits=2), ")]",
        )
    end
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    output_path = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(@__DIR__), "output", "xcat_nrb_inspect.png")

    surfaces = parse_xcat_nrb(nrb_path)

    selected = xcat_select_objects(
        surfaces;
        names=[
            "dias_pericardium",
            "dias_aorta",
            "dias_lad1",
            "dias_lad2",
            "dias_lad3",
            "dias_lcx",
            "dias_rca1",
            "dias_rca2",
        ],
    )

    print_surface_summary(selected)

    fig = make_xcat_inspection_figure(selected)
    save(output_path, fig)
    println()
    println("Saved inspection figure to: $output_path")
end

main()
