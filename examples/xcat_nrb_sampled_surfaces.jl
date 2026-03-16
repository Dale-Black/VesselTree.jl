using CairoMakie
using StaticArrays
using LinearAlgebra

include(joinpath(dirname(@__DIR__), "src", "xcat_nrb.jl"))

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function plot_sampled_surface_grid!(
    ax,
    surface::XCATNurbsSurface;
    color=:steelblue,
    linewidth=1.0,
    alpha=0.8,
    n_u::Int=60,
    n_v::Int=24,
)
    points, _, _, _ = xcat_sample_surface(surface; n_u=n_u, n_v=n_v, orient_outward=occursin("pericardium", lowercase(surface.name)))

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

    center = let acc = SVector(0.0, 0.0, 0.0)
        for p in points
            acc += p
        end
        acc / length(points)
    end
    text!(ax, center[1], center[2], center[3]; text=surface.name, color=color, fontsize=12)
    return ax
end

function make_xcat_surface_figure(surfaces::AbstractVector{XCATNurbsSurface})
    fig = Figure(size=(1800, 1300))
    ax = Axis3(
        fig[1, 1];
        title="XCAT sampled NURBS surfaces",
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
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
        is_outer = occursin("pericardium", lowercase(surface.name))
        plot_sampled_surface_grid!(
            ax,
            surface;
            color=color,
            linewidth=is_outer ? 0.8 : 1.3,
            alpha=is_outer ? 0.22 : 0.85,
            n_u=is_outer ? 80 : 60,
            n_v=is_outer ? 40 : 24,
        )
    end

    return fig
end

function print_surface_summary(surfaces::AbstractVector{XCATNurbsSurface})
    println("Sampled $(length(surfaces)) XCAT objects")
    println()
    for row in xcat_summary_rows(surfaces)
        println(
            rpad(row.name, 20),
            " uv=(", row.nu, ", ", row.nv, ")",
            " degree=(", row.pu, ", ", row.pv, ")",
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
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(@__DIR__), "output")

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

    mkpath(output_dir)
    fig = make_xcat_surface_figure(selected)
    fig_path = joinpath(output_dir, "xcat_nrb_sampled_surfaces.png")
    save(fig_path, fig)

    object_map = xcat_object_dict(surfaces)
    pericardium = object_map["dias_pericardium"]
    points_path = joinpath(output_dir, "xcat_dias_pericardium_points.csv")
    normals_path = joinpath(output_dir, "xcat_dias_pericardium_normals.csv")
    xcat_export_sampled_surface_csv(
        pericardium,
        points_path,
        normals_path;
        n_u=160,
        n_v=96,
        orient_outward=true,
    )

    println()
    println("Saved sampled-surface inspection figure to: $fig_path")
    println("Saved sampled pericardium points to: $points_path")
    println("Saved sampled pericardium normals to: $normals_path")
end

main()
