using CairoMakie
using StaticArrays
using LinearAlgebra
using Random

include(joinpath(dirname(@__DIR__), "src", "domain.jl"))
include(joinpath(dirname(@__DIR__), "src", "xcat_nrb.jl"))

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function plot_sampled_surface_grid!(
    ax,
    surface::XCATNurbsSurface;
    color=:steelblue,
    linewidth=1.0,
    alpha=0.5,
    n_u::Int=50,
    n_v::Int=24,
)
    points, _, _, _ = xcat_sample_surface(surface; n_u=n_u, n_v=n_v, orient_outward=true)

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

function make_shell_domain_figure(
    domain::CSVShellDomain,
    surfaces::AbstractVector{XCATNurbsSurface};
    n_interior::Int=6000,
)
    object_map = xcat_object_dict(surfaces)
    fig = Figure(size=(1800, 1300))
    ax = Axis3(
        fig[1, 1];
        title="XCAT myocardial shell domain inspection",
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        aspect=:data,
    )

    plot_sampled_surface_grid!(ax, object_map["dias_pericardium"]; color=:black, linewidth=0.8, alpha=0.18, n_u=80, n_v=40)

    cavity_palette = [:crimson, :darkorange, :goldenrod, :forestgreen, :dodgerblue, :mediumpurple, :deeppink, :brown]
    for (idx, cavity_name) in enumerate(domain.cavity_names)
        plot_sampled_surface_grid!(
            ax,
            object_map[cavity_name];
            color=cavity_palette[mod1(idx, length(cavity_palette))],
            linewidth=0.9,
            alpha=0.45,
            n_u=42,
            n_v=22,
        )
    end

    total = size(domain.interior_points, 1)
    keep = min(n_interior, total)
    sample_idx = randperm(total)[1:keep]
    scatter!(
        ax,
        domain.interior_points[sample_idx, 1],
        domain.interior_points[sample_idx, 2],
        domain.interior_points[sample_idx, 3];
        markersize=2.0,
        color=(:slategray, 0.18),
    )

    fig
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(@__DIR__), "output")

    mkpath(output_dir)
    domain, surfaces = xcat_myocardial_shell_domain(
        nrb_path;
        phase="dias",
        rng=MersenneTwister(42),
        target_interior=60_000,
        max_candidates=240_000,
        batch_size=8192,
    )

    println("Built XCAT shell domain")
    println("  center = $(domain.center)")
    println("  min_corner = $(domain.min_corner)")
    println("  max_corner = $(domain.max_corner)")
    println("  interior samples = $(size(domain.interior_points, 1))")
    println("  estimated volume = $(round(domain.volume; digits=2)) mm^3")
    println("  cavities = $(join(domain.cavity_names, ", "))")

    fig = make_shell_domain_figure(domain, surfaces)
    fig_path = joinpath(output_dir, "xcat_myocardial_shell_domain.png")
    save(fig_path, fig)
    println("Saved shell-domain inspection figure to: $fig_path")
end

main()
