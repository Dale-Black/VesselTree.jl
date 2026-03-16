using CairoMakie
using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function plot_tree_segments!(
    ax,
    tree::VascularTree;
    color=:steelblue,
    linewidth=2.0,
    alpha=0.9,
    min_radius_mm=0.0,
)
    seg = tree.segments
    for i in 1:seg.n
        seg.radius[i] < min_radius_mm && continue
        lines!(
            ax,
            [seg.proximal_x[i], seg.distal_x[i]],
            [seg.proximal_y[i], seg.distal_y[i]],
            [seg.proximal_z[i], seg.distal_z[i]];
            color=(color, alpha),
            linewidth=linewidth,
        )
    end
    return ax
end

function make_overview_figure(forest::CoronaryForest)
    fig = Figure(size=(1800, 1400))
    palette = Dict(
        "LAD" => :forestgreen,
        "LCX" => :mediumpurple,
        "RCA" => :crimson,
    )

    ax_full = Axis3(
        fig[1, 1];
        title="XCAT Kassab coronary forest",
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        aspect=:data,
    )
    ax_macro = Axis3(
        fig[1, 2];
        title="Macro view (radius >= 0.15 mm)",
        xlabel="X (mm)",
        ylabel="Y (mm)",
        zlabel="Z (mm)",
        aspect=:data,
    )

    for name in ["LAD", "LCX", "RCA"]
        tree = forest.trees[name]
        plot_tree_segments!(ax_full, tree; color=palette[name], linewidth=1.2, alpha=0.4)
        plot_tree_segments!(ax_macro, tree; color=palette[name], linewidth=2.6, alpha=0.95, min_radius_mm=0.15)
    end

    return fig
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(@__DIR__), "output")
    mkpath(output_dir)

    rng = MersenneTwister(123)
    params = kassab_coronary_params()
    result = generate_xcat_kassab_coronary(
        nrb_path,
        params;
        phase="dias",
        rng=rng,
        verbose=true,
        handoff_order=6,
        subdivision_max_order=6,
        target_interior=80_000,
        max_candidates=500_000,
        batch_size=8192,
        tree_capacity=20_000,
        target_terminals=Dict(
            "LAD" => 140,
            "LCX" => 90,
            "RCA" => 110,
        ),
        territory_fractions=Dict(
            "LAD" => 0.40,
            "LCX" => 0.25,
            "RCA" => 0.35,
        ),
    )

    println("Validation summary")
    for name in ["LAD", "LCX", "RCA"]
        tree = result.forest.trees[name]
        println(
            rpad(name, 6),
            " segments=", tree.segments.n,
            " terminals=", tree.n_terminals,
            " bifurcations=", tree.n_bifurcations,
        )
    end

    fig = make_overview_figure(result.forest)
    fig_path = joinpath(output_dir, "xcat_kassab_coronary.png")
    save(fig_path, fig)
    println("Saved XCAT Kassab figure to: $fig_path")
end

main()
