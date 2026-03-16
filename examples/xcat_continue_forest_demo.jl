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
    segment_range=1:tree.segments.n,
)
    seg = tree.segments
    for i in segment_range
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

function make_forest_figure(forest::CoronaryForest, initial_segment_counts::Dict{String, Int})
    fig = Figure(size=(1700, 1300))
    ax = Axis3(
        fig[1, 1];
        title="XCAT fixed-tree continuation demo\nbold = imported trunk, thin = continued growth",
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
        tree = forest.trees[name]
        initial_n = get(initial_segment_counts, name, 0)
        if initial_n > 0
            plot_tree_segments!(
                ax,
                tree;
                color=palette[name],
                linewidth=3.6,
                alpha=0.95,
                segment_range=1:min(initial_n, tree.segments.n),
            )
        end
        if tree.segments.n > initial_n
            plot_tree_segments!(
                ax,
                tree;
                color=palette[name],
                linewidth=1.4,
                alpha=0.55,
                segment_range=(initial_n + 1):tree.segments.n,
            )
        end
        rid = tree.root_segment_id
        seg = tree.segments
        text!(
            ax,
            seg.proximal_x[rid],
            seg.proximal_y[rid],
            seg.proximal_z[rid];
            text=name,
            color=palette[name],
            fontsize=15,
        )
    end

    return fig
end

function main()
    nrb_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_XCAT_NRB
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(@__DIR__), "output")
    mkpath(output_dir)

    rng = MersenneTwister(42)
    targets = Dict(
        "LAD" => 10,
        "LCX" => 8,
        "RCA" => 8,
    )
    params = kassab_coronary_params()
    result = generate_xcat_coronary_forest(
        nrb_path,
        params;
        rng=rng,
        verbose=true,
        kassab=true,
        target_interior=80_000,
        max_candidates=500_000,
        batch_size=8192,
        tree_capacity=5000,
        target_terminals=targets,
        territory_fractions=Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35),
    )
    forest = result.forest
    initial_segment_counts = Dict(name => tree.segments.n for (name, tree) in result.fixed_trees)
    initial_terminal_counts = Dict(name => tree.n_terminals for (name, tree) in result.fixed_trees)

    println("Continuation summary")
    for name in ["LAD", "LCX", "RCA"]
        tree = forest.trees[name]
        println(
            rpad(name, 6),
            " terminals=", tree.n_terminals,
            " (Δ", tree.n_terminals - initial_terminal_counts[name], ")",
            " segments=", tree.segments.n,
            " (Δ", tree.segments.n - initial_segment_counts[name], ")",
            " bifurcations=", tree.n_bifurcations,
        )
    end

    fig = make_forest_figure(forest, initial_segment_counts)
    fig_path = joinpath(output_dir, "xcat_continue_forest_demo.png")
    save(fig_path, fig)
    println("Saved continuation demo figure to: $fig_path")
end

main()
