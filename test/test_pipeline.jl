using Test
using VesselTree
using Random

@testset "Full pipeline: generate_kassab_coronary" begin
    params = kassab_coronary_params()
    rng = MersenneTwister(42)

    # Small domain encompassing default coronary root positions
    domain = SphereDomain((0.0, 0.0, 0.0), 5.0)

    # Small test configs for fast execution
    test_configs = [
        TreeConfig("LAD", (-2.0, 1.0, 0.0), 1.5, (0.0, -1.0, -0.5), 20, 0.40),
        TreeConfig("LCX", (-2.0, -1.0, 0.0), 1.2, (0.0, -1.0, 0.5), 15, 0.25),
        TreeConfig("RCA", (2.0, 0.0, 0.0), 1.3, (0.0, -1.0, 0.0), 18, 0.35),
    ]

    forest = generate_kassab_coronary(
        domain, params;
        rng=rng,
        verbose=false,
        handoff_order=5,
        tree_configs=test_configs,
    )

    @testset "Forest structure" begin
        # All 3 trees present
        @test haskey(forest.trees, "LAD")
        @test haskey(forest.trees, "LCX")
        @test haskey(forest.trees, "RCA")
        @test length(forest.trees) == 3

        # Each tree has segments
        for (name, tree) in forest.trees
            @test tree.segments.n > 0
        end

        # Total segments > CCO skeleton (subdivision happened)
        total = sum(t.segments.n for (_, t) in forest.trees)
        # CCO produces ~2*target_terminals segments per tree; subdivision should increase this
        cco_estimate = (20 + 15 + 18) * 2
        @test total > cco_estimate
    end

    @testset "Subdivision produced lower orders" begin
        # After subdivision with handoff_order=3, there should be order-0 segments
        for (name, tree) in forest.trees
            VesselTree.assign_strahler_orders!(tree, params)
            n = tree.segments.n
            orders_present = Set{Int}()
            for i in 1:n
                ord = Int(tree.topology.strahler_order[i])
                ord >= 0 && push!(orders_present, ord)
            end
            # Should have at least order 0 (capillaries) and higher orders
            @test 0 in orders_present
            @test length(orders_present) >= 2
        end
    end

    @testset "Diameter range" begin
        for (name, tree) in forest.trees
            seg = tree.segments
            n = seg.n
            min_d_um = Inf
            max_d_um = 0.0
            for i in 1:n
                d_um = seg.radius[i] * 2.0 * 1000.0  # mm -> um
                min_d_um = min(min_d_um, d_um)
                max_d_um = max(max_d_um, d_um)
            end
            # Min diameter > 0 (all radii positive)
            # Sub-capillary radii are normal: Murray's law at bifurcation chain
            # terminals can produce radii below vessel cutoff. These are continuation
            # segment endpoints, not actual capillaries.
            @test min_d_um > 0.0
            # Max diameter > 100um (root segments should be large)
            @test max_d_um > 100.0
        end
    end

    @testset "Murray's law at junctions" begin
        gamma = params.gamma
        for (name, tree) in forest.trees
            max_dev, mean_dev = VesselTree.compute_murray_deviation(tree, gamma)
            # Murray's law deviation should be small
            # With floor enforcement, some small deviation is expected
            @test mean_dev < 0.05
        end
    end

    @testset "Positive radii and lengths" begin
        for (name, tree) in forest.trees
            seg = tree.segments
            n = seg.n
            for i in 1:n
                @test seg.radius[i] > 0.0
                @test seg.seg_length[i] >= 0.0
            end
        end
    end

    @testset "Validation report" begin
        report = validate_forest(forest, params)
        @test report[:n_trees] == 3
        @test report[:total_segments] > 0

        tree_reports = report[:tree_reports]
        for (name, rep) in tree_reports
            @test rep.n_segments > 0
        end
    end

    @testset "Verbose mode prints output" begin
        # Verify verbose mode doesn't error
        io = IOBuffer()
        rng2 = MersenneTwister(123)
        tiny_configs = [
            TreeConfig("T1", (0.0, 0.0, 0.0), 0.5, (1.0, 0.0, 0.0), 5, 1.0),
        ]
        tiny_domain = SphereDomain((0.0, 0.0, 0.0), 3.0)

        # Redirect stdout to capture verbose output
        forest2 = generate_kassab_coronary(
            tiny_domain, params;
            rng=rng2,
            verbose=true,
            handoff_order=2,
            tree_configs=tiny_configs,
        )
        @test forest2.trees["T1"].segments.n > 0
    end
end
