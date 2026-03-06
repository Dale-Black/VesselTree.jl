using Test
using VesselTree
using Random

@testset "Multi-Tree Forest Growth" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    # --- CoronaryForest struct ---

    @testset "CoronaryForest — construction" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        configs = coronary_tree_configs()
        tmap = initialize_territories(domain, configs)
        trees = Dict{String, VascularTree}()
        for cfg in configs
            trees[cfg.name] = VascularTree(cfg.name, 100)
        end
        forest = CoronaryForest(trees, tmap, params)
        @test forest isa CoronaryForest
        @test length(forest.trees) == 3
        @test haskey(forest.trees, "LAD")
    end

    # --- generate_coronary_forest ---

    @testset "generate_coronary_forest — small forest" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        configs = [
            TreeConfig("A", (-5.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 10, 0.5),
            TreeConfig("B", (5.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 10, 0.5),
        ]

        rng = MersenneTwister(42)
        forest = generate_coronary_forest(domain, params; tree_configs=configs, rng=rng)

        @test forest isa CoronaryForest
        @test haskey(forest.trees, "A")
        @test haskey(forest.trees, "B")

        # Both trees should have some segments
        @test forest.trees["A"].segments.n > 1
        @test forest.trees["B"].segments.n > 1
    end

    @testset "generate_coronary_forest — all segments valid" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        configs = [
            TreeConfig("T1", (-4.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 15, 0.5),
            TreeConfig("T2", (4.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 15, 0.5),
        ]

        rng = MersenneTwister(99)
        forest = generate_coronary_forest(domain, params; tree_configs=configs, rng=rng)

        for (name, tree) in forest.trees
            for i in 1:tree.segments.n
                @test tree.segments.seg_length[i] > 0.0
                @test tree.segments.radius[i] > 0.0
            end
        end
    end

    @testset "generate_coronary_forest — round-robin balances trees" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        configs = [
            TreeConfig("X", (-4.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 20, 0.5),
            TreeConfig("Y", (4.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 20, 0.5),
        ]

        rng = MersenneTwister(42)
        forest = generate_coronary_forest(domain, params; tree_configs=configs, rng=rng)

        n_x = forest.trees["X"].segments.n
        n_y = forest.trees["Y"].segments.n
        # Trees should be roughly balanced (within 2x)
        @test n_x > 0
        @test n_y > 0
        ratio = max(n_x, n_y) / max(1, min(n_x, n_y))
        @test ratio < 5.0
    end

    # --- validate_forest ---

    @testset "validate_forest — basic" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        configs = [
            TreeConfig("A", (-4.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 10, 0.5),
            TreeConfig("B", (4.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 10, 0.5),
        ]

        rng = MersenneTwister(42)
        forest = generate_coronary_forest(domain, params; tree_configs=configs, rng=rng)

        report = validate_forest(forest, params)
        @test report isa Dict
        @test haskey(report, :total_segments)
        @test report[:total_segments] > 0
    end

    # --- Default coronary configs ---

    @testset "generate with coronary_tree_configs — runs without error" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        configs = [
            TreeConfig("LAD", (-4.0, 2.0, 0.0), 1.0, (1.0, -1.0, 0.0), 10, 0.40),
            TreeConfig("LCX", (-4.0, -2.0, 0.0), 0.8, (1.0, -1.0, 0.0), 8, 0.25),
            TreeConfig("RCA", (4.0, 0.0, 0.0), 0.9, (-1.0, -1.0, 0.0), 8, 0.35),
        ]

        rng = MersenneTwister(42)
        forest = generate_coronary_forest(domain, params; tree_configs=configs, rng=rng)

        @test length(forest.trees) == 3
        for (name, tree) in forest.trees
            @test tree.segments.n >= 1
        end
    end

end
