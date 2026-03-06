using Test
using VesselTree
using Random

@testset "Performance: generate_kassab_coronary" begin
    params = kassab_coronary_params()
    rng = MersenneTwister(99)

    domain = SphereDomain((0.0, 0.0, 0.0), 5.0)
    test_configs = [
        TreeConfig("PERF", (0.0, 0.0, 0.0), 1.5, (0.0, -1.0, 0.0), 100, 1.0),
    ]

    # Measure wall-clock time
    t0 = time()
    forest = generate_kassab_coronary(
        domain, params;
        rng=rng, verbose=false, handoff_order=5,
        tree_configs=test_configs,
    )
    elapsed = time() - t0
    tree = forest.trees["PERF"]
    n_segments = tree.segments.n

    @testset "Segment count scales with terminals" begin
        # 100 CCO terminals with handoff_order=5 should produce thousands of segments
        @test n_segments > 1000
    end

    @testset "All radii positive" begin
        @test all(tree.segments.radius[i] > 0 for i in 1:n_segments)
    end

    @testset "Timing within budget" begin
        # 100-terminal tree should complete in < 30s
        @test elapsed < 30.0
    end

    @testset "Murray's law O(n) update_radii!" begin
        # update_radii! should handle large trees in < 1s
        t1 = time()
        VesselTree.update_radii!(tree, params.gamma)
        t_murray = time() - t1
        @test t_murray < 1.0
    end
end
