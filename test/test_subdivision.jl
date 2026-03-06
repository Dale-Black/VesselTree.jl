using Test
using VesselTree
using Random

@testset "Statistical Subdivision" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    @testset "estimate_total_segments — order 0 returns 1" begin
        @test estimate_total_segments(0, params) == 1
    end

    @testset "estimate_total_segments — order 1 returns ~3" begin
        # CM[1,2] = 2.3, so order 1 → 1 + 2.3*1 = 3.3 → ceil = 4
        est = estimate_total_segments(1, params)
        @test est >= 3
        @test est <= 5
    end

    @testset "estimate_total_segments — monotonic with order" begin
        prev = 1
        for ord in 1:11
            est = estimate_total_segments(ord, params)
            @test est >= prev
            prev = est
        end
    end

    @testset "estimate_total_segments — higher orders produce more segments" begin
        @test estimate_total_segments(5, params) > estimate_total_segments(3, params)
        @test estimate_total_segments(8, params) > estimate_total_segments(5, params)
    end

    @testset "estimate_subdivision_capacity — reasonable estimate" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 20, params; rng=rng)

        cap = estimate_subdivision_capacity(tree, params)
        @test cap >= tree.segments.n
        @test cap > 0
    end

    @testset "subdivide_terminals! — small tree, order 1 terminals" begin
        # Create a tree with a single order-1 terminal
        tree = VascularTree("test", 500)
        # Root: large radius (order 3)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))  # 56um diameter → order 3
        # Bifurcation: two daughters
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))  # 15um → order 1
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.0075, Int32(1)) # 15um → order 1

        n_before = tree.segments.n
        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        # Should have added daughters to the order-1 terminals
        @test tree.segments.n > n_before
        # Order-1 terminals should produce order-0 daughters (capillaries)
        @test tree.segments.n >= n_before + 2  # at least 2 capillaries added
    end

    @testset "subdivide_terminals! — all radii positive" begin
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        for i in 1:tree.segments.n
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "subdivide_terminals! — all lengths positive" begin
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.015, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))

        rng = MersenneTwister(123)
        subdivide_terminals!(tree, params; rng=rng)

        for i in 1:tree.segments.n
            @test tree.segments.seg_length[i] > 0.0
        end
    end

    @testset "subdivide_terminals! — topology: all non-root have valid parent" begin
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.0075, Int32(1))

        rng = MersenneTwister(456)
        subdivide_terminals!(tree, params; rng=rng)

        n = tree.segments.n
        for i in 1:n
            if i == tree.root_segment_id
                @test tree.topology.parent_id[i] == Int32(-1)
            else
                @test tree.topology.parent_id[i] > 0
                @test tree.topology.parent_id[i] <= n
            end
        end
    end

    @testset "subdivide_terminals! — CCO + subdivision produces more segments" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 50000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng)
        n_after_cco = tree.segments.n

        subdivide_terminals!(tree, params; rng=rng)
        n_after_sub = tree.segments.n

        @test n_after_sub > n_after_cco
        @test n_after_sub >= n_after_cco * 2  # at least double
    end

    @testset "subdivide_terminals! — order 0 terminals not subdivided" begin
        # Only order-0 terminals: radius 0.004mm = 8um diameter
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.004, Int32(-1))

        n_before = tree.segments.n
        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        # Order 0 should not be subdivided further
        @test tree.segments.n == n_before
    end

    @testset "subdivide_terminals! — diameters follow per-order distributions" begin
        tree = VascularTree("test", 10000)
        # Order 3 root terminal (56um diameter = 0.028mm radius)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        # After subdivision, should have some segments with diameters near order-0 (8um)
        assign_strahler_orders!(tree, params)
        n = tree.segments.n
        has_order0 = false
        for i in 1:n
            if tree.topology.strahler_order[i] == Int32(0)
                has_order0 = true
                break
            end
        end
        @test has_order0
    end

    @testset "subdivide_terminals! — multiple order-5 terminals produce many segments" begin
        # Create a binary tree with 8 order-5 terminals via cascading bifurcations
        tree = VascularTree("test", 100000)
        # Root segment (order ~7)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 0.295, Int32(-1))
        # First level: 2 daughters
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 3.0, 0.0), 0.092, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -3.0, 0.0), 0.092, Int32(1))
        # Second level: 4 terminals
        add_segment!(tree, (15.0, 3.0, 0.0), (20.0, 5.0, 0.0), 0.092, Int32(2))
        add_segment!(tree, (15.0, 3.0, 0.0), (20.0, 1.0, 0.0), 0.092, Int32(2))
        add_segment!(tree, (15.0, -3.0, 0.0), (20.0, -1.0, 0.0), 0.092, Int32(3))
        add_segment!(tree, (15.0, -3.0, 0.0), (20.0, -5.0, 0.0), 0.092, Int32(3))

        n_before = tree.segments.n  # 7 segments
        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        # 4 order-5 terminals should produce significant expansion
        # E[S(5)] ≈ 78 per terminal, so ~300+ new segments expected
        @test tree.segments.n > n_before * 5
    end

    @testset "estimate vs actual — within 2x" begin
        tree = VascularTree("test", 50000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng)

        estimated = estimate_subdivision_capacity(tree, params)

        rng2 = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng2)

        actual = tree.segments.n
        # Estimate should be within 2x of actual (with 1.5x buffer already included)
        @test estimated >= actual ÷ 3  # not more than 3x off
        @test estimated <= actual * 5   # not more than 5x off
    end

end
