using Test
using VesselTree
using Random

@testset "Statistical Subdivision" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    @testset "estimate_total_segments — order 0 returns 1" begin
        @test estimate_total_segments(0, params) == 1
    end

    @testset "estimate_total_segments — order 1 returns reasonable count" begin
        # With bifurcation chains: each daughter creates 1 daughter + 1 continuation
        # CM[1,2] = 2.3, so order 1 → 1 + 2.3*(1+1) = 5.6 → ceil = 6 or 7
        est = estimate_total_segments(1, params)
        @test est >= 5
        @test est <= 10
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
        # Use proper capacity estimation for bifurcation chain approach
        est_tree = VascularTree("est", 500)
        add_segment!(est_tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        grow_tree!(est_tree, domain, 50, params; rng=MersenneTwister(42))
        cap = estimate_subdivision_capacity(est_tree, params)

        tree = VascularTree("test", max(cap, 1000))
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

    @testset "subdivide_terminals! — all junctions are bifurcations (no trifurcations)" begin
        tree = VascularTree("test", 10000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        n = tree.segments.n
        topo = tree.topology
        n_tri = 0
        n_single = 0
        for i in 1:n
            c1 = Int(topo.child1_id[i]); c2 = Int(topo.child2_id[i]); c3 = Int(topo.child3_id[i])
            nc = (c1 > 0 ? 1 : 0) + (c2 > 0 ? 1 : 0) + (c3 > 0 ? 1 : 0)
            nc >= 3 && (n_tri += 1)
            nc == 1 && (n_single += 1)
        end
        @test n_tri == 0
        @test n_single == 0
    end

    @testset "subdivide_terminals! — Murray's law exact at new junctions" begin
        tree = VascularTree("test", 10000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))
        n_cco = tree.segments.n

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        gamma = params.gamma
        max_dev = 0.0
        n_checked = 0
        topo = tree.topology
        for i in (n_cco + 1):tree.segments.n
            c1 = Int(topo.child1_id[i]); c2 = Int(topo.child2_id[i])
            if c1 > 0 && c2 > 0
                rp = tree.segments.radius[i]
                r1 = tree.segments.radius[c1]
                r2 = tree.segments.radius[c2]
                lhs = rp^gamma
                rhs = r1^gamma + r2^gamma
                dev = abs(lhs - rhs) / max(lhs, 1e-30)
                dev > max_dev && (max_dev = dev)
                n_checked += 1
            end
        end
        @test n_checked > 5  # we checked at least some bifurcations
        @test max_dev < 0.01  # < 1% deviation
    end

    @testset "subdivide_terminals! — S/E ratios from bifurcation chains" begin
        tree = VascularTree("test", 10000)
        # Order 3 terminal with enough capacity for subdivision
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        VesselTree.assign_strahler_orders!(tree, params)
        elements = group_into_elements(tree, params)
        se = compute_se_ratios(elements)

        # With bifurcation chains, elements with branches should have S/E > 1
        for (ord, ratio) in se
            @test ratio >= 1.0
        end
        # At least some elements should have S/E > 1 (continuation segments)
        @test any(ratio > 1.0 for (_, ratio) in se)
    end

    @testset "subdivide_terminals! — asymmetry built in (r_daughter < r_parent)" begin
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))
        n_cco = tree.segments.n

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        topo = tree.topology
        # All children of new bifurcations should have radius < parent
        violations = 0
        for i in (n_cco + 1):tree.segments.n
            c1 = Int(topo.child1_id[i]); c2 = Int(topo.child2_id[i])
            if c1 > 0 && c2 > 0
                rp = tree.segments.radius[i]
                r1 = tree.segments.radius[c1]
                r2 = tree.segments.radius[c2]
                (r1 > rp * 1.001 || r2 > rp * 1.001) && (violations += 1)
            end
        end
        @test violations == 0
    end

    @testset "estimate vs actual — within range" begin
        # First grow a CCO tree to estimate capacity
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        est_tree = VascularTree("est", 500)
        add_segment!(est_tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        grow_tree!(est_tree, domain, 30, params; rng=MersenneTwister(42))
        estimated = estimate_subdivision_capacity(est_tree, params)

        # Now do the real run with estimated capacity
        tree = VascularTree("test", max(estimated, 1000))
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        grow_tree!(tree, domain, 30, params; rng=MersenneTwister(42))
        subdivide_terminals!(tree, params; rng=MersenneTwister(42))

        actual = tree.segments.n
        # Estimate should be >= actual (it includes 1.5x buffer)
        @test estimated >= actual ÷ 3  # not more than 3x off
        @test estimated <= actual * 5   # not more than 5x off
    end

end
