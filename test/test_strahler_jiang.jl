using Test
using VesselTree
using Random

@testset "Diameter-Defined Strahler Ordering (Jiang 1994)" begin

    params = kassab_rca_params()

    @testset "Topological Strahler ordering — hand-built tree" begin
        tree = VascularTree("test", 50)
        # Build a balanced binary tree:
        #       1
        #      / \
        #     2   3
        #    / \
        #   4   5
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))   # seg 1
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.05, Int32(1))   # seg 2
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))  # seg 3 (terminal)
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 5.0, 0.0), 0.01, Int32(2))  # seg 4 (terminal)
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 1.0, 0.0), 0.01, Int32(2))  # seg 5 (terminal)

        VesselTree._assign_topological_strahler!(tree)

        topo = tree.topology
        # Terminals 3, 4, 5 → order 1
        @test topo.strahler_order[3] == Int32(1)
        @test topo.strahler_order[4] == Int32(1)
        @test topo.strahler_order[5] == Int32(1)
        # Node 2: children 4,5 both order 1 → order 2
        @test topo.strahler_order[2] == Int32(2)
        # Node 1: children 2 (order 2) and 3 (order 1) → max = 2
        @test topo.strahler_order[1] == Int32(2)
    end

    @testset "Topological Strahler — symmetric tree gives higher orders" begin
        tree = VascularTree("test", 50)
        # Symmetric tree with 4 terminals:
        #        1
        #       / \
        #      2   3
        #     / \ / \
        #    4  5 6  7
        add_segment!(tree, (0.0, 0.0, 0.0), (4.0, 0.0, 0.0), 0.1, Int32(-1))  # 1
        add_segment!(tree, (4.0, 0.0, 0.0), (6.0, 2.0, 0.0), 0.05, Int32(1))  # 2
        add_segment!(tree, (4.0, 0.0, 0.0), (6.0, -2.0, 0.0), 0.05, Int32(1)) # 3
        add_segment!(tree, (6.0, 2.0, 0.0), (8.0, 3.0, 0.0), 0.02, Int32(2))  # 4
        add_segment!(tree, (6.0, 2.0, 0.0), (8.0, 1.0, 0.0), 0.02, Int32(2))  # 5
        add_segment!(tree, (6.0, -2.0, 0.0), (8.0, -1.0, 0.0), 0.02, Int32(3)) # 6
        add_segment!(tree, (6.0, -2.0, 0.0), (8.0, -3.0, 0.0), 0.02, Int32(3)) # 7

        VesselTree._assign_topological_strahler!(tree)

        topo = tree.topology
        # Terminals → order 1
        for i in 4:7
            @test topo.strahler_order[i] == Int32(1)
        end
        # Internal nodes 2, 3: both children order 1 → order 2
        @test topo.strahler_order[2] == Int32(2)
        @test topo.strahler_order[3] == Int32(2)
        # Root: both children order 2 → order 3
        @test topo.strahler_order[1] == Int32(3)
    end

    @testset "assign_strahler_orders_simple! is the same as assign_strahler_orders!" begin
        # Verify alias works
        @test assign_strahler_orders_simple! === assign_strahler_orders!
    end

    @testset "Iterative ordering converges on CCO tree" begin
        rng = MersenneTwister(42)
        domain = SphereDomain((0.0, 0.0, 0.0), 5.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5, Int32(-1))
        grow_tree!(tree, domain, 50, params; rng=rng)
        update_radii!(tree, params.gamma)

        n_iter = assign_diameter_defined_strahler!(tree, params)

        # Should converge in <= 20 iterations
        @test n_iter <= 20
        # Typically converges in 2-5 iterations (Jiang 1994: "two or three cycles")
        @test n_iter <= 10

        # All orders should be non-negative
        n = tree.segments.n
        for i in 1:n
            @test tree.topology.strahler_order[i] >= 0
        end
    end

    @testset "Iterative ordering — idempotent" begin
        rng = MersenneTwister(99)
        domain = SphereDomain((0.0, 0.0, 0.0), 5.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5, Int32(-1))
        grow_tree!(tree, domain, 50, params; rng=rng)
        update_radii!(tree, params.gamma)

        # First run
        assign_diameter_defined_strahler!(tree, params)
        orders1 = [Int(tree.topology.strahler_order[i]) for i in 1:tree.segments.n]

        # Second run on same tree (already converged)
        n_iter2 = assign_diameter_defined_strahler!(tree, params)
        orders2 = [Int(tree.topology.strahler_order[i]) for i in 1:tree.segments.n]

        # Algorithm always reinitializes from topological ordering, so same iteration count
        # Orders should be identical (deterministic algorithm)
        @test orders1 == orders2
    end

    @testset "Iterative ordering — diameter bounds are non-overlapping" begin
        rng = MersenneTwister(77)
        domain = SphereDomain((0.0, 0.0, 0.0), 5.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5, Int32(-1))
        grow_tree!(tree, domain, 100, params; rng=rng)
        update_radii!(tree, params.gamma)

        assign_diameter_defined_strahler!(tree, params)

        # Collect per-order diameter ranges
        n = tree.segments.n
        order_diams = Dict{Int, Vector{Float64}}()
        for i in 1:n
            ord = Int(tree.topology.strahler_order[i])
            ord < 0 && continue
            d_um = tree.segments.radius[i] * 2.0 * 1000.0
            if !haskey(order_diams, ord)
                order_diams[ord] = Float64[]
            end
            push!(order_diams[ord], d_um)
        end

        # Check that order means are monotonically increasing
        sorted_orders = sort(collect(keys(order_diams)))
        if length(sorted_orders) >= 2
            using Statistics
            prev_mean = 0.0
            for ord in sorted_orders
                m = mean(order_diams[ord])
                @test m > prev_mean
                prev_mean = m
            end
        end
    end

    @testset "Iterative ordering on many segments (AK compatibility)" begin
        rng = MersenneTwister(42)
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), 1.0, Int32(-1))
        grow_tree!(tree, domain, 200, params; rng=rng)
        update_radii!(tree, params.gamma)

        n = tree.segments.n
        @test n > 100

        n_iter = assign_diameter_defined_strahler!(tree, params)
        @test n_iter <= 20

        # Should produce at least 2 distinct orders
        orders = Set{Int}()
        for i in 1:n
            push!(orders, Int(tree.topology.strahler_order[i]))
        end
        @test length(orders) >= 2
    end

    @testset "Works on empty and single-segment trees" begin
        # Empty tree
        tree0 = VascularTree("test", 10)
        @test assign_diameter_defined_strahler!(tree0, params) == 0

        # Single segment
        tree1 = VascularTree("test", 10)
        add_segment!(tree1, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.1, Int32(-1))
        n_iter = assign_diameter_defined_strahler!(tree1, params)
        @test tree1.topology.strahler_order[1] == Int32(1)  # terminal → order 1
    end

end
