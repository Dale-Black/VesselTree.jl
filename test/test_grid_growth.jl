using Test
using VesselTree
using Random

@testset "Grid-Accelerated Growth" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    @testset "grid query vs brute force — nearest segments match" begin
        # Grow a moderate tree, then verify grid-based nearest matches brute force
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 1000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng)

        n = tree.segments.n
        @test n >= 20  # some terminals were added

        # Build grid and compare nearest for random test points
        cell_size = 2.0
        grid = build_grid(tree.segments, n, domain, cell_size)

        rng2 = MersenneTwister(123)
        distances_buf = zeros(n)
        K = 5
        matches = 0
        n_tests = 20

        for _ in 1:n_tests
            pt = sample_point(domain, rng2)
            tx, ty, tz = pt

            # Brute force: compute all distances
            compute_all_distances!(distances_buf, tree.segments, tx, ty, tz, n)
            brute_nearest = partialsortperm(@view(distances_buf[1:n]), 1:min(K, n))
            brute_top = Set(brute_nearest)

            # Grid query with generous radius
            nearby = query_nearby(grid, tx, ty, tz, 15.0)
            if length(nearby) >= K
                nearby_dists = [VesselTree.point_segment_distance(
                    tx, ty, tz,
                    tree.segments.proximal_x[si], tree.segments.proximal_y[si], tree.segments.proximal_z[si],
                    tree.segments.distal_x[si], tree.segments.distal_y[si], tree.segments.distal_z[si],
                ) for si in nearby]
                perm = partialsortperm(nearby_dists, 1:min(K, length(nearby)))
                grid_nearest = Set(nearby[perm[j]] for j in 1:min(K, length(perm)))
                # The closest segment from brute force should appear in grid results
                if brute_nearest[1] in grid_nearest
                    matches += 1
                end
            end
        end

        # Grid should find the true nearest in at least 80% of cases
        @test matches >= n_tests * 0.8
    end

    @testset "grow_tree! with grid — 50 terminals, valid tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng)

        n = tree.segments.n
        @test n >= 10
        @test tree.n_terminals >= 5

        # All radii positive
        for i in 1:n
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "grow_tree! with grid — Murray's law preserved" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(123)
        grow_tree!(tree, domain, 50, params; rng=rng)

        seg = tree.segments
        topo = tree.topology
        violations = 0
        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                if c1 > 0 && c2 > 0
                    lhs = seg.radius[i]^gamma
                    rhs = seg.radius[c1]^gamma + seg.radius[c2]^gamma
                    if !isapprox(lhs, rhs, rtol=1e-6)
                        violations += 1
                    end
                end
            end
        end
        @test violations == 0
    end

    @testset "grow_tree! with grid — connectivity check" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(456)
        grow_tree!(tree, domain, 30, params; rng=rng)

        n = tree.segments.n
        topo = tree.topology
        for i in 1:n
            if i == tree.root_segment_id
                @test topo.parent_id[i] == Int32(-1)
            else
                @test topo.parent_id[i] > 0
                @test topo.parent_id[i] <= n
            end
        end
    end

    @testset "grow_tree! with grid — 200 terminals in box" begin
        domain = BoxDomain((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0))
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(555)
        grow_tree!(tree, domain, 200, params; rng=rng)

        @test n_segments(tree) >= 100
        @test tree.n_bifurcations >= 30
    end

    @testset "grow_tree! with grid + kassab — asymmetric radii" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true)

        seg = tree.segments
        topo = tree.topology
        ratios = Float64[]
        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                if c1 > 0 && c2 > 0
                    r1 = seg.radius[c1]
                    r2 = seg.radius[c2]
                    push!(ratios, min(r1, r2) / max(r1, r2))
                end
            end
        end
        if length(ratios) > 2
            @test !all(r ≈ ratios[1] for r in ratios)
        end
    end

    @testset "_find_nearest_via_grid — falls back to brute force" begin
        # With very small search radius, grid returns too few => falls back to full scan
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 200)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(88)
        grow_tree!(tree, domain, 20, params; rng=rng)

        n = tree.segments.n
        cell_size = 2.0
        grid = build_grid(tree.segments, n, domain, cell_size)
        distances_buf = zeros(n)

        # Very small search radius — should fall back to brute force
        nearest = VesselTree._find_nearest_via_grid(
            grid, tree.segments, distances_buf,
            0.0, 0.0, 0.0, 0.001, 5, n
        )
        @test length(nearest) >= 1
        @test length(nearest) <= 5

        # Large search radius — should use grid path
        nearest2 = VesselTree._find_nearest_via_grid(
            grid, tree.segments, distances_buf,
            0.0, 0.0, 0.0, 20.0, 5, n
        )
        @test length(nearest2) >= 1
        @test length(nearest2) <= 5

        # Both should include the same closest segment
        @test nearest[1] == nearest2[1]
    end

    @testset "generate_coronary_forest with grid — basic" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        configs = [
            TreeConfig("A", (0.0, 0.0, 5.0), 0.5, (0.0, 0.0, -1.0), 10, 0.5),
            TreeConfig("B", (0.0, 0.0, -5.0), 0.5, (0.0, 0.0, 1.0), 10, 0.5),
        ]
        rng = MersenneTwister(42)
        forest = generate_coronary_forest(domain, params; tree_configs=configs, rng=rng)

        total = sum(t.segments.n for (_, t) in forest.trees)
        @test total >= 4  # at least root segments + some growth
        @test length(forest.trees) == 2
    end

end
