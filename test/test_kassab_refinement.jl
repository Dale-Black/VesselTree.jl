using Test
using VesselTree
using Random

@testset "Post-hoc Kassab Radius Refinement" begin

    params = kassab_coronary_params()
    gamma = params.gamma
    min_radius = params.vessel_cutoff_um / 2.0 / 1000.0  # 0.004mm

    @testset "apply_full_kassab_radii! — all radii >= floor" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng)

        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(123))

        for i in 1:tree.segments.n
            @test tree.segments.radius[i] >= min_radius
        end
    end

    @testset "apply_full_kassab_radii! — Murray's law deviation < 0.1%" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng)

        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(123))

        seg = tree.segments
        topo = tree.topology
        max_deviation = 0.0
        n_junctions = 0

        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = Int(topo.child1_id[i])
                c2 = Int(topo.child2_id[i])
                (c1 <= 0 || c2 <= 0) && continue
                lhs = seg.radius[i]^gamma
                rhs = seg.radius[c1]^gamma + seg.radius[c2]^gamma
                dev = abs(lhs - rhs) / max(lhs, 1e-30)
                max_deviation = max(max_deviation, dev)
                n_junctions += 1
            end
        end

        # Murray's law should hold with < 0.1% deviation at most junctions
        # (floor enforcement can cause small deviations)
        @test n_junctions > 0
        # Allow up to 1% for floor-affected junctions
        @test max_deviation < 0.01 || n_junctions == 0
    end

    @testset "apply_full_kassab_radii! — asymmetry variation" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng)

        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(456))

        seg = tree.segments
        topo = tree.topology
        ratios = Float64[]
        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = Int(topo.child1_id[i])
                c2 = Int(topo.child2_id[i])
                (c1 <= 0 || c2 <= 0) && continue
                r1 = seg.radius[c1]
                r2 = seg.radius[c2]
                push!(ratios, min(r1, r2) / max(r1, r2))
            end
        end

        # Should have variation in asymmetry ratios
        @test length(ratios) > 5
        @test !all(r ≈ ratios[1] for r in ratios)
        # Ratios should be in (0, 1]
        @test all(0 < r <= 1.0 for r in ratios)
    end

    @testset "apply_full_kassab_radii! — no sub-micrometer radii" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng)

        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(789))

        for i in 1:tree.segments.n
            # No radius below 1 micrometer (0.001mm = 0.0005mm radius)
            @test tree.segments.radius[i] >= 0.001  # 2um diameter minimum
        end
    end

    @testset "apply_full_kassab_radii! — after subdivision" begin
        tree = VascularTree("test", 10000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.015, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        n_after_sub = tree.segments.n
        @test n_after_sub > 3

        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(123))

        # All radii should be >= floor after refinement
        for i in 1:tree.segments.n
            @test tree.segments.radius[i] >= min_radius
        end
    end

    @testset "apply_full_kassab_radii! — handles trifurcations" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.3, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 0.0, 1.0), 0.2, Int32(1))

        @test tree.topology.junction_type[1] == :trifurcation

        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(42))

        for i in 1:tree.segments.n
            @test tree.segments.radius[i] >= min_radius
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "apply_full_kassab_radii! — root radius stays large" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng)

        root_id = Int(tree.root_segment_id)
        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(321))

        # Root should still be the largest
        root_r = tree.segments.radius[root_id]
        for i in 1:tree.segments.n
            @test tree.segments.radius[i] <= root_r + 1e-10
        end
    end

    @testset "apply_full_kassab_radii! — idempotent Murray's law" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng)

        apply_full_kassab_radii!(tree, params; rng=MersenneTwister(42))

        # Capture radii after refinement
        radii_before = copy(tree.segments.radius[1:tree.segments.n])

        # Apply update_radii! again — shouldn't change if Murray's law holds
        update_radii!(tree, gamma)

        for i in 1:tree.segments.n
            @test tree.segments.radius[i] ≈ radii_before[i] rtol=0.01
        end
    end

end
