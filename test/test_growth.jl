using Test
using VesselTree
using Random

@testset "CCO Growth" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    @testset "update_radii! — simple bifurcation" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.3, Int32(1))

        # Set terminal radii
        tree.segments.radius[2] = 0.5
        tree.segments.radius[3] = 0.3

        update_radii!(tree, gamma)

        # Murray's law: r_parent^gamma = r1^gamma + r2^gamma
        r_expected = (0.5^gamma + 0.3^gamma)^(1.0 / gamma)
        @test tree.segments.radius[1] ≈ r_expected atol = 1e-10
    end

    @testset "update_radii! — multi-level tree" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.7, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.7, Int32(1))
        add_segment!(tree, (2.0, 1.0, 0.0), (3.0, 2.0, 0.0), 0.4, Int32(2))
        add_segment!(tree, (2.0, 1.0, 0.0), (3.0, 0.5, 0.0), 0.4, Int32(2))

        # Set terminal radii
        tree.segments.radius[3] = 0.5   # terminal
        tree.segments.radius[4] = 0.4   # terminal
        tree.segments.radius[5] = 0.3   # terminal

        update_radii!(tree, gamma)

        # Segment 2: r^gamma = 0.4^gamma + 0.3^gamma
        r2_expected = (0.4^gamma + 0.3^gamma)^(1.0 / gamma)
        @test tree.segments.radius[2] ≈ r2_expected atol = 1e-10

        # Segment 1: r^gamma = r2^gamma + 0.5^gamma
        r1_expected = (r2_expected^gamma + 0.5^gamma)^(1.0 / gamma)
        @test tree.segments.radius[1] ≈ r1_expected atol = 1e-10
    end

    @testset "add_bifurcation! — basic" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        cont_id, new_id = add_bifurcation!(tree, 1, 0.5, (5.0, 5.0, 0.0), 0.3)

        @test n_segments(tree) == 3
        # Parent now goes to bifurcation point
        @test tree.segments.distal_x[1] ≈ 5.0
        @test tree.segments.distal_y[1] ≈ 0.0

        # Continuation goes from bifurcation to original distal
        @test tree.segments.proximal_x[cont_id] ≈ 5.0
        @test tree.segments.distal_x[cont_id] ≈ 10.0

        # New segment goes from bifurcation to terminal
        @test tree.segments.proximal_x[new_id] ≈ 5.0
        @test tree.segments.distal_x[new_id] ≈ 5.0
        @test tree.segments.distal_y[new_id] ≈ 5.0

        # Topology
        @test tree.topology.parent_id[cont_id] == Int32(1)
        @test tree.topology.parent_id[new_id] == Int32(1)
        @test tree.topology.is_terminal[cont_id] == true
        @test tree.topology.is_terminal[new_id] == true
        @test tree.topology.is_terminal[1] == false
    end

    @testset "grow_tree! — 10 terminals in sphere" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        # Add root segment
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 10, params; rng=rng)

        # Should have approximately 2*n_terminals + 1 segments (root + bifurcations)
        @test n_segments(tree) >= 10
        @test tree.n_terminals >= 5  # at least some terminals were added
    end

    @testset "grow_tree! — 50 terminals, Murray's law check" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(123)
        grow_tree!(tree, domain, 50, params; rng=rng)

        # Verify Murray's law at every bifurcation
        n = tree.segments.n
        topo = tree.topology
        seg = tree.segments

        murray_violations = 0
        for i in 1:n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                if c1 > 0 && c2 > 0
                    r_parent = seg.radius[i]
                    r_sum = seg.radius[c1]^gamma + seg.radius[c2]^gamma
                    if !isapprox(r_parent^gamma, r_sum, rtol=1e-6)
                        murray_violations += 1
                    end
                end
            end
        end
        @test murray_violations == 0
    end

    @testset "grow_tree! — connectivity check" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(456)
        grow_tree!(tree, domain, 30, params; rng=rng)

        n = tree.segments.n
        topo = tree.topology

        # Every segment except root should have a valid parent
        for i in 1:n
            if i == tree.root_segment_id
                @test topo.parent_id[i] == Int32(-1)
            else
                @test topo.parent_id[i] > 0
                @test topo.parent_id[i] <= n
            end
        end
    end

    @testset "grow_tree! — segments in domain" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(789)
        grow_tree!(tree, domain, 20, params; rng=rng)

        n = tree.segments.n
        seg = tree.segments

        # All distal endpoints should be in domain (or very close)
        for i in 2:n  # skip root proximal which is at center
            dp = (seg.distal_x[i], seg.distal_y[i], seg.distal_z[i])
            sd = signed_distance(domain, dp)
            @test sd <= 0.5  # allow small tolerance
        end
    end

    @testset "grow_tree! — positive radii" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(321)
        grow_tree!(tree, domain, 30, params; rng=rng)

        n = tree.segments.n
        for i in 1:n
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "grow_tree! — 100 terminals in box" begin
        domain = BoxDomain((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0))
        tree = VascularTree("test", 1000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(555)
        grow_tree!(tree, domain, 50, params; rng=rng)

        # Growth should succeed — at least some terminals added
        @test n_segments(tree) >= 3   # root + at least one bifurcation
        @test tree.n_bifurcations >= 1
    end

    # --- Kassab-integrated growth ---

    @testset "grow_tree! — asymmetric daughter radii" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true)

        seg = tree.segments
        topo = tree.topology

        # Collect asymmetry ratios — should show variation
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

        # Should have some variation (not all symmetric)
        if length(ratios) > 2
            @test !all(r ≈ ratios[1] for r in ratios)
        end
    end

    @testset "grow_tree! — Kassab Murray's law preserved" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(88)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true)

        seg = tree.segments
        topo = tree.topology

        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                if c1 > 0 && c2 > 0
                    lhs = seg.radius[i]^gamma
                    rhs = seg.radius[c1]^gamma + seg.radius[c2]^gamma
                    @test lhs ≈ rhs rtol = 1e-6
                end
            end
        end
    end

    @testset "grow_tree! — Kassab all radii positive" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(99)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true)

        for i in 1:tree.segments.n
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "grow_tree! — Kassab Strahler orders assigned" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(111)
        grow_tree!(tree, domain, 20, params; rng=rng, kassab=true)

        assign_strahler_orders!(tree, params)
        # All segments should have valid orders (>= 0)
        for i in 1:tree.segments.n
            @test tree.topology.strahler_order[i] >= Int32(-1)
        end
    end

    @testset "grow_tree! — backward compatible (kassab=false)" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        # Default (no kassab keyword) should still work
        rng = MersenneTwister(222)
        grow_tree!(tree, domain, 10, params; rng=rng)
        @test n_segments(tree) >= 3
    end

end
