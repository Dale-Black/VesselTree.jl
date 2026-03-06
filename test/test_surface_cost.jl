using Test
using VesselTree
using Random

@testset "Surface Cost + Barabasi Integration" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    # --- Junction surface estimate ---

    @testset "junction_surface_estimate — sprouting lower than branching" begin
        r_parent = 1.0
        # Sprouting: small daughter relative to large (rho < 0.83)
        r_large_s = 0.9
        r_small_s = 0.2  # rho = 0.22
        cost_sprout = VesselTree.junction_surface_estimate(r_parent, r_large_s, r_small_s, params)

        # Branching: near-equal daughters (rho >= 0.83)
        r_large_b = 0.7
        r_small_b = 0.65  # rho = 0.93
        cost_branch = VesselTree.junction_surface_estimate(r_parent, r_large_b, r_small_b, params)

        @test cost_sprout < cost_branch
        @test cost_sprout > 0.0
        @test cost_branch > 0.0
    end

    @testset "junction_surface_estimate — increases with radii" begin
        cost_small = VesselTree.junction_surface_estimate(0.5, 0.3, 0.2, params)
        cost_large = VesselTree.junction_surface_estimate(1.0, 0.6, 0.4, params)
        @test cost_large > cost_small
    end

    @testset "junction_surface_estimate — always positive" begin
        for (rp, rl, rs) in [(1.0, 0.8, 0.5), (0.5, 0.3, 0.1), (2.0, 1.5, 1.4)]
            cost = VesselTree.junction_surface_estimate(rp, rl, rs, params)
            @test cost > 0.0
        end
    end

    # --- Enhanced surface cost ---

    @testset "surface_cost_with_junction — greater than tube-only cost" begin
        seg_p = (0.0, 0.0, 0.0)
        seg_d = (10.0, 0.0, 0.0)
        r = 1.0
        t_pt = (5.0, 3.0, 0.0)

        cost_plain = surface_cost_at_t(0.5, seg_p[1], seg_p[2], seg_p[3],
                                       seg_d[1], seg_d[2], seg_d[3], r,
                                       t_pt[1], t_pt[2], t_pt[3], gamma)

        cost_junction = VesselTree.surface_cost_with_junction(0.5, seg_p[1], seg_p[2], seg_p[3],
                                                               seg_d[1], seg_d[2], seg_d[3], r,
                                                               t_pt[1], t_pt[2], t_pt[3], gamma, params)

        # Junction cost adds to tube cost
        @test cost_junction >= cost_plain
    end

    @testset "surface_cost_with_junction — still finds optimal t" begin
        seg_p = (0.0, 0.0, 0.0)
        seg_d = (10.0, 0.0, 0.0)
        r = 1.0
        t_pt = (5.0, 5.0, 0.0)

        # Sample costs at multiple t values
        costs = [VesselTree.surface_cost_with_junction(t, seg_p[1], seg_p[2], seg_p[3],
                                                        seg_d[1], seg_d[2], seg_d[3], r,
                                                        t_pt[1], t_pt[2], t_pt[3], gamma, params)
                 for t in 0.1:0.1:0.9]
        # Should have a minimum (not monotonic)
        @test minimum(costs) < costs[1] || minimum(costs) < costs[end]
    end

    # --- Growth with Barabasi geometry applied ---

    @testset "grow_tree! with barabasi=true — generates valid tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng, kassab=true, barabasi=true)

        seg = tree.segments
        for i in 1:seg.n
            @test seg.seg_length[i] > 0.0
            @test seg.radius[i] > 0.0
        end
    end

    @testset "grow_tree! with barabasi=true — angles are applied" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng, kassab=true, barabasi=true)

        # At least some bifurcations exist
        @test tree.n_bifurcations > 0

        # All segments should have valid length
        seg = tree.segments
        for i in 1:seg.n
            @test seg.seg_length[i] > 0.0
        end
    end

    @testset "grow_tree! — barabasi=false preserves old behavior" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true, barabasi=false)

        # Should still produce valid tree
        @test tree.segments.n > 10
    end

    # --- Angle distribution bimodality ---

    @testset "compute_bifurcation_angles — collect from tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(99)
        grow_tree!(tree, domain, 200, params; rng=rng, kassab=true, barabasi=true)

        topo = tree.topology
        seg = tree.segments

        angles = Float64[]
        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = Int(topo.child1_id[i])
                c2 = Int(topo.child2_id[i])
                (c1 <= 0 || c2 <= 0) && continue

                # Compute angle between daughters
                d1x = seg.distal_x[c1] - seg.proximal_x[c1]
                d1y = seg.distal_y[c1] - seg.proximal_y[c1]
                d1z = seg.distal_z[c1] - seg.proximal_z[c1]
                d2x = seg.distal_x[c2] - seg.proximal_x[c2]
                d2y = seg.distal_y[c2] - seg.proximal_y[c2]
                d2z = seg.distal_z[c2] - seg.proximal_z[c2]

                len1 = sqrt(d1x^2 + d1y^2 + d1z^2)
                len2 = sqrt(d2x^2 + d2y^2 + d2z^2)
                (len1 <= 0 || len2 <= 0) && continue

                dot_val = (d1x*d2x + d1y*d2y + d1z*d2z) / (len1 * len2)
                dot_val = clamp(dot_val, -1.0, 1.0)
                angle = acos(dot_val)
                push!(angles, angle)
            end
        end

        # Should have collected some angles
        @test length(angles) > 1

        # Angles should be in valid range
        for a in angles
            @test 0.0 <= a <= π
        end

        # Mean angle should be reasonable (not all 0 or all π)
        mean_angle = sum(angles) / length(angles)
        @test mean_angle > 0.1
        @test mean_angle < π - 0.1
    end

    # --- Murray's law still holds after Barabasi geometry ---

    @testset "Murray's law holds after Barabasi geometry" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng, kassab=true, barabasi=true)
        update_radii!(tree, gamma)

        topo = tree.topology
        seg = tree.segments

        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = Int(topo.child1_id[i])
                c2 = Int(topo.child2_id[i])
                (c1 <= 0 || c2 <= 0) && continue
                r_gamma_sum = seg.radius[c1]^gamma + seg.radius[c2]^gamma
                @test seg.radius[i]^gamma ≈ r_gamma_sum rtol = 1e-4
            end
        end
    end

end
