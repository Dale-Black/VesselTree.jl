using Test
using VesselTree
using Random

@testset "Trifurcation Detection + Handling" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    # --- compute_trifurcation_angles ---

    @testset "compute_trifurcation_angles — symmetric" begin
        # Three equal daughters → equal angles, spaced 120 deg apart from parent dir
        a1, a2, a3 = compute_trifurcation_angles(1.0, 1.0, 1.0, params)
        # All three should be equal for symmetric case
        @test a1 ≈ a2 atol = 0.01
        @test a2 ≈ a3 atol = 0.01
        # All should be positive
        @test a1 > 0.0
        @test a2 > 0.0
        @test a3 > 0.0
        # Each should be less than π
        @test a1 < π
    end

    @testset "compute_trifurcation_angles — asymmetric" begin
        # One large, two small → large continues near straight, smalls deflect more
        a1, a2, a3 = compute_trifurcation_angles(0.8, 0.3, 0.2, params)
        # Largest daughter (0.8) should have smallest deflection
        @test a1 <= a2
        @test a1 <= a3
        # All positive
        @test a1 >= 0.0
        @test a2 > 0.0
        @test a3 > 0.0
    end

    @testset "compute_trifurcation_angles — all positive bounded" begin
        for (r1, r2, r3) in [(1.0, 1.0, 1.0), (0.8, 0.5, 0.3), (0.9, 0.9, 0.1)]
            a1, a2, a3 = compute_trifurcation_angles(r1, r2, r3, params)
            @test a1 >= 0.0
            @test a2 >= 0.0
            @test a3 >= 0.0
            @test a1 <= π
            @test a2 <= π
            @test a3 <= π
        end
    end

    # --- check_trifurcation_merge ---

    @testset "check_trifurcation_merge — finds nearby bifurcation" begin
        tree = VascularTree("test", 100)
        # Root along x-axis
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        # Two daughters (makes segment 1 a bifurcation)
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))

        @test tree.topology.junction_type[1] == :bifurcation

        # New bifurcation point very close to existing distal (bifurcation) point of segment 1
        # With large parent radius → chi = 2*pi*r / d should be > 0.83
        result = check_trifurcation_merge(tree, (10.01, 0.0, 0.0), params)
        @test result !== nothing
        @test result == Int32(1)
    end

    @testset "check_trifurcation_merge — returns nothing for distant point" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))

        # Very far away — chi will be tiny
        result = check_trifurcation_merge(tree, (100.0, 100.0, 100.0), params)
        @test result === nothing
    end

    @testset "check_trifurcation_merge — skips terminals and trifurcations" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        # Only 1 child → still :bifurcation but with child2 = -1
        # Actually, 1 child makes it :bifurcation already (from add_segment! logic)
        # Let's make 3 children so it's :trifurcation
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 0.0, 1.0), 0.5, Int32(1))

        @test tree.topology.junction_type[1] == :trifurcation

        # Even though close, should skip because already trifurcation
        result = check_trifurcation_merge(tree, (10.01, 0.0, 0.0), params)
        @test result === nothing
    end

    # --- merge_to_trifurcation! ---

    @testset "merge_to_trifurcation! — basic" begin
        tree = VascularTree("test", 100)
        # Root along x-axis
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        # Two daughters (bifurcation)
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))

        @test tree.topology.junction_type[1] == :bifurcation
        @test tree.n_bifurcations == 1
        @test tree.n_trifurcations == 0

        # Merge a new terminal onto the existing bifurcation
        new_terminal = (15.0, 0.0, 2.0)
        terminal_radius = 0.3
        new_id = merge_to_trifurcation!(tree, Int32(1), new_terminal, terminal_radius, params)

        @test new_id > 0
        @test tree.topology.junction_type[1] == :trifurcation
        @test tree.n_trifurcations == 1
        @test tree.n_bifurcations == 0

        # Third child should be set
        @test tree.topology.child3_id[1] == new_id
        @test tree.topology.parent_id[new_id] == Int32(1)
        @test tree.topology.is_terminal[new_id] == true

        # New segment starts from parent's distal point
        @test tree.segments.proximal_x[new_id] ≈ tree.segments.distal_x[1] atol = 1e-10
        @test tree.segments.proximal_y[new_id] ≈ tree.segments.distal_y[1] atol = 1e-10
        @test tree.segments.proximal_z[new_id] ≈ tree.segments.distal_z[1] atol = 1e-10

        # Segment endpoint should be the terminal point
        @test tree.segments.distal_x[new_id] ≈ 15.0 atol = 1e-10
        @test tree.segments.distal_y[new_id] ≈ 0.0 atol = 1e-10
        @test tree.segments.distal_z[new_id] ≈ 2.0 atol = 1e-10

        @test tree.segments.radius[new_id] ≈ terminal_radius
    end

    @testset "merge_to_trifurcation! — Murray's law holds" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.4, Int32(1))

        merge_to_trifurcation!(tree, Int32(1), (15.0, 0.0, 2.0), 0.3, params)

        # Now update radii
        update_radii!(tree, gamma)

        # Murray's law at trifurcation: r_parent^gamma = r1^gamma + r2^gamma + r3^gamma
        seg = tree.segments
        topo = tree.topology
        c1 = Int(topo.child1_id[1])
        c2 = Int(topo.child2_id[1])
        c3 = Int(topo.child3_id[1])

        r_gamma_sum = seg.radius[c1]^gamma + seg.radius[c2]^gamma + seg.radius[c3]^gamma
        @test seg.radius[1]^gamma ≈ r_gamma_sum rtol = 1e-6
    end

    # --- apply_trifurcation_geometry! ---

    @testset "apply_trifurcation_geometry! — basic" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 0.0, 1.0), 0.3, Int32(1))

        seg = tree.segments
        c1 = Int(tree.topology.child1_id[1])
        c2 = Int(tree.topology.child2_id[1])
        c3 = Int(tree.topology.child3_id[1])

        len1 = seg.seg_length[c1]
        len2 = seg.seg_length[c2]
        len3 = seg.seg_length[c3]

        apply_trifurcation_geometry!(tree, Int32(1), params)

        # Lengths should be preserved
        @test seg.seg_length[c1] ≈ len1 rtol = 1e-6
        @test seg.seg_length[c2] ≈ len2 rtol = 1e-6
        @test seg.seg_length[c3] ≈ len3 rtol = 1e-6

        # All daughters should start from the bifurcation point
        @test seg.proximal_x[c1] ≈ seg.distal_x[1] atol = 1e-10
        @test seg.proximal_x[c2] ≈ seg.distal_x[1] atol = 1e-10
        @test seg.proximal_x[c3] ≈ seg.distal_x[1] atol = 1e-10
    end

    @testset "apply_trifurcation_geometry! — preserves segment length" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (14.0, 2.0, 0.0), 0.6, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (14.0, -2.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (13.0, 0.0, 3.0), 0.3, Int32(1))

        seg = tree.segments
        c1 = Int(tree.topology.child1_id[1])
        c2 = Int(tree.topology.child2_id[1])
        c3 = Int(tree.topology.child3_id[1])

        len1 = seg.seg_length[c1]
        len2 = seg.seg_length[c2]
        len3 = seg.seg_length[c3]

        apply_trifurcation_geometry!(tree, Int32(1), params)

        @test seg.seg_length[c1] ≈ len1 rtol = 1e-6
        @test seg.seg_length[c2] ≈ len2 rtol = 1e-6
        @test seg.seg_length[c3] ≈ len3 rtol = 1e-6
    end

    # --- Integration: grow_tree! with trifurcation ---

    @testset "grow_tree! — trifurcations created in large tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 100, params; rng=rng, kassab=true, trifurcation=true)

        # Should have created some segments (trifurcations add fewer segments per terminal)
        @test tree.segments.n > 10

        # All segments should have positive length and radius
        for i in 1:tree.segments.n
            @test tree.segments.seg_length[i] > 0.0
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "grow_tree! — trifurcation percentage in range" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(123)
        grow_tree!(tree, domain, 300, params; rng=rng, kassab=true, trifurcation=true)

        total_junctions = tree.n_bifurcations + tree.n_trifurcations
        if total_junctions > 0
            tri_pct = tree.n_trifurcations / total_junctions * 100.0
            # Barabasi: ~8% trifurcation rate in blood vessels; we accept 0-25%
            @test tri_pct >= 0.0
            @test tri_pct <= 25.0
        end
    end

    @testset "grow_tree! — Murray's law at trifurcations" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(77)
        grow_tree!(tree, domain, 80, params; rng=rng, kassab=true, trifurcation=true)

        # Update radii to ensure consistency
        update_radii!(tree, gamma)

        topo = tree.topology
        seg = tree.segments

        for i in 1:seg.n
            if topo.junction_type[i] == :trifurcation
                c1 = Int(topo.child1_id[i])
                c2 = Int(topo.child2_id[i])
                c3 = Int(topo.child3_id[i])
                if c1 > 0 && c2 > 0 && c3 > 0
                    r_gamma_sum = seg.radius[c1]^gamma + seg.radius[c2]^gamma + seg.radius[c3]^gamma
                    @test seg.radius[i]^gamma ≈ r_gamma_sum rtol = 1e-4
                end
            end
        end
    end

    @testset "grow_tree! — trifurcation=false preserves old behavior" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true, trifurcation=false)

        # No trifurcations should be created
        @test tree.n_trifurcations == 0
    end

end
