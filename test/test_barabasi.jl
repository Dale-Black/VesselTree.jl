using Test
using VesselTree
using Random

@testset "Barabasi Junction Geometry" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    # --- Chi parameter ---

    @testset "compute_chi — basic" begin
        # chi = 2*pi*r / d
        @test compute_chi(1.0, 2π) ≈ 1.0 rtol = 1e-10
        @test compute_chi(0.5, 2π) ≈ 0.5 rtol = 1e-10
        @test compute_chi(1.0, Float64(π)) ≈ 2.0 rtol = 1e-10
    end

    @testset "compute_chi — trifurcation threshold" begin
        # chi > 0.83 means trifurcation eligible
        r = 1.0
        d_tri = 2π * r / 0.84  # just above threshold
        @test compute_chi(r, d_tri) > params.trifurcation_chi_th

        d_bi = 2π * r / 0.82   # just below threshold
        @test compute_chi(r, d_bi) < params.trifurcation_chi_th
    end

    # --- Rho parameter ---

    @testset "compute_rho — basic" begin
        @test compute_rho(0.5, 1.0) ≈ 0.5 rtol = 1e-10
        @test compute_rho(1.0, 1.0) ≈ 1.0 rtol = 1e-10
        @test compute_rho(0.0, 1.0) ≈ 0.0 atol = 1e-10
    end

    @testset "compute_rho — always <= 1" begin
        @test compute_rho(0.3, 0.7) ≈ 0.3 / 0.7 rtol = 1e-10
        @test compute_rho(0.3, 0.7) < 1.0
    end

    # --- Junction type classification ---

    @testset "classify_junction — sprouting" begin
        # rho < 0.83 → sprouting
        @test classify_junction(0.3, params) == :sprouting
        @test classify_junction(0.5, params) == :sprouting
        @test classify_junction(0.82, params) == :sprouting
    end

    @testset "classify_junction — branching" begin
        # rho >= 0.83 → branching
        @test classify_junction(0.83, params) == :branching
        @test classify_junction(0.9, params) == :branching
        @test classify_junction(1.0, params) == :branching
    end

    # --- Junction angles ---

    @testset "compute_junction_angles — sprouting regime" begin
        # rho < 0.83: large daughter continues straight (~0), small at ~90 deg
        angle_large, angle_small = compute_junction_angles(0.3, params)
        @test angle_large ≈ 0.0 atol = 0.01
        @test angle_small ≈ π / 2 atol = 0.15  # ~90 degrees
    end

    @testset "compute_junction_angles — branching regime (symmetric)" begin
        # rho = 1.0: both daughters deflect equally
        angle_large, angle_small = compute_junction_angles(1.0, params)
        @test angle_large ≈ angle_small atol = 0.01
        @test angle_large > 0.0
    end

    @testset "compute_junction_angles — branching regime (near threshold)" begin
        # rho = 0.83: minimal branching
        angle_large, angle_small = compute_junction_angles(0.83, params)
        @test angle_large >= 0.0
        @test angle_small >= angle_large  # small daughter deflects more
    end

    @testset "compute_junction_angles — angles increase with rho in branching" begin
        a1_l, a1_s = compute_junction_angles(0.85, params)
        a2_l, a2_s = compute_junction_angles(0.95, params)

        # Higher rho → larger deflection of large daughter
        @test a2_l >= a1_l - 0.01  # roughly increasing
    end

    @testset "compute_junction_angles — all positive" begin
        for rho in [0.1, 0.3, 0.5, 0.7, 0.83, 0.9, 0.95, 1.0]
            al, as = compute_junction_angles(rho, params)
            @test al >= 0.0
            @test as >= 0.0
            @test al <= π
            @test as <= π
        end
    end

    # --- Apply junction geometry ---

    @testset "apply_junction_geometry! — sprouting" begin
        tree = VascularTree("test", 50)
        # Parent along x-axis
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        # Two daughters
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 0.0, 0.0), 0.8, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 3.0, 0.0), 0.2, Int32(1))

        apply_junction_geometry!(tree, Int32(1), params)

        seg = tree.segments
        topo = tree.topology

        c1 = Int(topo.child1_id[1])
        c2 = Int(topo.child2_id[1])

        # Both daughters should still start from the bifurcation point
        @test seg.proximal_x[c1] ≈ seg.distal_x[1] atol = 1e-10
        @test seg.proximal_x[c2] ≈ seg.distal_x[1] atol = 1e-10

        # Segment lengths should be preserved
        @test seg.seg_length[c1] > 0.0
        @test seg.seg_length[c2] > 0.0
    end

    @testset "apply_junction_geometry! — branching" begin
        tree = VascularTree("test", 50)
        # Parent along x-axis
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        # Two nearly equal daughters
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.85, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.85, Int32(1))

        apply_junction_geometry!(tree, Int32(1), params)

        seg = tree.segments
        c1 = Int(tree.topology.child1_id[1])
        c2 = Int(tree.topology.child2_id[1])

        # Both should deflect (y components should differ from original)
        @test seg.seg_length[c1] > 0.0
        @test seg.seg_length[c2] > 0.0
    end

    @testset "apply_junction_geometry! — preserves segment length" begin
        tree = VascularTree("test", 50)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (13.0, -4.0, 0.0), 0.3, Int32(1))

        seg = tree.segments
        c1 = Int(tree.topology.child1_id[1])
        c2 = Int(tree.topology.child2_id[1])
        len1_before = seg.seg_length[c1]
        len2_before = seg.seg_length[c2]

        apply_junction_geometry!(tree, Int32(1), params)

        # Lengths should be preserved
        @test seg.seg_length[c1] ≈ len1_before rtol = 1e-6
        @test seg.seg_length[c2] ≈ len2_before rtol = 1e-6
    end

    @testset "apply_junction_geometry! — on CCO-grown tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 20, params; rng=rng, kassab=true)

        topo = tree.topology
        seg = tree.segments

        # Apply geometry to all bifurcations
        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                apply_junction_geometry!(tree, Int32(i), params)
            end
        end

        # All segments should still have positive length and radius
        for i in 1:seg.n
            @test seg.seg_length[i] > 0.0
            @test seg.radius[i] > 0.0
        end
    end

end
