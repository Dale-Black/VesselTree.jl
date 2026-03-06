using Test
using VesselTree
using Random

@testset "Hemodynamics" begin

    params = kassab_coronary_params()
    gamma = params.gamma
    mu = params.blood_viscosity

    # --- Poiseuille Resistance ---

    @testset "compute_resistances! — Poiseuille formula" begin
        tree = VascularTree("test", 10)
        # Single straight tube: length=10, radius=1
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        compute_resistances!(tree, mu)

        # R = 8 * mu * L / (pi * r^4)
        L = 10.0
        r = 1.0
        R_expected = 8.0 * mu * L / (π * r^4)
        @test tree.segments.resistance[1] ≈ R_expected rtol = 1e-10
    end

    @testset "compute_resistances! — radius scaling" begin
        tree = VascularTree("test", 10)
        # Two segments with different radii
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 2.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(1))

        compute_resistances!(tree, mu)

        R1 = 8.0 * mu * 5.0 / (π * 2.0^4)
        R2 = 8.0 * mu * 5.0 / (π * 1.0^4)
        @test tree.segments.resistance[1] ≈ R1 rtol = 1e-10
        @test tree.segments.resistance[2] ≈ R2 rtol = 1e-10

        # Halving radius should increase resistance by 2^4 = 16x (same length)
        @test tree.segments.resistance[2] / tree.segments.resistance[1] ≈ 16.0 rtol = 1e-10
    end

    @testset "compute_resistances! — AK kernel on many segments" begin
        # Build a chain of 500 segments to test AK resistance computation
        n = 500
        tree = VascularTree("test", n + 10)
        rng = MersenneTwister(42)

        # Build a chain: each segment connects to previous
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5 + rand(rng), Int32(-1))
        for i in 2:n
            r = 0.1 + rand(rng) * 2.0
            px = Float64(i - 1)
            dx = Float64(i) + rand(rng)
            add_segment!(tree, (px, 0.0, 0.0), (dx, rand(rng), 0.0), r, Int32(i - 1))
        end

        compute_resistances!(tree, mu)

        # Verify against scalar reference
        seg = tree.segments
        for i in 1:seg.n
            R_ref = 8.0 * mu * seg.seg_length[i] / (π * seg.radius[i]^4)
            @test seg.resistance[i] ≈ R_ref rtol = 1e-10
        end
    end

    @testset "compute_resistances! — all positive" begin
        tree = VascularTree("test", 50)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (8.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (8.0, -3.0, 0.0), 0.4, Int32(1))

        compute_resistances!(tree, mu)

        for i in 1:tree.segments.n
            @test tree.segments.resistance[i] > 0.0
        end
    end

    # --- Flow Computation ---

    @testset "compute_flows! — single segment" begin
        tree = VascularTree("test", 10)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        compute_resistances!(tree, mu)
        compute_flows!(tree, params)

        # Single terminal: flow = total_flow = (P_root - P_term) / R_total
        R = tree.segments.resistance[1]
        expected_flow = (params.root_pressure - params.terminal_pressure) / R
        @test tree.segments.flow[1] ≈ expected_flow rtol = 1e-6
        @test tree.segments.flow[1] > 0.0
    end

    @testset "compute_flows! — flow conservation at bifurcation" begin
        tree = VascularTree("test", 50)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.4, Int32(1))

        # Set radii via Murray's law
        update_radii!(tree, gamma)
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)

        seg = tree.segments
        topo = tree.topology

        # Flow conservation: Q_parent = Q_child1 + Q_child2
        c1 = topo.child1_id[1]
        c2 = topo.child2_id[1]
        @test seg.flow[1] ≈ seg.flow[c1] + seg.flow[c2] rtol = 1e-6
        @test seg.flow[1] > 0.0
    end

    @testset "compute_flows! — all flows positive" begin
        tree = VascularTree("test", 50)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.4, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 5.0, 0.0), 0.3, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 1.0, 0.0), 0.2, Int32(2))

        update_radii!(tree, gamma)
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)

        for i in 1:tree.segments.n
            @test tree.segments.flow[i] > 0.0
        end
    end

    @testset "compute_flows! — flow conservation at every bifurcation" begin
        tree = VascularTree("test", 100)
        # Build a 3-level tree
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 5.0, 0.0), 0.3, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 1.0, 0.0), 0.3, Int32(2))
        add_segment!(tree, (10.0, -3.0, 0.0), (14.0, -1.0, 0.0), 0.3, Int32(3))
        add_segment!(tree, (10.0, -3.0, 0.0), (14.0, -5.0, 0.0), 0.3, Int32(3))

        update_radii!(tree, gamma)
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)

        topo = tree.topology
        seg = tree.segments

        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                @test seg.flow[i] ≈ seg.flow[c1] + seg.flow[c2] rtol = 1e-6
            end
        end
    end

    # --- Pressure Computation ---

    @testset "compute_pressures! — single segment" begin
        tree = VascularTree("test", 10)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        # Root proximal = root_pressure
        @test tree.segments.pressure_proximal[1] ≈ params.root_pressure rtol = 1e-6
        # Distal = proximal - Q*R
        P_distal = params.root_pressure - tree.segments.flow[1] * tree.segments.resistance[1]
        @test tree.segments.pressure_distal[1] ≈ P_distal rtol = 1e-6
    end

    @testset "compute_pressures! — monotone decrease root to terminal" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.4, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 5.0, 0.0), 0.3, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 1.0, 0.0), 0.2, Int32(2))

        update_radii!(tree, gamma)
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        seg = tree.segments
        topo = tree.topology

        for i in 1:seg.n
            # Pressure drops along each segment
            @test seg.pressure_proximal[i] >= seg.pressure_distal[i]

            # Child proximal = parent distal
            if topo.parent_id[i] > 0
                parent = topo.parent_id[i]
                @test seg.pressure_proximal[i] ≈ seg.pressure_distal[parent] rtol = 1e-6
            end
        end
    end

    @testset "compute_pressures! — root pressure correct" begin
        tree = VascularTree("test", 50)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.4, Int32(1))

        update_radii!(tree, gamma)
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        @test tree.segments.pressure_proximal[tree.root_segment_id] ≈ params.root_pressure rtol = 1e-6
    end

    @testset "compute_pressures! — all pressures in valid range" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 5.0, 0.0), 0.3, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (14.0, 1.0, 0.0), 0.3, Int32(2))
        add_segment!(tree, (10.0, -3.0, 0.0), (14.0, -1.0, 0.0), 0.3, Int32(3))
        add_segment!(tree, (10.0, -3.0, 0.0), (14.0, -5.0, 0.0), 0.3, Int32(3))

        update_radii!(tree, gamma)
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        seg = tree.segments
        for i in 1:seg.n
            @test seg.pressure_proximal[i] >= 0.0
            @test seg.pressure_proximal[i] <= params.root_pressure * 1.01  # small tolerance
            @test seg.pressure_distal[i] >= 0.0
            @test seg.pressure_distal[i] <= params.root_pressure * 1.01
        end
    end

    # --- Validation ---

    @testset "validate_hemodynamics — valid tree" begin
        tree = VascularTree("test", 50)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.4, Int32(1))

        update_radii!(tree, gamma)
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        @test validate_hemodynamics(tree, params) == true
    end

    @testset "validate_hemodynamics — CCO-grown tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng)

        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        @test validate_hemodynamics(tree, params) == true
    end

    @testset "validate_hemodynamics — invalid: negative flow" begin
        tree = VascularTree("test", 10)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        # Corrupt flow
        tree.segments.flow[1] = -1.0
        @test validate_hemodynamics(tree, params) == false
    end

    @testset "validate_hemodynamics — invalid: negative pressure" begin
        tree = VascularTree("test", 10)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        # Corrupt pressure
        tree.segments.pressure_distal[1] = -100.0
        @test validate_hemodynamics(tree, params) == false
    end

    # --- Flow-weighted radius assignment ---

    @testset "assign_terminal_flows! — non-uniform flow" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 20, params; rng=rng)

        compute_resistances!(tree, mu)
        assign_terminal_flows!(tree, domain, params)

        seg = tree.segments
        topo = tree.topology

        # All terminal flows should be positive
        terminal_flows = Float64[]
        for i in 1:seg.n
            if topo.is_terminal[i]
                @test seg.flow[i] > 0.0
                push!(terminal_flows, seg.flow[i])
            end
        end

        # Flows should not all be identical (territory-weighted)
        if length(terminal_flows) > 2
            @test !all(f ≈ terminal_flows[1] for f in terminal_flows)
        end
    end

    @testset "assign_terminal_flows! — total flow conservation" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(77)
        grow_tree!(tree, domain, 30, params; rng=rng)

        compute_resistances!(tree, mu)
        assign_terminal_flows!(tree, domain, params)

        seg = tree.segments
        topo = tree.topology

        # Sum of terminal flows should equal root flow
        terminal_sum = 0.0
        for i in 1:seg.n
            if topo.is_terminal[i]
                terminal_sum += seg.flow[i]
            end
        end
        @test seg.flow[tree.root_segment_id] ≈ terminal_sum rtol = 1e-6
    end

    @testset "recompute_radii_from_flow! — Murray's law holds" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(88)
        grow_tree!(tree, domain, 20, params; rng=rng)

        compute_resistances!(tree, mu)
        assign_terminal_flows!(tree, domain, params)
        recompute_radii_from_flow!(tree, params)

        seg = tree.segments
        topo = tree.topology

        # Murray's law: r_parent^gamma = r_c1^gamma + r_c2^gamma
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

    @testset "recompute_radii_from_flow! — all radii positive" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(99)
        grow_tree!(tree, domain, 20, params; rng=rng)

        compute_resistances!(tree, mu)
        assign_terminal_flows!(tree, domain, params)
        recompute_radii_from_flow!(tree, params)

        for i in 1:tree.segments.n
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "Flow-weighted pipeline — validate_hemodynamics" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(111)
        grow_tree!(tree, domain, 30, params; rng=rng)

        compute_resistances!(tree, mu)
        assign_terminal_flows!(tree, domain, params)
        recompute_radii_from_flow!(tree, params)

        # Recompute everything with new radii
        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        @test validate_hemodynamics(tree, params) == true
    end

    # --- Integration with CCO-grown trees ---

    @testset "Full pipeline on CCO tree — 50 terminals" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(999)
        grow_tree!(tree, domain, 50, params; rng=rng)

        compute_resistances!(tree, mu)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)

        seg = tree.segments
        topo = tree.topology
        n = seg.n

        # All flows positive
        for i in 1:n
            @test seg.flow[i] > 0.0
        end

        # Flow conservation at every bifurcation
        for i in 1:n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                if c1 > 0 && c2 > 0
                    @test seg.flow[i] ≈ seg.flow[c1] + seg.flow[c2] rtol = 1e-4
                end
            end
        end

        # Pressures decrease monotonically
        for i in 1:n
            @test seg.pressure_proximal[i] >= seg.pressure_distal[i] - 1e-6
        end

        # Root pressure is correct
        @test seg.pressure_proximal[tree.root_segment_id] ≈ params.root_pressure rtol = 1e-6

        # Validation passes
        @test validate_hemodynamics(tree, params) == true
    end

end
