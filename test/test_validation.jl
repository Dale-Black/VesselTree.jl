using Test
using VesselTree
using Random

@testset "Validation Framework" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    # --- ValidationReport struct ---

    @testset "ValidationReport — construction" begin
        report = ValidationReport()
        @test report isa ValidationReport
        @test length(report.diameter_ks_pvalues) == 0
        @test report.asymmetry_median == 0.0
        @test report.ld_ratio_mean == 0.0
        @test report.angle_mean == 0.0
    end

    # --- Diameter extraction via AK ---

    @testset "compute_diameters! — basic" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 0.5, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.3, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.2, Int32(1))

        diameters_buf = zeros(100)
        VesselTree.compute_diameters!(diameters_buf, tree.segments)

        @test diameters_buf[1] ≈ 1.0 rtol = 1e-10  # 2 * 0.5
        @test diameters_buf[2] ≈ 0.6 rtol = 1e-10  # 2 * 0.3
        @test diameters_buf[3] ≈ 0.4 rtol = 1e-10  # 2 * 0.2
    end

    # --- Asymmetry ratio computation ---

    @testset "compute_asymmetry_ratios — basic" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.8, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.4, Int32(1))

        ratios = VesselTree.compute_asymmetry_ratios(tree)
        @test length(ratios) == 1
        @test ratios[1] ≈ 0.4 / 0.8 rtol = 1e-10  # r_small / r_large = 0.5
    end

    @testset "compute_asymmetry_ratios — multiple bifurcations" begin
        tree = VascularTree("test", 200)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.8, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.6, Int32(1))
        # Second bifurcation
        add_segment!(tree, (15.0, 1.0, 0.0), (20.0, 2.0, 0.0), 0.5, Int32(2))
        add_segment!(tree, (15.0, 1.0, 0.0), (20.0, 0.0, 0.0), 0.3, Int32(2))

        ratios = VesselTree.compute_asymmetry_ratios(tree)
        @test length(ratios) == 2
        # All ratios should be in [0, 1]
        for r in ratios
            @test 0.0 <= r <= 1.0
        end
    end

    # --- L/D ratio ---

    @testset "compute_ld_ratios — basic" begin
        tree = VascularTree("test", 100)
        # Segment with length 10, radius 1 -> L/D = 10/2 = 5
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        ld_buf = zeros(100)
        VesselTree.compute_ld_ratios!(ld_buf, tree.segments)
        @test ld_buf[1] ≈ 5.0 rtol = 1e-10  # 10.0 / (2*1.0) = 5.0
    end

    # --- Branching angles ---

    @testset "compute_branching_angles — basic" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        # Two daughters at 45 degrees from parent
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 5.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -5.0, 0.0), 0.5, Int32(1))

        angles = VesselTree.compute_branching_angles(tree)
        @test length(angles) == 1
        # Angle between the two daughter vectors (45deg + 45deg = 90deg)
        @test angles[1] ≈ π / 2 atol = 0.01
    end

    @testset "compute_branching_angles — straight continuation" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (20.0, 0.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (10.0, 5.0, 0.0), 0.5, Int32(1))

        angles = VesselTree.compute_branching_angles(tree)
        @test length(angles) == 1
        # 90 degree angle
        @test angles[1] ≈ π / 2 atol = 0.01
    end

    # --- validate_tree integration ---

    @testset "validate_tree — on small tree" begin
        tree = VascularTree("test", 200)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.8, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.4, Int32(1))

        report = validate_tree(tree, params)

        @test report isa ValidationReport
        @test report.n_segments == 3
        @test report.n_bifurcations >= 1
        @test report.ld_ratio_mean > 0.0
    end

    @testset "validate_tree — on CCO-grown tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng, kassab=true)

        report = validate_tree(tree, params)

        @test report.n_segments > 10
        @test report.n_bifurcations > 0
        @test report.ld_ratio_mean > 0.0
        @test report.asymmetry_median > 0.0
        @test report.asymmetry_median <= 1.0
        @test length(report.branching_angles) > 0
    end

    # --- print_report ---

    @testset "print_report — does not error" begin
        report = ValidationReport()
        report = ValidationReport(
            n_segments=10,
            n_bifurcations=5,
            n_trifurcations=0,
            diameter_ks_pvalues=Dict(0 => 0.5, 1 => 0.3),
            asymmetry_ratios=Float64[0.5, 0.6, 0.7],
            asymmetry_median=0.6,
            ld_ratios=Float64[5.0, 6.0],
            ld_ratio_mean=5.5,
            ld_ratio_std=0.7,
            branching_angles=Float64[1.0, 1.2],
            angle_mean=1.1,
            angle_std=0.14,
        )
        # Should not throw
        buf = IOBuffer()
        print_report(buf, report)
        output = String(take!(buf))
        @test length(output) > 0
        @test contains(output, "Segments")
    end

    # --- Bulk AK computations ---

    @testset "compute_diameters! — uses AK (large array)" begin
        tree = VascularTree("test", 1000)
        for i in 1:100
            r = 0.1 + 0.01 * i
            add_segment!(tree, (Float64(i), 0.0, 0.0), (Float64(i) + 1.0, 0.0, 0.0), r, Int32(-1))
            # Note: adding all as roots is fine for testing diameter computation
        end
        # Fix: actually tree only allows one root, so build a chain
        tree2 = VascularTree("test", 1000)
        add_segment!(tree2, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5, Int32(-1))
        for i in 2:100
            r = 0.1 + 0.005 * i
            add_segment!(tree2, (Float64(i-1), 0.0, 0.0), (Float64(i), 0.0, 0.0), r, Int32(Int32(i-1)))
        end

        diameters_buf = zeros(1000)
        VesselTree.compute_diameters!(diameters_buf, tree2.segments)

        # All diameters should be 2*radius
        for i in 1:tree2.segments.n
            @test diameters_buf[i] ≈ 2.0 * tree2.segments.radius[i] rtol = 1e-10
        end
    end

end
