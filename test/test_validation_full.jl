using Test
using VesselTree
using Random

@testset "Full Validation Suite" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    # --- Murray's law deviation ---

    @testset "compute_murray_deviation — perfect tree" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.7, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))

        # Force Murray consistency
        update_radii!(tree, gamma)

        max_dev, mean_dev = VesselTree.compute_murray_deviation(tree, gamma)
        @test max_dev < 1e-6
        @test mean_dev < 1e-6
    end

    @testset "compute_murray_deviation — imperfect tree" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.7, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))
        # Don't update radii — parent r=1.0 may not match sum of children

        max_dev, mean_dev = VesselTree.compute_murray_deviation(tree, gamma)
        @test max_dev >= 0.0
        @test mean_dev >= 0.0
    end

    # --- Trifurcation percentage ---

    @testset "compute_trifurcation_pct — no trifurcations" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))

        pct = VesselTree.compute_trifurcation_pct(tree)
        @test pct ≈ 0.0 atol = 1e-10
    end

    @testset "compute_trifurcation_pct — with trifurcation" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 0.0, 1.0), 0.3, Int32(1))

        pct = VesselTree.compute_trifurcation_pct(tree)
        @test pct ≈ 100.0 atol = 1e-10  # 1 trifurcation, 0 bifurcations
    end

    # --- P(lambda) distribution ---

    @testset "compute_lambda_distribution — basic" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.5, Int32(1))

        lambdas = VesselTree.compute_lambda_distribution(tree)
        @test length(lambdas) == 1  # one bifurcation
        @test lambdas[1] > 0.0
    end

    @testset "compute_lambda_distribution — values positive" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 1000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true)

        lambdas = VesselTree.compute_lambda_distribution(tree)
        @test length(lambdas) > 0
        for l in lambdas
            @test l > 0.0
        end
    end

    # --- Connectivity matrix in validation ---

    @testset "validate_connectivity_in_report — computed" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng, kassab=true)

        report = validate_tree(tree, params)
        @test haskey(report.extra, :connectivity_chi2)
        @test haskey(report.extra, :trifurcation_pct)
        @test haskey(report.extra, :murray_max_deviation)
        @test haskey(report.extra, :murray_mean_deviation)
    end

    # --- Full validate_tree with all metrics ---

    @testset "validate_tree — comprehensive" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 50, params; rng=rng, kassab=true)
        update_radii!(tree, gamma)

        report = validate_tree(tree, params)

        # Murray's law should be satisfied after update_radii!
        @test report.extra[:murray_max_deviation] < 1e-4
        @test report.extra[:murray_mean_deviation] < 1e-4

        # Should have lambda values
        @test haskey(report.extra, :lambda_values)
        @test length(report.extra[:lambda_values]) > 0

        # Trifurcation pct should be in valid range
        @test report.extra[:trifurcation_pct] >= 0.0
        @test report.extra[:trifurcation_pct] <= 100.0
    end

    # --- print_report with extra metrics ---

    @testset "print_report — includes new metrics" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 2000)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 30, params; rng=rng, kassab=true)

        report = validate_tree(tree, params)
        buf = IOBuffer()
        print_report(buf, report)
        output = String(take!(buf))

        @test contains(output, "Murray")
        @test contains(output, "Trifurcation")
    end

end
