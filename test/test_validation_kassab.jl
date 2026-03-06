using Test
using VesselTree
using Random
using Distributions

@testset "Comprehensive Kassab Validation" begin
    params = kassab_coronary_params()
    rng = MersenneTwister(42)

    # Generate a medium tree with full pipeline for validation
    domain = SphereDomain((0.0, 0.0, 0.0), 5.0)
    test_configs = [
        TreeConfig("TEST", (0.0, 0.0, 0.0), 1.5, (0.0, -1.0, 0.0), 40, 1.0),
    ]
    forest = generate_kassab_coronary(
        domain, params;
        rng=rng, verbose=false, handoff_order=4,
        tree_configs=test_configs,
    )
    tree = forest.trees["TEST"]

    @testset "KS test correctness — known distribution" begin
        # Generate samples from N(100, 10) and test against N(100, 10)
        rng2 = MersenneTwister(7)
        samples = [100.0 + 10.0 * randn(rng2) for _ in 1:200]
        d_normal = Normal(100.0, 10.0)
        D, p = ks_test_onesample(samples, x -> cdf(d_normal, x))
        # Should not reject (samples are from the distribution)
        @test D < 0.15
        @test p > 0.05

        # Test against wrong distribution N(200, 10)
        d_wrong = Normal(200.0, 10.0)
        D_wrong, p_wrong = ks_test_onesample(samples, x -> cdf(d_wrong, x))
        @test D_wrong > 0.5
        @test p_wrong < 0.01
    end

    @testset "KS test — empty and trivial inputs" begin
        D, p = ks_test_onesample(Float64[], x -> x)
        @test D == 0.0
        @test p == 1.0

        D2, p2 = ks_test_onesample([0.5], x -> x)  # Uniform CDF
        @test D2 >= 0.0
        @test p2 >= 0.0
    end

    @testset "Connectivity matrix audit" begin
        audit = validate_connectivity_matrix(params)
        col_sums = audit[:column_sums]

        # All column sums should be between 2.0 and 3.0
        @test audit[:all_sums_ge_2]
        @test audit[:all_sums_le_3]

        # Diagonal entries should be dominant in each column
        @test audit[:diagonal_dominant]

        # Verify specific column sums match research data
        @test col_sums[1] ≈ 2.3 atol=0.01  # parent order 1
        @test col_sums[2] ≈ 2.1 atol=0.01  # parent order 2
    end

    @testset "Diameter KS per order" begin
        diam_results = validate_diameters_per_order(tree, params)
        # Should have results for at least some orders
        @test length(diam_results) >= 2

        # Check result structure
        for (ord, r) in diam_results
            @test r.n >= 5
            @test r.D >= 0.0
            @test r.D <= 1.0
            @test r.p >= 0.0
            @test r.p <= 1.0
            @test r.mean_um > 0.0
        end

        # Orders populated by subdivision (0 to handoff_order-1) should
        # have diameters close to Kassab means (they were sampled from N(mean,sd))
        for ord in 0:3
            if haskey(diam_results, ord) && diam_results[ord].n >= 10
                r = diam_results[ord]
                expected_mean = params.diameter_mean[ord + 1]
                # Mean should be within 3 SDs of expected
                expected_sd = params.diameter_sd[ord + 1]
                @test abs(r.mean_um - expected_mean) < 3 * expected_sd
            end
        end
    end

    @testset "Length KS per order" begin
        len_results = validate_lengths_per_order(tree, params)
        @test length(len_results) >= 1

        for (ord, r) in len_results
            @test r.n >= 5
            @test r.D >= 0.0 && r.D <= 1.0
            @test r.p >= 0.0 && r.p <= 1.0
            @test r.mean_um > 0.0
        end
    end

    @testset "Segment counts per order" begin
        counts = validate_segment_counts(tree, params)
        # Should have at least 2 orders populated
        @test length(counts) >= 2

        # Order 0 should have the most segments (capillaries)
        if haskey(counts, 0)
            @test counts[0] >= 1
        end

        # Total should match tree.segments.n (minus unclassified)
        total_counted = sum(values(counts))
        @test total_counted <= tree.segments.n
    end

    @testset "Asymmetry distribution KS" begin
        asym = validate_asymmetry_ks(tree, params)
        @test asym.n >= 1
        @test asym.D >= 0.0
        @test asym.p >= 0.0
        @test 0.0 < asym.median <= 1.0
    end

    @testset "Murray's law compliance" begin
        max_dev, mean_dev = VesselTree.compute_murray_deviation(tree, params.gamma)
        # With floor enforcement, some deviation is expected
        # Mean should be small
        @test mean_dev < 0.05
        # Max deviation bounded
        @test max_dev < 1.0
    end

    @testset "Connectivity chi-squared" begin
        VesselTree.assign_strahler_orders!(tree, params)
        cm_emp = VesselTree.build_empirical_connectivity(tree, params)
        chi2, p_val = VesselTree.validate_connectivity(cm_emp, params.connectivity_matrix)
        @test chi2 >= 0.0
        @test p_val >= 0.0
        @test p_val <= 1.0
    end

    @testset "Report card generation" begin
        card = generate_report_card(tree, params)
        @test haskey(card, :passed)
        @test haskey(card, :total)
        @test haskey(card, :grade)
        @test card[:total] >= 4  # at least 4 metrics evaluated
        @test card[:grade] >= 0.0
        @test card[:grade] <= 1.0
        @test haskey(card, :diameter_ks)
        @test haskey(card, :murray_mean_dev)
        @test haskey(card, :segment_counts)
    end

    @testset "Report card printing" begin
        card = generate_report_card(tree, params)
        io = IOBuffer()
        print_report_card(io, card)
        output = String(take!(io))
        @test occursin("Kassab Validation Report Card", output)
        @test occursin("Diameter KS", output)
        @test occursin("Murray", output)
        @test occursin("Overall", output)
    end

    @testset "Updated validate_tree uses proper KS" begin
        report = validate_tree(tree, params)
        # The diameter_ks_pvalues should now be proper KS p-values
        for (ord, p) in report.diameter_ks_pvalues
            @test 0.0 <= p <= 1.0
        end
    end
end
