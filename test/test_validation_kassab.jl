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

        # Verify specific column sums match real Kassab 1993 Table 6 (RCA)
        @test col_sums[1] ≈ 2.881 atol=0.01  # parent order 1: CM[1,2]+CM[2,2]
        @test col_sums[2] ≈ 2.884 atol=0.01  # parent order 2: CM[1,3]+CM[2,3]+CM[3,3]
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

        # Orders populated by subdivision should have element diameters
        # close to Kassab element-level means
        for ord in 0:3
            if haskey(diam_results, ord) && diam_results[ord].n >= 10
                r = diam_results[ord]
                expected_mean = params.diameter_mean_elem[ord + 1]
                expected_sd = params.diameter_sd_elem[ord + 1]
                @test abs(r.mean_um - expected_mean) < 5 * expected_sd
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

    @testset "Report card generation — 9 metrics" begin
        card = generate_report_card(tree, params)
        @test haskey(card, :passed)
        @test haskey(card, :total)
        @test haskey(card, :grade)
        @test card[:total] >= 9  # all 9 metrics evaluated
        @test card[:grade] >= 0.0
        @test card[:grade] <= 1.0
        @test haskey(card, :diameter_ks)
        @test haskey(card, :length_ks)
        @test haskey(card, :connectivity_chi2)
        @test haskey(card, :asymmetry)
        @test haskey(card, :murray_mean_dev)
        @test haskey(card, :trifurcation_pct)
        @test haskey(card, :se_ratios)
        @test haskey(card, :element_counts)
        @test haskey(card, :segment_counts)
    end

    @testset "Report card printing — all 9 metrics shown" begin
        card = generate_report_card(tree, params)
        io = IOBuffer()
        print_report_card(io, card)
        output = String(take!(io))
        @test occursin("Kassab Validation Report Card", output)
        @test occursin("Diameter KS", output)
        @test occursin("Length KS", output)
        @test occursin("Connectivity CM", output)
        @test occursin("Asymmetry", output)
        @test occursin("Murray", output)
        @test occursin("Trifurcation", output)
        @test occursin("S/E ratios", output)
        @test occursin("Element counts", output)
        @test occursin("Orders populated", output)
        @test occursin("Overall", output)
    end

    @testset "Element-level diameter KS" begin
        diam_results = validate_diameters_per_order(tree, params; element_level=true)
        @test length(diam_results) >= 1
        for (ord, r) in diam_results
            @test r.n >= 5
            @test 0.0 <= r.D <= 1.0
            @test 0.0 <= r.p <= 1.0
            @test r.mean_um > 0.0
        end
    end

    @testset "Element-level vs segment-level diameter KS" begin
        seg_results = validate_diameters_per_order(tree, params; element_level=false)
        elem_results = validate_diameters_per_order(tree, params; element_level=true)
        # Both should produce results
        @test length(seg_results) >= 1
        @test length(elem_results) >= 1
    end

    @testset "Element-level length KS" begin
        len_results = validate_lengths_per_order(tree, params; element_level=true)
        @test length(len_results) >= 1
        for (ord, r) in len_results
            @test r.n >= 5
            @test 0.0 <= r.D <= 1.0
            @test 0.0 <= r.p <= 1.0
            @test r.mean_um > 0.0
        end
    end

    @testset "Element count validation" begin
        elem_counts = validate_element_counts(tree, params)
        @test length(elem_counts) >= 1
        for (ord, count) in elem_counts
            @test ord >= 0
            @test count >= 1
        end
    end

    @testset "S/E ratio validation" begin
        se_results = validate_se_ratios(tree, params)
        @test length(se_results) >= 1
        for (ord, r) in se_results
            @test r.measured >= 1.0  # S/E always >= 1
            @test r.reference >= 1.0
            @test r.ratio_error >= 0.0
        end
    end

    @testset "Element-level asymmetry" begin
        asym = validate_asymmetry_ks(tree, params)
        @test asym.n >= 0
        if asym.n >= 5
            @test 0.0 < asym.median <= 1.0
        end
    end

    @testset "Updated validate_tree uses proper KS" begin
        report = validate_tree(tree, params)
        for (ord, p) in report.diameter_ks_pvalues
            @test 0.0 <= p <= 1.0
        end
    end

    @testset "Report card — trifurcation metric present" begin
        card = generate_report_card(tree, params)
        @test haskey(card, :trifurcation_pct)
        @test card[:trifurcation_pct] >= 0.0
    end

    @testset "Report card — connectivity uses element CM" begin
        card = generate_report_card(tree, params)
        @test haskey(card, :connectivity_chi2)
        @test card[:connectivity_chi2] >= 0.0
        @test 0.0 <= card[:connectivity_pval] <= 1.0
    end

    @testset "Report card — orders populated metric" begin
        card = generate_report_card(tree, params)
        @test card[:n_orders_populated] >= 2
    end

    @testset "validate_diameters_per_order — element_level=false backward compat" begin
        results = validate_diameters_per_order(tree, params; element_level=false)
        for (ord, r) in results
            expected_mean = params.diameter_mean[ord + 1]
            expected_sd = params.diameter_sd[ord + 1]
            @test abs(r.mean_um - expected_mean) < 5 * expected_sd
        end
    end

    @testset "validate_lengths_per_order — element_level=false backward compat" begin
        results = validate_lengths_per_order(tree, params; element_level=false)
        for (ord, r) in results
            @test r.mean_um > 0.0
        end
    end
end
