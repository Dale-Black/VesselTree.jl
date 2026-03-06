using Test
using VesselTree
using Distributions

@testset "Morphometric Parameters" begin

    @testset "kassab_coronary_params construction" begin
        p = kassab_coronary_params()
        @test p isa MorphometricParams
        @test p.n_orders == 12
        @test p.gamma == 7 / 3
        @test p.vessel_cutoff_um == 8.0
        @test p.artery_name == "RCA"
    end

    @testset "Per-artery param functions" begin
        rca = kassab_rca_params()
        lad = kassab_lad_params()
        lcx = kassab_lcx_params()

        @test rca.artery_name == "RCA"
        @test lad.artery_name == "LAD"
        @test lcx.artery_name == "LCX"

        @test rca.n_orders == 12  # orders 0-11
        @test lad.n_orders == 12
        @test lcx.n_orders == 11  # orders 0-10 (LCX has fewer)

        # kassab_coronary_params is alias for RCA
        @test kassab_coronary_params().diameter_mean == rca.diameter_mean
    end

    @testset "RCA segment diameters (Table 1)" begin
        p = kassab_rca_params()
        @test length(p.diameter_mean) == 12
        @test length(p.diameter_sd) == 12

        # Order 0 = capillary = 8 um (estimated)
        @test p.diameter_mean[1] == 8.0
        # Order 1 = 9.6 um (Kassab 1993 Table 1)
        @test p.diameter_mean[2] == 9.6
        # Order 11 = stem = 3218 um (Kassab 1993 Table 1)
        @test p.diameter_mean[12] == 3218.0

        # Diameters increase monotonically
        for i in 2:12
            @test p.diameter_mean[i] > p.diameter_mean[i-1]
        end

        # SD is positive
        @test all(p.diameter_sd .> 0)
    end

    @testset "RCA segment lengths (Table 1)" begin
        p = kassab_rca_params()
        @test length(p.length_mean) == 12
        @test length(p.length_sd) == 12

        # Order 1 segment length = 69 um (Kassab 1993 Table 1: 0.069 mm)
        @test p.length_mean[2] == 69.0
        # Order 11 segment length = 3240 um
        @test p.length_mean[12] == 3240.0

        # All positive
        @test all(p.length_mean .> 0)
    end

    @testset "Element-level data" begin
        p = kassab_rca_params()

        # Element data arrays exist and are correct size
        @test length(p.diameter_mean_elem) == 12
        @test length(p.diameter_sd_elem) == 12
        @test length(p.length_mean_elem) == 12
        @test length(p.length_sd_elem) == 12

        # Element diameters are slightly smaller (averaging reduces outliers)
        # RCA order 1 element D=9.3 vs segment D=9.6
        @test p.diameter_mean_elem[2] == 9.3

        # Element lengths are longer (sum of constituent segments)
        # RCA order 1: element L=125 vs segment L=69
        @test p.length_mean_elem[2] == 125.0
        @test p.length_mean_elem[2] > p.length_mean[2]
    end

    @testset "S/E ratios (Table 5)" begin
        p = kassab_rca_params()
        @test length(p.se_ratio) == 12

        # Order 0 capillaries = 1 segment per element
        @test p.se_ratio[1] == 1.0
        # RCA order 1 S/E = 1.88
        @test p.se_ratio[2] == 1.88
        # RCA order 11 S/E = 26
        @test p.se_ratio[12] == 26.0
        # All >= 1
        @test all(p.se_ratio .>= 1.0)
    end

    @testset "Element count targets (Table 9)" begin
        p = kassab_rca_params()
        @test length(p.element_count_target) == 12

        # Order 0 not measured
        @test p.element_count_target[1] == 0.0
        # RCA order 1 = 393,294 elements
        @test p.element_count_target[2] == 393294.0
        # RCA order 11 = 1 element (the stem)
        @test p.element_count_target[12] == 1.0
        # Counts decrease with order (excluding order 0)
        for i in 3:12
            @test p.element_count_target[i] < p.element_count_target[i-1]
        end
    end

    @testset "Diameter boundaries" begin
        p = kassab_coronary_params()
        @test length(p.diameter_bounds) == 13  # n_orders + 1

        # First bound is 0
        @test p.diameter_bounds[1] == 0.0
        # Last bound is Inf
        @test p.diameter_bounds[13] == Inf

        # Bounds increase monotonically
        for i in 2:12
            @test p.diameter_bounds[i] > p.diameter_bounds[i-1]
        end
    end

    @testset "Connectivity matrix (Table 6 — RCA)" begin
        p = kassab_rca_params()
        @test size(p.connectivity_matrix) == (12, 12)

        # Column 1 (parent order 0) should be all zeros — order 0 has no children
        @test all(p.connectivity_matrix[:, 1] .== 0.0)

        # Verify specific values from Kassab 1993 Table 6
        @test p.connectivity_matrix[1, 2] == 2.75   # order 0 from order 1
        @test p.connectivity_matrix[2, 3] == 2.13   # order 1 from order 2
        @test p.connectivity_matrix[6, 7] == 2.43   # order 5 from order 6

        # No negative entries
        @test all(p.connectivity_matrix .>= 0.0)

        # Column sums for parent orders 1-11 should be > 1
        for col in 2:12
            col_sum = sum(p.connectivity_matrix[:, col])
            @test col_sum > 1.0
        end
    end

    @testset "LCX connectivity matrix (Table 8)" begin
        p = kassab_lcx_params()
        @test size(p.connectivity_matrix) == (11, 11)

        # Column 1 (parent order 0) should be all zeros
        @test all(p.connectivity_matrix[:, 1] .== 0.0)

        # Verify specific LCX values from Table 8
        @test p.connectivity_matrix[1, 2] == 3.18   # order 0 from order 1 (pooled with LAD)
        @test p.connectivity_matrix[6, 7] == 2.51   # order 5 from order 6
    end

    @testset "Asymmetry distribution" begin
        p = kassab_coronary_params()
        @test p.asymmetry_alpha == 2.5
        @test p.asymmetry_beta == 0.8

        # Mean of Beta(2.5, 0.8) ≈ 0.76
        d = Beta(p.asymmetry_alpha, p.asymmetry_beta)
        @test mean(d) ≈ 0.7576 atol = 0.01
    end

    @testset "Barabasi thresholds" begin
        p = kassab_coronary_params()
        @test p.trifurcation_chi_th == 0.83
        @test p.sprouting_rho_th == 0.83
    end

    @testset "Hemodynamic parameters" begin
        p = kassab_coronary_params()
        @test p.blood_viscosity == 0.0035
        @test p.root_pressure ≈ 13332.0 atol = 1.0
        @test p.terminal_pressure ≈ 3999.6 atol = 1.0
        @test p.root_pressure > p.terminal_pressure
    end

    @testset "classify_order" begin
        p = kassab_rca_params()

        # Order 0: 0-8.815 um (RCA bounds)
        @test classify_order(p, 8.0) == 0
        @test classify_order(p, 5.0) == 0

        # Order 1: 8.815-11.085 um
        @test classify_order(p, 9.6) == 1

        # Order 5: 44.7-92.5 um
        @test classify_order(p, 64.4) == 5

        # Order 11: >= 2319.5 um
        @test classify_order(p, 3218.0) == 11
        @test classify_order(p, 10000.0) == 11
    end

end
