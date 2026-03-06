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
    end

    @testset "Diameter tables" begin
        p = kassab_coronary_params()
        @test length(p.diameter_mean) == 12
        @test length(p.diameter_sd) == 12

        # Order 0 = capillary = 8 um
        @test p.diameter_mean[1] == 8.0
        # Order 11 = stem = 4500 um
        @test p.diameter_mean[12] == 4500.0

        # Diameters increase monotonically
        for i in 2:12
            @test p.diameter_mean[i] > p.diameter_mean[i-1]
        end

        # SD is positive
        @test all(p.diameter_sd .> 0)
    end

    @testset "Length tables" begin
        p = kassab_coronary_params()
        @test length(p.length_mean) == 12
        @test length(p.length_sd) == 12

        # Order 0 capillary length = 20 um
        @test p.length_mean[1] == 20.0
        # Order 11 stem length = 150000 um
        @test p.length_mean[12] == 150000.0

        # Lengths increase monotonically
        for i in 2:12
            @test p.length_mean[i] > p.length_mean[i-1]
        end
    end

    @testset "Diameter boundaries" begin
        p = kassab_coronary_params()
        @test length(p.diameter_bounds) == 13  # n_orders + 1

        # First bound is 5 um
        @test p.diameter_bounds[1] == 5.0
        # Last bound is Inf
        @test p.diameter_bounds[13] == Inf

        # Bounds increase monotonically
        for i in 2:12
            @test p.diameter_bounds[i] > p.diameter_bounds[i-1]
        end
    end

    @testset "Connectivity matrix" begin
        p = kassab_coronary_params()
        @test size(p.connectivity_matrix) == (12, 12)

        # Column 1 (parent order 0) should be all zeros — order 0 has no children
        @test all(p.connectivity_matrix[:, 1] .== 0.0)

        # Diagonal-adjacent entries should be positive for higher orders
        @test p.connectivity_matrix[1, 2] > 0  # order 0 daughters of order 1 parent
        @test p.connectivity_matrix[2, 3] > 0  # order 1 daughters of order 2 parent

        # No negative entries
        @test all(p.connectivity_matrix .>= 0.0)

        # Column sums for parent orders 1-11 should be roughly 2-3
        for col in 2:12
            col_sum = sum(p.connectivity_matrix[:, col])
            @test col_sum > 1.0   # at least one daughter
            @test col_sum < 5.0   # not too many
        end
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
        p = kassab_coronary_params()

        # Order 0: 5-11 um
        @test classify_order(p, 8.0) == 0
        @test classify_order(p, 5.0) == 0
        @test classify_order(p, 10.0) == 0

        # Order 1: 11-21 um
        @test classify_order(p, 15.0) == 1
        @test classify_order(p, 11.0) == 1

        # Order 5: 140-250 um
        @test classify_order(p, 184.0) == 5

        # Order 11: >= 4500 um
        @test classify_order(p, 4500.0) == 11
        @test classify_order(p, 10000.0) == 11

        # Below minimum
        @test classify_order(p, 4.0) == -1
    end

end
