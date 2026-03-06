using Test
using VesselTree
import AcceleratedKernels as AK

@testset "VesselTree.jl" begin

    @testset "Module loading" begin
        @test VesselTree isa Module
        @test isdefined(VesselTree, :MURRAY_GAMMA)
        @test isdefined(VesselTree, :AK)
    end

    @testset "Scientific constants" begin
        # Murray exponent is 7/3, NOT 3.0
        @test VesselTree.MURRAY_GAMMA == 7 / 3
        @test VesselTree.MURRAY_GAMMA != 3.0
        @test VesselTree.MURRAY_GAMMA ≈ 2.333 atol = 0.001

        # Vessel cutoff is 8 um (capillary level)
        @test VesselTree.VESSEL_CUTOFF_UM == 8.0

        # Blood viscosity is 3.5 cP = 0.0035 Pa*s
        @test VesselTree.BLOOD_VISCOSITY == 0.0035

        # Root pressure is 100 mmHg
        @test VesselTree.ROOT_PRESSURE ≈ 13332.0 atol = 1.0

        # Terminal pressure is 30 mmHg
        @test VesselTree.TERMINAL_PRESSURE ≈ 3999.6 atol = 1.0

        # Pressure drop is positive
        @test VesselTree.ROOT_PRESSURE > VesselTree.TERMINAL_PRESSURE

        # Barabasi thresholds
        @test VesselTree.TRIFURCATION_CHI_TH == 0.83
        @test VesselTree.SPROUTING_RHO_TH == 0.83

        # Asymmetry distribution parameters
        @test VesselTree.ASYMMETRY_ALPHA == 2.5
        @test VesselTree.ASYMMETRY_BETA == 0.8
    end

    @testset "AcceleratedKernels availability" begin
        # Verify AK is accessible through VesselTree
        v = rand(100)
        out = similar(v)

        # Test foreachindex works
        AK.foreachindex(out) do i
            out[i] = v[i] * 2.0
        end
        @test all(out .≈ v .* 2.0)

        # Test reduce works
        s = AK.reduce(+, v; init=0.0)
        @test s ≈ sum(v)

        # Test minimum/maximum work
        @test AK.minimum(v) == minimum(v)
        @test AK.maximum(v) == maximum(v)

        # Test sort works
        sorted = AK.sort(v)
        @test issorted(sorted)

        # Test any/all work
        bools = Vector{Bool}(v .> 0.5)
        @test AK.any(identity, bools) == any(bools)
    end

end

include("test_types.jl")
include("test_parameters.jl")
include("test_domain.jl")
include("test_distance.jl")
include("test_intersection.jl")
include("test_kamiya.jl")
include("test_growth.jl")
include("test_spatial.jl")
include("test_hemodynamics.jl")
include("test_kassab.jl")
include("test_barabasi.jl")
include("test_trifurcation.jl")
include("test_surface_cost.jl")
include("test_validation.jl")
include("test_validation_full.jl")
include("test_forest.jl")
include("test_forest_growth.jl")
include("test_grid_growth.jl")
include("test_subdivision.jl")
include("test_kassab_refinement.jl")
