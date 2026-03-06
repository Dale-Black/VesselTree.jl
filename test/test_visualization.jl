using Test
using VesselTree

@testset "Visualization (stubs + extension)" begin

    @testset "Function stubs exist" begin
        @test isdefined(VesselTree, :plot_tree)
        @test isdefined(VesselTree, :plot_tree_2d)
        @test isdefined(VesselTree, :plot_validation_report)
        @test isdefined(VesselTree, :plot_forest)
    end

    # Extension tests only run if CairoMakie is available
    has_makie = try
        @eval using CairoMakie
        true
    catch
        false
    end

    if has_makie
        @testset "plot_tree renders" begin
            tree = VascularTree("viz", 100)
            add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
            add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
            add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.3, Int32(1))
            fig = plot_tree(tree)
            @test fig isa Figure
        end

        @testset "plot_tree_2d renders" begin
            tree = VascularTree("viz2d", 100)
            add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
            add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
            fig = plot_tree_2d(tree)
            @test fig isa Figure
            fig2 = plot_tree_2d(tree; projection=:xz)
            @test fig2 isa Figure
        end

        @testset "plot_validation_report renders" begin
            tree = VascularTree("vizval", 100)
            add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
            add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
            add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.3, Int32(1))
            params = kassab_coronary_params()
            fig = plot_validation_report(tree, params)
            @test fig isa Figure
        end
    else
        @info "CairoMakie not available — skipping extension tests"
    end
end
