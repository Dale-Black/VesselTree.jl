using Test
using VesselTree

@testset "Contrast transport" begin
    @testset "single segment step response" begin
        tree = VascularTree("single", 4)
        length_mm = 1000.0 / π
        seg_id = add_segment!(tree, (0.0, 0.0, 0.0), (length_mm, 0.0, 0.0), 1.0, Int32(-1))
        tree.segments.flow[seg_id] = 1e-6  # 1 mL/s

        times = collect(0.0:0.1:3.0)
        root_curve = fill(10.0, length(times))
        result = simulate_contrast_transport(
            tree,
            times;
            root_input=root_curve,
            recompute_hemodynamics=false,
            return_outlet=true,
            storage_type=Float64,
        )

        mid_idx = findfirst(==(0.5), times)
        late_idx = findfirst(==(1.5), times)
        @test result.transit_time_s[1] ≈ 1.0 atol = 1e-6
        @test result.segment_volume_mL[1] ≈ 1.0 atol = 1e-6
        @test result.concentration[1, mid_idx] ≈ 5.0 atol = 1e-6
        @test result.concentration[1, late_idx] ≈ 10.0 atol = 1e-6
        @test result.outlet_concentration[1, findfirst(==(0.9), times)] ≈ 0.0 atol = 1e-6
        @test result.outlet_concentration[1, findfirst(==(1.0), times)] ≈ 10.0 atol = 1e-6
        @test all(result.concentration .>= 0.0)
    end

    @testset "bifurcation daughters receive identical concentration history" begin
        tree = VascularTree("bifurcation", 8)
        root_len = 500.0 / π
        child_len = 1000.0 / π
        root = add_segment!(tree, (0.0, 0.0, 0.0), (root_len, 0.0, 0.0), 1.0, Int32(-1))
        child1 = add_segment!(tree, (root_len, 0.0, 0.0), (root_len + child_len, 1.0, 0.0), 1.0, root)
        child2 = add_segment!(tree, (root_len, 0.0, 0.0), (root_len + child_len, -1.0, 0.0), 1.0, root)

        tree.segments.flow[root] = 2e-6
        tree.segments.flow[child1] = 1e-6
        tree.segments.flow[child2] = 1e-6

        times = collect(0.0:0.1:3.0)
        result = simulate_contrast_transport(
            tree,
            times;
            root_input=fill(8.0, length(times)),
            recompute_hemodynamics=false,
            storage_type=Float64,
        )

        @test result.transit_time_s[root] ≈ 0.25 atol = 1e-6
        @test result.transit_time_s[child1] ≈ 1.0 atol = 1e-4
        @test result.transit_time_s[child2] ≈ 1.0 atol = 1e-4
        @test result.concentration[child1, :] ≈ result.concentration[child2, :]
        @test result.concentration[child1, findfirst(==(0.5), times)] < result.concentration[root, findfirst(==(0.5), times)]
    end

    @testset "forest wrapper returns one result per tree" begin
        tree = VascularTree("LAD", 4)
        length_mm = 1000.0 / π
        seg_id = add_segment!(tree, (0.0, 0.0, 0.0), (length_mm, 0.0, 0.0), 1.0, Int32(-1))
        tree.segments.flow[seg_id] = 1e-6

        tree2 = deepcopy(tree)
        tree2.name = "RCA"

        cfg = TreeConfig("LAD", (0.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 1, 0.5)
        cfg2 = TreeConfig("RCA", (1.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 1, 0.5)
        domain = BoxDomain((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))
        tmap = initialize_territories(domain, [cfg, cfg2])
        forest = CoronaryForest(Dict("LAD" => tree, "RCA" => tree2), tmap, kassab_lad_params())

        times = collect(0.0:0.2:1.0)
        results = simulate_forest_contrast(
            forest,
            times;
            root_inputs=Dict("LAD" => 5.0, "RCA" => 7.0),
            recompute_hemodynamics=false,
            storage_type=Float64,
        )

        @test haskey(results, "LAD")
        @test haskey(results, "RCA")
        @test size(results["LAD"].concentration) == (1, length(times))
        @test size(results["RCA"].concentration) == (1, length(times))

        threaded_results = simulate_forest_contrast(
            forest,
            times;
            root_inputs=Dict("LAD" => 5.0, "RCA" => 7.0),
            recompute_hemodynamics=false,
            storage_type=Float64,
            threaded=true,
        )

        @test threaded_results["LAD"].concentration ≈ results["LAD"].concentration
        @test threaded_results["RCA"].concentration ≈ results["RCA"].concentration
    end
end
