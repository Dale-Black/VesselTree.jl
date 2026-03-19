using Test
using JLD2
using StaticArrays
using VesselTree

@testset "Persistence and export" begin
    tree = VascularTree("LAD", 8)
    seg_id = add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 0.5, Int32(-1))
    tree.segments.flow[seg_id] = 1e-9

    cfg = TreeConfig("LAD", (0.0, 0.0, 0.0), 0.5, (1.0, 0.0, 0.0), 1, 1.0)
    domain = BoxDomain((-5.0, -5.0, -5.0), (15.0, 5.0, 5.0))
    territory = initialize_territories(domain, [cfg])
    forest = CoronaryForest(Dict("LAD" => tree), territory, kassab_lad_params())

    times = collect(0.0:0.5:2.0)
    results = simulate_forest_contrast(
        forest,
        times;
        root_inputs=Dict("LAD" => fill(5.0, length(times))),
        recompute_hemodynamics=false,
        storage_type=Float32,
    )

    tmp = mktempdir()
    ts = "20260318T120000"

    tree_snapshot = save_forest_snapshot(tmp, forest; timestamp=ts, formats=[:jld2, :csv, :graph_json, :wenbo_txt])
    @test isfile(tree_snapshot.manifest_path)
    @test isfile(tree_snapshot.tree_paths["LAD"]["jld2"])
    @test endswith(tree_snapshot.tree_paths["LAD"]["jld2"], "LAD-$(ts).jld2")
    @test isfile(tree_snapshot.tree_paths["LAD"]["csv"])
    @test isfile(tree_snapshot.tree_paths["LAD"]["graph_json"])
    @test isfile(tree_snapshot.tree_paths["LAD"]["wenbo_txt"])

    contrast_snapshot = save_contrast_snapshot(
        tmp,
        forest,
        results;
        timestamp=ts,
        voxel_spacing_mm=(2.0, 2.0, 2.0),
        voxel_supersample=1,
        dense_mode=:always,
    )
    @test isfile(contrast_snapshot.segment_path)
    @test isfile(contrast_snapshot.volume_path)
    contrast_file = JLD2.load(contrast_snapshot.volume_path)
    @test haskey(contrast_file, "concentration_4d_mg_mL")
    @test size(contrast_file["concentration_4d_mg_mL"], 4) == length(times)
    @test size(contrast_file["voxel_indices_ijk"], 2) == 3

    surface = XCATNurbsSurface(
        "base_surface",
        2,
        2,
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        reshape(SVector{3, Float64}[
            SVector(0.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(1.0, 1.0, 0.0),
        ], 2, 2),
    )
    base_nrb = joinpath(tmp, "base.nrb")
    write_xcat_nrb(base_nrb, [surface])

    fused_nrb = joinpath(tmp, "fused.nrb")
    export_fused_nrb_model(base_nrb, fused_nrb, forest; name_prefix="grown_test")
    parsed = parse_xcat_nrb(fused_nrb)
    @test length(parsed) == 1 + tree.segments.n
    @test any(startswith(s.name, "grown_test_LAD_") for s in parsed)

    bundle = save_xcat_run_artifacts(
        joinpath(tmp, "bundle"),
        base_nrb,
        forest,
        results;
        timestamp=ts,
        voxel_spacing_mm=(2.0, 2.0, 2.0),
        voxel_supersample=1,
        dense_mode=:never,
    )
    @test isdir(bundle.run_dir)
    @test isfile(bundle.fused_nrb_path)
    @test isfile(bundle.manifest_path)
end
