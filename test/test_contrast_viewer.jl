using Test
using VesselTree

@testset "Contrast viewer" begin
    tree = VascularTree("LAD", 4)
    length_mm = 1000.0 / π
    seg_id = add_segment!(tree, (0.0, 0.0, 0.0), (length_mm, 0.0, 0.0), 1.0, Int32(-1))
    tree.segments.flow[seg_id] = 1e-6

    cfg = TreeConfig("LAD", (0.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 1, 1.0)
    domain = BoxDomain((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))
    tmap = initialize_territories(domain, [cfg])
    forest = CoronaryForest(Dict("LAD" => tree), tmap, kassab_lad_params())

    times = collect(0.0:0.2:1.0)
    results = simulate_forest_contrast(
        forest,
        times;
        root_inputs=Dict("LAD" => fill(5.0, length(times))),
        recompute_hemodynamics=false,
        storage_type=Float64,
    )

    viewer = prepare_contrast_viewer_data(
        forest,
        results;
        min_radius_um=0.0,
        max_segments_per_tree=10,
        time_stride=2,
    )

    @test length(viewer.trees) == 1
    @test viewer.trees[1].name == "LAD"
    @test length(viewer.times) == 3
    @test size(viewer.trees[1].concentration_by_time) == (3, 1)
    @test length(viewer.trees[1].segment_length_mm) == 1
    @test length(viewer.trees[1].segment_diameter_um) == 1
    @test viewer.cmax > 0

    mktempdir() do dir
        html_path = export_contrast_viewer_html(joinpath(dir, "index.html"), viewer; title="Test Contrast Viewer")
        @test isfile(html_path)
        html = read(html_path, String)
        @test occursin("Plotly.newPlot", html)
        @test occursin("Test Contrast Viewer", html)
        @test occursin("Dynamic segment-average iodine concentration", html)
        @test occursin("customdata[1]", html)
        @test occursin("customdata[2]", html)
    end
end
