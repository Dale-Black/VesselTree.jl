using Test
using VesselTree
using Random

@testset "Export & Persistence" begin

    params = kassab_coronary_params()

    # --- JLD2 save/load roundtrip ---

    @testset "save_tree / load_tree — roundtrip" begin
        tree = VascularTree("test_roundtrip", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.3, Int32(1))

        tmpfile = tempname() * ".jld2"
        try
            save_tree(tmpfile, tree)
            @test isfile(tmpfile)

            loaded = load_tree(tmpfile)
            @test loaded.name == tree.name
            @test loaded.segments.n == tree.segments.n
            @test loaded.n_bifurcations == tree.n_bifurcations

            for i in 1:tree.segments.n
                @test loaded.segments.radius[i] ≈ tree.segments.radius[i] rtol = 1e-10
                @test loaded.segments.proximal_x[i] ≈ tree.segments.proximal_x[i] rtol = 1e-10
                @test loaded.segments.distal_x[i] ≈ tree.segments.distal_x[i] rtol = 1e-10
            end
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    @testset "save_tree / load_tree — CCO tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("cco_save", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 20, params; rng=rng, kassab=true)

        tmpfile = tempname() * ".jld2"
        try
            save_tree(tmpfile, tree)
            loaded = load_tree(tmpfile)
            @test loaded.segments.n == tree.segments.n
            @test loaded.n_terminals == tree.n_terminals
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    # --- VTP export ---

    @testset "export_centerlines_vtp — creates file" begin
        tree = VascularTree("vtp_test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.3, Int32(1))

        tmpfile = tempname()
        try
            outpath = export_centerlines_vtp(tree, tmpfile)
            @test isfile(outpath)
            @test endswith(outpath, ".vtp")
        finally
            isfile(tmpfile * ".vtp") && rm(tmpfile * ".vtp")
        end
    end

    @testset "export_centerlines_vtp — CCO tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("cco_vtp", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 20, params; rng=rng, kassab=true)

        tmpfile = tempname()
        try
            outpath = export_centerlines_vtp(tree, tmpfile)
            @test isfile(outpath)
            # File should have some content
            @test filesize(outpath) > 100
        finally
            isfile(tmpfile * ".vtp") && rm(tmpfile * ".vtp")
        end
    end

    # --- Forest export ---

    @testset "export_forest_vtp — creates per-tree files" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 15.0)
        configs = [
            TreeConfig("A", (-4.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 5, 0.5),
            TreeConfig("B", (4.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 5, 0.5),
        ]
        rng = MersenneTwister(42)
        forest = generate_coronary_forest(domain, params; tree_configs=configs, rng=rng)

        tmpdir = mktempdir()
        try
            export_forest_vtp(forest, tmpdir)
            files = readdir(tmpdir)
            @test length(files) >= 2
            @test any(f -> endswith(f, ".vtp"), files)
        finally
            rm(tmpdir; recursive=true)
        end
    end

    # --- STL mesh export ---

    @testset "export_stl — creates valid binary STL" begin
        tree = VascularTree("stl_test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.3, Int32(1))

        tmpfile = tempname() * ".stl"
        try
            outpath = export_stl(tree, tmpfile)
            @test isfile(outpath)
            # Binary STL: 80 header + 4 bytes triangle count + 50 bytes per triangle
            n_segs = 3
            nf = 16  # default circumferential_resolution
            expected_triangles = 2 * nf * n_segs
            expected_size = 80 + 4 + 50 * expected_triangles
            @test filesize(outpath) == expected_size

            # Verify triangle count in header
            data = read(outpath)
            tri_count = reinterpret(UInt32, data[81:84])[1]
            @test tri_count == expected_triangles
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    @testset "export_stl — custom resolution" begin
        tree = VascularTree("stl_res", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.5, Int32(-1))

        tmpfile = tempname() * ".stl"
        try
            export_stl(tree, tmpfile; circumferential_resolution=8)
            data = read(tmpfile)
            tri_count = reinterpret(UInt32, data[81:84])[1]
            @test tri_count == 2 * 8 * 1  # 2 tris per facet * 8 facets * 1 segment
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    # --- Graph JSON export ---

    @testset "export_graph_json — valid JSON structure" begin
        tree = VascularTree("json_test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))

        tmpfile = tempname() * ".json"
        try
            export_graph_json(tree, tmpfile)
            @test isfile(tmpfile)
            content = read(tmpfile, String)
            # Basic JSON structure checks
            @test startswith(content, "{")
            @test endswith(content, "}")
            @test occursin("\"nodes\"", content)
            @test occursin("\"edges\"", content)
            @test occursin("\"name\"", content)
            @test occursin("json_test", content)
            @test occursin("\"n_segments\"", content)
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    # --- CSV export ---

    @testset "export_csv — valid CSV with headers" begin
        tree = VascularTree("csv_test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (10.0, 0.0, 0.0), (15.0, -1.0, 0.0), 0.3, Int32(1))

        tmpfile = tempname() * ".csv"
        try
            export_csv(tree, tmpfile)
            @test isfile(tmpfile)
            lines = readlines(tmpfile)
            @test length(lines) == 4  # header + 3 segments
            # Check header
            @test startswith(lines[1], "id,proximal_x")
            @test occursin("radius", lines[1])
            @test occursin("strahler_order", lines[1])
            # Check data rows have correct number of columns
            header_cols = length(split(lines[1], ","))
            for row in lines[2:end]
                @test length(split(row, ",")) == header_cols
            end
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

end
