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

end
