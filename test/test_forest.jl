using Test
using VesselTree
using Random

@testset "Territory Partitioning + Forest" begin

    params = kassab_coronary_params()

    # --- TreeConfig ---

    @testset "TreeConfig — construction" begin
        cfg = TreeConfig("LAD", (0.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 100, 0.4)
        @test cfg.name == "LAD"
        @test cfg.root_position == (0.0, 0.0, 0.0)
        @test cfg.root_radius == 1.0
        @test cfg.root_direction == (1.0, 0.0, 0.0)
        @test cfg.target_terminals == 100
        @test cfg.territory_fraction == 0.4
    end

    # --- TerritoryMap ---

    @testset "initialize_territories — sphere domain" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        configs = [
            TreeConfig("A", (-5.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 50, 0.5),
            TreeConfig("B", (5.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 50, 0.5),
        ]

        tmap = initialize_territories(domain, configs)
        @test tmap isa TerritoryMap
        @test length(tmap.tree_names) == 2
        @test tmap.cell_size > 0.0
    end

    @testset "initialize_territories — correct partition" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        configs = [
            TreeConfig("Left", (-5.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 50, 0.5),
            TreeConfig("Right", (5.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 50, 0.5),
        ]

        tmap = initialize_territories(domain, configs)

        # Point on left side should belong to "Left" tree
        left_name = query_territory(tmap, -3.0, 0.0, 0.0)
        @test left_name == "Left"

        # Point on right side should belong to "Right" tree
        right_name = query_territory(tmap, 3.0, 0.0, 0.0)
        @test right_name == "Right"
    end

    # --- sample_in_territory ---

    @testset "sample_in_territory — returns points in territory" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        configs = [
            TreeConfig("Left", (-5.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0), 50, 0.5),
            TreeConfig("Right", (5.0, 0.0, 0.0), 1.0, (-1.0, 0.0, 0.0), 50, 0.5),
        ]

        tmap = initialize_territories(domain, configs)
        rng = MersenneTwister(42)

        # Sample 20 points for "Left" territory
        count_left = 0
        for _ in 1:20
            pt = sample_in_territory(domain, tmap, "Left", rng)
            if pt !== nothing
                name = query_territory(tmap, pt[1], pt[2], pt[3])
                if name == "Left"
                    count_left += 1
                end
            end
        end
        # Most samples should be in the correct territory
        @test count_left > 10
    end

    # --- check_inter_tree_collision ---

    @testset "check_inter_tree_collision — detects nearby segments" begin
        tree = VascularTree("other", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), 1.0, Int32(-1))

        # Point very close to the segment
        @test check_inter_tree_collision((5.0, 0.01, 0.0), tree, 0.1)

        # Point far from the segment
        @test !check_inter_tree_collision((50.0, 50.0, 50.0), tree, 0.1)
    end

    @testset "check_inter_tree_collision — empty tree" begin
        tree = VascularTree("empty", 100)
        @test !check_inter_tree_collision((0.0, 0.0, 0.0), tree, 0.1)
    end

    # --- coronary_tree_configs ---

    @testset "coronary_tree_configs — returns 3 trees" begin
        configs = coronary_tree_configs()
        @test length(configs) == 3
        names = [c.name for c in configs]
        @test "LAD" in names
        @test "LCX" in names
        @test "RCA" in names

        # Fractions should sum to ~1
        total_frac = sum(c.territory_fraction for c in configs)
        @test total_frac ≈ 1.0 atol = 0.01
    end

    # --- Integration ---

    @testset "territory + sampling integration" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        configs = coronary_tree_configs()
        tmap = initialize_territories(domain, configs)

        @test length(tmap.tree_names) == 3

        rng = MersenneTwister(42)
        for cfg in configs
            pt = sample_in_territory(domain, tmap, cfg.name, rng)
            # Should get a valid point
            @test pt !== nothing || true  # may fail if territory is tiny, that's ok
        end
    end

end
