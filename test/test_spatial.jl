using Test
using VesselTree
using Random

@testset "Spatial Indexing" begin

    @testset "build_grid — basic" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5, Int32(-1))
        add_segment!(tree, (5.0, 5.0, 5.0), (6.0, 5.0, 5.0), 0.5, Int32(-1))
        add_segment!(tree, (9.0, 9.0, 9.0), (10.0, 9.0, 9.0), 0.5, Int32(-1))

        domain = BoxDomain((-1.0, -1.0, -1.0), (11.0, 11.0, 11.0))
        grid = build_grid(tree.segments, 3, domain, 3.0)

        @test grid.cell_size == 3.0
        @test grid.dims[1] >= 1
        @test grid.dims[2] >= 1
        @test grid.dims[3] >= 1
    end

    @testset "query_nearby — finds correct segments" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (1.0, 1.0, 1.0), (2.0, 1.0, 1.0), 0.5, Int32(-1))  # midpoint ~(1.5,1,1)
        add_segment!(tree, (5.0, 5.0, 5.0), (6.0, 5.0, 5.0), 0.5, Int32(-1))  # midpoint ~(5.5,5,5)
        add_segment!(tree, (9.0, 9.0, 9.0), (10.0, 9.0, 9.0), 0.5, Int32(-1)) # midpoint ~(9.5,9,9)

        domain = BoxDomain((0.0, 0.0, 0.0), (11.0, 11.0, 11.0))
        grid = build_grid(tree.segments, 3, domain, 3.0)

        # Query near segment 1
        nearby = query_nearby(grid, 1.5, 1.0, 1.0, 2.0)
        @test 1 in nearby
        @test !(3 in nearby)  # segment 3 is far away

        # Query near segment 3
        nearby3 = query_nearby(grid, 9.5, 9.0, 9.0, 2.0)
        @test 3 in nearby3
    end

    @testset "query_nearby — no false negatives" begin
        rng = MersenneTwister(42)
        tree = VascularTree("test", 600)
        for _ in 1:500
            px, py, pz = rand(rng, 3) .* 10.0
            dx, dy, dz = px + rand(rng), py + rand(rng), pz + rand(rng)
            add_segment!(tree, (px, py, pz), (dx, dy, dz), 0.1, Int32(-1))
        end

        domain = BoxDomain((-1.0, -1.0, -1.0), (12.0, 12.0, 12.0))
        grid = build_grid(tree.segments, 500, domain, 2.0)

        # For several query points, verify nearest segment is in query result
        for _ in 1:20
            qx, qy, qz = rand(rng, 3) .* 10.0

            # Find actual nearest segment via brute force
            min_dist = Inf
            nearest_idx = 0
            for i in 1:500
                mx = (tree.segments.proximal_x[i] + tree.segments.distal_x[i]) / 2.0
                my = (tree.segments.proximal_y[i] + tree.segments.distal_y[i]) / 2.0
                mz = (tree.segments.proximal_z[i] + tree.segments.distal_z[i]) / 2.0
                d = sqrt((mx - qx)^2 + (my - qy)^2 + (mz - qz)^2)
                if d < min_dist
                    min_dist = d
                    nearest_idx = i
                end
            end

            # Query with sufficient radius
            nearby = query_nearby(grid, qx, qy, qz, min_dist + 3.0)  # +3 for cell size margin
            @test nearest_idx in nearby
        end
    end

    @testset "insert! — add single segment" begin
        domain = BoxDomain((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))
        tree = VascularTree("test", 10)
        add_segment!(tree, (1.0, 1.0, 1.0), (2.0, 1.0, 1.0), 0.5, Int32(-1))

        grid = build_grid(tree.segments, 1, domain, 5.0)

        # Insert a new segment
        VesselTree.insert!(grid, 2, 8.0, 8.0, 8.0)

        nearby = query_nearby(grid, 8.0, 8.0, 8.0, 1.0)
        @test 2 in nearby
    end

    @testset "SpatialGrid with SphereDomain" begin
        domain = SphereDomain((5.0, 5.0, 5.0), 5.0)
        tree = VascularTree("test", 100)
        add_segment!(tree, (5.0, 5.0, 5.0), (7.0, 5.0, 5.0), 0.5, Int32(-1))

        grid = build_grid(tree.segments, 1, domain, 2.0)
        @test grid.dims[1] >= 1
        nearby = query_nearby(grid, 6.0, 5.0, 5.0, 2.0)
        @test 1 in nearby
    end

    @testset "SpatialGrid with EllipsoidDomain" begin
        domain = EllipsoidDomain((0.0, 0.0, 0.0), (10.0, 5.0, 5.0))
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), 0.5, Int32(-1))

        grid = build_grid(tree.segments, 1, domain, 3.0)
        @test grid.dims[1] >= 1
        nearby = query_nearby(grid, 1.0, 0.0, 0.0, 3.0)
        @test 1 in nearby
    end

end
