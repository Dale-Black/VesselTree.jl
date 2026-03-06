using Test
using VesselTree
import AcceleratedKernels as AK

@testset "Distance Kernels" begin

    @testset "point_segment_distance — point on segment" begin
        # Point is the midpoint of segment (0,0,0)→(2,0,0)
        d = point_segment_distance(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0)
        @test d ≈ 0.0 atol = 1e-12
    end

    @testset "point_segment_distance — point at proximal" begin
        d = point_segment_distance(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
        @test d ≈ 0.0 atol = 1e-12
    end

    @testset "point_segment_distance — point at distal" begin
        d = point_segment_distance(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
        @test d ≈ 0.0 atol = 1e-12
    end

    @testset "point_segment_distance — perpendicular from midpoint" begin
        # Segment along x-axis, point offset in y
        d = point_segment_distance(0.5, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
        @test d ≈ 3.0 atol = 1e-12
    end

    @testset "point_segment_distance — closest to proximal endpoint" begin
        # Point is behind the proximal end
        d = point_segment_distance(-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0)
        @test d ≈ sqrt(2.0) atol = 1e-12
    end

    @testset "point_segment_distance — closest to distal endpoint" begin
        # Point is beyond the distal end
        d = point_segment_distance(3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0)
        @test d ≈ sqrt(2.0) atol = 1e-12
    end

    @testset "point_segment_distance — 3D diagonal segment" begin
        # Segment from (0,0,0) to (1,1,1), point at (0,0,0)
        d = point_segment_distance(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        @test d ≈ 0.0 atol = 1e-12
    end

    @testset "point_segment_distance — zero-length segment" begin
        # Degenerate segment: both endpoints at (1,1,1)
        d = point_segment_distance(4.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        @test d ≈ 5.0 atol = 1e-12  # sqrt(9+16+0)
    end

    @testset "point_segment_distance_sq — consistency" begin
        cx, cy, cz = 2.5, 3.0, -1.0
        px, py, pz = 0.0, 0.0, 0.0
        dx, dy, dz = 4.0, 0.0, 0.0
        d = point_segment_distance(cx, cy, cz, px, py, pz, dx, dy, dz)
        d_sq = point_segment_distance_sq(cx, cy, cz, px, py, pz, dx, dy, dz)
        @test d_sq ≈ d * d atol = 1e-12
    end

    @testset "compute_all_distances! — basic" begin
        tree = VascularTree("test", 100)
        # Three segments along x-axis
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5, Int32(-1))
        add_segment!(tree, (2.0, 0.0, 0.0), (3.0, 0.0, 0.0), 0.5, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (6.0, 0.0, 0.0), 0.5, Int32(-1))

        distances = zeros(100)
        compute_all_distances!(distances, tree.segments, 1.5, 0.0, 0.0, 3)

        @test distances[1] ≈ 0.5 atol = 1e-10   # closest to seg 1 distal
        @test distances[2] ≈ 0.5 atol = 1e-10   # closest to seg 2 proximal
        @test distances[3] ≈ 3.5 atol = 1e-10   # closest to seg 3 proximal
    end

    @testset "compute_all_distances! — matches scalar loop" begin
        tree = VascularTree("test", 1100)
        rng = Random.MersenneTwister(42)
        for i in 1:1000
            px, py, pz = rand(rng, 3) .* 10.0
            dx, dy, dz = px + rand(rng) * 2.0, py + rand(rng) * 2.0, pz + rand(rng) * 2.0
            add_segment!(tree, (px, py, pz), (dx, dy, dz), 0.1, Int32(-1))
        end

        cx, cy, cz = 5.0, 5.0, 5.0
        n = tree.segments.n

        # AK version
        distances_ak = zeros(1100)
        compute_all_distances!(distances_ak, tree.segments, cx, cy, cz, n)

        # Scalar loop version
        distances_scalar = zeros(n)
        for i in 1:n
            distances_scalar[i] = point_segment_distance(
                cx, cy, cz,
                tree.segments.proximal_x[i], tree.segments.proximal_y[i], tree.segments.proximal_z[i],
                tree.segments.distal_x[i], tree.segments.distal_y[i], tree.segments.distal_z[i],
            )
        end

        for i in 1:n
            @test distances_ak[i] ≈ distances_scalar[i] atol = 1e-12
        end
    end

    @testset "find_nearest_segments — correct indices" begin
        distances = zeros(100)
        # Set up known distances
        distances[1] = 10.0
        distances[2] = 5.0
        distances[3] = 1.0
        distances[4] = 8.0
        distances[5] = 3.0

        nearest = find_nearest_segments(distances, 5, 3)
        @test length(nearest) == 3
        @test nearest[1] == 3   # distance 1.0
        @test nearest[2] == 5   # distance 3.0
        @test nearest[3] == 2   # distance 5.0
    end

    @testset "find_nearest_segments — k > n" begin
        distances = [2.0, 1.0, 3.0, 0.0, 0.0]
        nearest = find_nearest_segments(distances, 3, 10)
        @test length(nearest) == 3
    end

    @testset "find_nearest_segments — 1000 segments" begin
        rng = Random.MersenneTwister(99)
        tree = VascularTree("test", 1100)
        for i in 1:1000
            px, py, pz = rand(rng, 3) .* 100.0
            dx, dy, dz = px + rand(rng), py + rand(rng), pz + rand(rng)
            add_segment!(tree, (px, py, pz), (dx, dy, dz), 0.1, Int32(-1))
        end

        distances = zeros(1100)
        compute_all_distances!(distances, tree.segments, 50.0, 50.0, 50.0, 1000)

        nearest = find_nearest_segments(distances, 1000, 5)
        @test length(nearest) == 5

        # Nearest should have the smallest distance
        for i in 2:5
            @test distances[nearest[i]] >= distances[nearest[i-1]]
        end
    end

end
