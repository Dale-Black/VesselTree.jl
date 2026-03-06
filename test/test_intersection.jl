using Test
using VesselTree
import AcceleratedKernels as AK

@testset "Intersection Testing" begin

    @testset "segments_intersect — crossing X pattern" begin
        # Two segments crossing in 2D (z=0 plane)
        # Segment A: (-1,0,0) → (1,0,0)
        # Segment B: (0,-1,0) → (0,1,0)
        d_sq = VesselTree.segments_min_distance_sq(
            -1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, -1.0, 0.0, 0.0, 1.0, 0.0,
        )
        @test d_sq ≈ 0.0 atol = 1e-12
        @test segments_intersect(
            -1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.1,
        )
    end

    @testset "segments_intersect — parallel non-overlapping" begin
        # Two parallel segments separated by 2 units in y
        d_sq = VesselTree.segments_min_distance_sq(
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 2.0, 0.0, 1.0, 2.0, 0.0,
        )
        @test d_sq ≈ 4.0 atol = 1e-12
        @test !segments_intersect(
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0,
        )
    end

    @testset "segments_intersect — skew segments in 3D" begin
        # Segment A along x-axis, segment B along z-axis offset by 1 in y
        d_sq = VesselTree.segments_min_distance_sq(
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
            1.0, 1.0, -1.0, 1.0, 1.0, 1.0,
        )
        @test d_sq ≈ 1.0 atol = 1e-12  # separated by 1 in y
    end

    @testset "segments_intersect — shared endpoint" begin
        # Two segments sharing an endpoint
        d_sq = VesselTree.segments_min_distance_sq(
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 2.0, 1.0, 0.0,
        )
        @test d_sq ≈ 0.0 atol = 1e-12
    end

    @testset "segments_intersect — collinear overlapping" begin
        # Two collinear segments that overlap
        d_sq = VesselTree.segments_min_distance_sq(
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 3.0, 0.0, 0.0,
        )
        @test d_sq ≈ 0.0 atol = 1e-12
    end

    @testset "segments_intersect — collinear non-overlapping" begin
        # Two collinear segments with a gap
        d_sq = VesselTree.segments_min_distance_sq(
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            2.0, 0.0, 0.0, 3.0, 0.0, 0.0,
        )
        @test d_sq ≈ 1.0 atol = 1e-12
    end

    @testset "segments_intersect — zero-length segments" begin
        # Both are points
        d_sq = VesselTree.segments_min_distance_sq(
            1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            4.0, 0.0, 0.0, 4.0, 0.0, 0.0,
        )
        @test d_sq ≈ 9.0 atol = 1e-12
    end

    @testset "check_intersections! — AK kernel" begin
        tree = VascularTree("test", 100)
        # Horizontal segments at different y positions
        add_segment!(tree, (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), 0.1, Int32(-1))  # y=0
        add_segment!(tree, (0.0, 1.0, 0.0), (2.0, 1.0, 0.0), 0.1, Int32(-1))  # y=1
        add_segment!(tree, (0.0, 3.0, 0.0), (2.0, 3.0, 0.0), 0.1, Int32(-1))  # y=3

        results = Vector{Bool}(undef, 100)
        # New vertical segment at x=1, from y=-1 to y=2 — crosses segs 1 and 2, not 3
        check_intersections!(results, tree.segments,
            1.0, -1.0, 0.0, 1.0, 2.0, 0.0, 0.1, 3)

        @test results[1] == true   # crosses y=0
        @test results[2] == true   # crosses y=1
        @test results[3] == false  # doesn't reach y=3
    end

    @testset "has_any_intersection — AK reduction" begin
        results = Vector{Bool}([false, false, false, true, false])
        @test has_any_intersection(results, 5) == true

        results_none = Vector{Bool}([false, false, false, false, false])
        @test has_any_intersection(results_none, 5) == false
    end

    @testset "check_intersections! — AK vs scalar on 1000 segments" begin
        rng = Random.MersenneTwister(42)
        tree = VascularTree("test", 1100)
        for _ in 1:1000
            px, py, pz = rand(rng, 3) .* 10.0
            dx, dy, dz = px + rand(rng), py + rand(rng), pz + rand(rng)
            add_segment!(tree, (px, py, pz), (dx, dy, dz), 0.1, Int32(-1))
        end

        npx, npy, npz = 5.0, 5.0, 5.0
        ndx, ndy, ndz = 6.0, 6.0, 6.0
        min_dist = 0.5
        n = tree.segments.n

        # AK version
        results_ak = Vector{Bool}(undef, 1100)
        check_intersections!(results_ak, tree.segments, npx, npy, npz, ndx, ndy, ndz, min_dist, n)

        # Scalar version
        results_scalar = Vector{Bool}(undef, n)
        for i in 1:n
            results_scalar[i] = segments_intersect(
                tree.segments.proximal_x[i], tree.segments.proximal_y[i], tree.segments.proximal_z[i],
                tree.segments.distal_x[i], tree.segments.distal_y[i], tree.segments.distal_z[i],
                npx, npy, npz, ndx, ndy, ndz, min_dist,
            )
        end

        for i in 1:n
            @test results_ak[i] == results_scalar[i]
        end
    end

    @testset "check_domain_crossing — inside sphere" begin
        d = SphereDomain((0.0, 0.0, 0.0), 10.0)
        # Segment fully inside
        @test !check_domain_crossing(d, (1.0, 1.0, 1.0), (2.0, 2.0, 2.0))
    end

    @testset "check_domain_crossing — crossing sphere boundary" begin
        d = SphereDomain((0.0, 0.0, 0.0), 5.0)
        # Segment from inside to outside
        @test check_domain_crossing(d, (0.0, 0.0, 0.0), (10.0, 0.0, 0.0))
    end

    @testset "check_domain_crossing — ellipsoid" begin
        d = EllipsoidDomain((0.0, 0.0, 0.0), (5.0, 3.0, 3.0))
        # Fully inside
        @test !check_domain_crossing(d, (1.0, 0.0, 0.0), (2.0, 0.0, 0.0))
        # Crossing boundary in y direction
        @test check_domain_crossing(d, (0.0, 0.0, 0.0), (0.0, 10.0, 0.0))
    end

end
