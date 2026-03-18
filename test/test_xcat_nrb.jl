using VesselTree
using Test
using StaticArrays
using LinearAlgebra

@testset "XCAT Centerline Import" begin
    centers = [SVector(Float64(i), 0.0, 0.0) for i in 0:10]
    radii = fill(1.0, length(centers))
    line = XCATCenterline("synthetic", centers, radii, :v)

    @test xcat_centerline_length_mm(line) ≈ 10.0 atol = 1e-9

    resampled = xcat_resample_centerline(
        line;
        min_spacing_mm=2.0,
        max_spacing_mm=2.0,
        radius_spacing_factor=1.0,
    )

    @test first(resampled.centers) == first(line.centers)
    @test last(resampled.centers) == last(line.centers)
    @test length(resampled.centers) < length(line.centers)

    spacings = [norm(resampled.centers[i + 1] - resampled.centers[i]) for i in 1:(length(resampled.centers) - 1)]
    @test all(s -> s > 0.9, spacings)
    @test maximum(spacings) <= 4.1

    short = XCATCenterline(
        "short",
        [SVector(0.0, 0.0, 0.0), SVector(0.2, 0.0, 0.0)],
        [1.0, 1.0],
        :v,
    )
    parent = XCATCenterline(
        "parent",
        [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(2.0, 0.0, 0.0), SVector(3.0, 0.0, 0.0)],
        fill(1.0, 4),
        :v,
    )
    child = XCATCenterline(
        "child",
        [SVector(2.1, 0.1, 0.0), SVector(4.0, 0.1, 0.0), SVector(5.0, 0.1, 0.0)],
        fill(1.0, 3),
        :v,
    )
    snapped_child, conn = VesselTree._make_tree_connection(parent, child)
    @test conn.parent_index == length(parent.centers)
    @test conn.child_index == 1
    @test first(snapped_child.centers) == child.centers[1]
end
