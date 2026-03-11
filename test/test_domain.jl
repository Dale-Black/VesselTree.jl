using Test
using VesselTree
using Random
using DelimitedFiles

@testset "Domain Types" begin

    @testset "SphereDomain — in_domain" begin
        d = SphereDomain((0.0, 0.0, 0.0), 1.0)
        @test in_domain(d, (0.0, 0.0, 0.0))       # center
        @test in_domain(d, (1.0, 0.0, 0.0))        # on surface
        @test in_domain(d, (0.5, 0.5, 0.0))         # inside
        @test !in_domain(d, (1.1, 0.0, 0.0))        # outside
        @test !in_domain(d, (0.8, 0.8, 0.0))        # outside (0.64+0.64 > 1)
    end

    @testset "SphereDomain — signed_distance" begin
        d = SphereDomain((0.0, 0.0, 0.0), 2.0)
        @test signed_distance(d, (0.0, 0.0, 0.0)) ≈ -2.0   # center
        @test signed_distance(d, (2.0, 0.0, 0.0)) ≈ 0.0     # surface
        @test signed_distance(d, (3.0, 0.0, 0.0)) ≈ 1.0     # outside
        @test signed_distance(d, (1.0, 0.0, 0.0)) ≈ -1.0    # inside
    end

    @testset "SphereDomain — sample_point" begin
        d = SphereDomain((1.0, 2.0, 3.0), 5.0)
        rng = MersenneTwister(42)
        for _ in 1:100
            p = sample_point(d, rng)
            @test in_domain(d, p)
        end
    end

    @testset "BoxDomain — in_domain" begin
        d = BoxDomain((0.0, 0.0, 0.0), (1.0, 2.0, 3.0))
        @test in_domain(d, (0.5, 1.0, 1.5))         # inside
        @test in_domain(d, (0.0, 0.0, 0.0))         # corner
        @test in_domain(d, (1.0, 2.0, 3.0))         # opposite corner
        @test !in_domain(d, (-0.1, 0.5, 0.5))       # outside x
        @test !in_domain(d, (0.5, 2.1, 0.5))        # outside y
        @test !in_domain(d, (0.5, 0.5, 3.1))        # outside z
    end

    @testset "BoxDomain — signed_distance" begin
        d = BoxDomain((0.0, 0.0, 0.0), (2.0, 2.0, 2.0))
        @test signed_distance(d, (1.0, 1.0, 1.0)) ≈ -1.0   # center
        @test signed_distance(d, (2.0, 1.0, 1.0)) ≈ 0.0     # on face
        @test signed_distance(d, (3.0, 1.0, 1.0)) ≈ 1.0     # outside
        @test signed_distance(d, (0.5, 0.5, 0.5)) ≈ -0.5    # inside near corner
    end

    @testset "BoxDomain — sample_point" begin
        d = BoxDomain((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
        rng = MersenneTwister(123)
        for _ in 1:100
            p = sample_point(d, rng)
            @test in_domain(d, p)
        end
    end

    @testset "EllipsoidDomain — in_domain" begin
        d = EllipsoidDomain((0.0, 0.0, 0.0), (2.0, 1.0, 1.0))
        @test in_domain(d, (0.0, 0.0, 0.0))         # center
        @test in_domain(d, (2.0, 0.0, 0.0))         # on surface along major axis
        @test in_domain(d, (0.0, 1.0, 0.0))         # on surface along minor axis
        @test !in_domain(d, (2.1, 0.0, 0.0))        # outside major
        @test !in_domain(d, (0.0, 1.1, 0.0))        # outside minor
        @test in_domain(d, (1.0, 0.5, 0.0))         # inside
    end

    @testset "EllipsoidDomain — signed_distance" begin
        d = EllipsoidDomain((0.0, 0.0, 0.0), (2.0, 2.0, 2.0))
        # Sphere-like ellipsoid — approximate SDF
        @test signed_distance(d, (0.0, 0.0, 0.0)) < 0.0   # inside
        @test signed_distance(d, (2.0, 0.0, 0.0)) ≈ 0.0 atol = 0.01  # surface
        @test signed_distance(d, (4.0, 0.0, 0.0)) > 0.0   # outside
    end

    @testset "EllipsoidDomain — sample_point" begin
        d = EllipsoidDomain((5.0, 5.0, 5.0), (3.0, 2.0, 1.0))
        rng = MersenneTwister(99)
        for _ in 1:100
            p = sample_point(d, rng)
            @test in_domain(d, p)
        end
    end

    @testset "SphereDomain — offset center" begin
        d = SphereDomain((10.0, 20.0, 30.0), 5.0)
        @test in_domain(d, (10.0, 20.0, 30.0))
        @test in_domain(d, (15.0, 20.0, 30.0))
        @test !in_domain(d, (15.1, 20.0, 30.0))
    end

    @testset "CSVVolumeDomain — v3-compatible interface" begin
        mktempdir() do tmp
            points = [
                1.0 0.0 0.0
                -1.0 0.0 0.0
                0.0 1.0 0.0
                0.0 -1.0 0.0
                0.0 0.0 1.0
                0.0 0.0 -1.0
            ]
            normals = copy(points)

            points_csv = joinpath(tmp, "points.csv")
            normals_csv = joinpath(tmp, "normals.csv")
            writedlm(points_csv, points, ',')
            writedlm(normals_csv, normals, ',')

            d = csv_volume_domain(
                points_csv,
                normals_csv;
                rng=MersenneTwister(7),
                target_interior=512,
                max_candidates=8192,
                batch_size=512,
                coordinate_scale=1.0,
            )

            @test d isa CSVVolumeDomain
            @test d.length_scale ≈ 1.0
            @test in_domain(d, (0.0, 0.0, 0.0))
            @test !in_domain(d, (1.5, 0.0, 0.0))

            p = sample_point(d, MersenneTwister(9))
            @test in_domain(d, p)

            proj = project_to_domain(d, (1.5, 0.0, 0.0))
            @test in_domain(d, proj)
            @test signed_distance(d, proj) <= 1e-4
        end
    end

end
