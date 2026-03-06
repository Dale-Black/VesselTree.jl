using Test
using VesselTree
import AcceleratedKernels as AK

@testset "Kamiya Optimization" begin

    gamma = VesselTree.MURRAY_GAMMA  # 7/3

    @testset "compute_radii_symmetric" begin
        r_parent = 1.0
        r1, r2 = compute_radii_symmetric(r_parent, gamma)
        @test r1 ≈ r2
        # Murray's law: r_parent^gamma = r1^gamma + r2^gamma
        @test r_parent^gamma ≈ r1^gamma + r2^gamma atol = 1e-10
    end

    @testset "compute_radii_asymmetric" begin
        r_parent = 1.0
        asymmetry = 0.5  # r_small = 0.5 * r_large
        r_large, r_small = compute_radii_asymmetric(r_parent, asymmetry, gamma)
        @test r_small ≈ asymmetry * r_large atol = 1e-10
        @test r_parent^gamma ≈ r_large^gamma + r_small^gamma atol = 1e-10

        # Extreme asymmetry
        r_l, r_s = compute_radii_asymmetric(1.0, 0.1, gamma)
        @test r_s < r_l
        @test 1.0^gamma ≈ r_l^gamma + r_s^gamma atol = 1e-10
    end

    @testset "surface_cost_at_t — endpoints" begin
        # Segment from (0,0,0) to (4,0,0), terminal at (2,3,0)
        seg_p = (0.0, 0.0, 0.0)
        seg_d = (4.0, 0.0, 0.0)
        r = 1.0

        cost_start = surface_cost_at_t(0.01, seg_p..., seg_d..., r, 2.0, 3.0, 0.0, gamma)
        cost_mid = surface_cost_at_t(0.5, seg_p..., seg_d..., r, 2.0, 3.0, 0.0, gamma)
        cost_end = surface_cost_at_t(0.99, seg_p..., seg_d..., r, 2.0, 3.0, 0.0, gamma)

        # All should be positive
        @test cost_start > 0.0
        @test cost_mid > 0.0
        @test cost_end > 0.0
    end

    @testset "surface_cost_at_t — optimizer beats endpoints" begin
        seg_p = (0.0, 0.0, 0.0)
        seg_d = (10.0, 0.0, 0.0)
        r = 1.0

        t_opt, cost_opt = optimize_bifurcation_point(
            seg_p..., seg_d..., r, 5.0, 5.0, 0.0, gamma,
        )
        cost_at_01 = surface_cost_at_t(0.01, seg_p..., seg_d..., r, 5.0, 5.0, 0.0, gamma)
        cost_at_99 = surface_cost_at_t(0.99, seg_p..., seg_d..., r, 5.0, 5.0, 0.0, gamma)

        # Optimized cost should be <= endpoint costs
        @test cost_opt <= cost_at_01
        @test cost_opt <= cost_at_99
    end

    @testset "optimize_bifurcation_point — symmetric geometry" begin
        # Segment along x-axis, terminal perpendicular at midpoint
        t_opt, cost = optimize_bifurcation_point(
            0.0, 0.0, 0.0, 10.0, 0.0, 0.0,
            1.0, 5.0, 5.0, 0.0, gamma,
        )

        # t_opt should be near 0.5 for symmetric perpendicular terminal
        @test 0.3 < t_opt < 0.7
        @test cost > 0.0
    end

    @testset "optimize_bifurcation_point — terminal near proximal" begin
        # Terminal is near the proximal end
        t_opt, cost = optimize_bifurcation_point(
            0.0, 0.0, 0.0, 10.0, 0.0, 0.0,
            1.0, 1.0, 2.0, 0.0, gamma,
        )

        # t_opt should be near the proximal end
        @test t_opt < 0.3
    end

    @testset "optimize_bifurcation_point — terminal near distal" begin
        # Terminal is near the distal end
        t_opt, cost = optimize_bifurcation_point(
            0.0, 0.0, 0.0, 10.0, 0.0, 0.0,
            1.0, 9.0, 2.0, 0.0, gamma,
        )

        # t_opt should be near the distal end
        @test t_opt > 0.7
    end

    @testset "evaluate_all_connections! — basic" begin
        tree = VascularTree("test", 100)
        # Segments at different distances from terminal at (2, 1, 0)
        add_segment!(tree, (0.0, 0.0, 0.0), (4.0, 0.0, 0.0), 1.0, Int32(-1))   # y=0, close
        add_segment!(tree, (0.0, 10.0, 0.0), (4.0, 10.0, 0.0), 1.0, Int32(-1))  # y=10, far
        add_segment!(tree, (0.0, 20.0, 0.0), (4.0, 20.0, 0.0), 1.0, Int32(-1))  # y=20, very far

        costs = zeros(100)
        bifurc_t = zeros(100)

        evaluate_all_connections!(costs, bifurc_t, tree.segments, 3,
            2.0, 1.0, 0.0, gamma)

        # All costs should be positive
        @test all(costs[1:3] .> 0.0)
        # All bifurc_t should be in (0, 1)
        @test all(0.0 .< bifurc_t[1:3] .< 1.0)
        # Nearer segment should have lower cost (closer terminal means shorter new branch)
        @test costs[1] < costs[2]
        @test costs[1] < costs[3]
    end

    @testset "select_best_connection" begin
        costs = [5.0, 2.0, 8.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        bifurc_t = [0.3, 0.5, 0.7, 0.4, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0]

        idx, t, cost = select_best_connection(costs, bifurc_t, 5)
        @test idx == 4
        @test t ≈ 0.4
        @test cost ≈ 1.0
    end

    @testset "evaluate_all_connections! — AK vs scalar on 100 segments" begin
        rng = Random.MersenneTwister(42)
        tree = VascularTree("test", 200)
        for _ in 1:100
            px, py, pz = rand(rng, 3) .* 10.0
            dx, dy, dz = px + rand(rng) * 2.0, py + rand(rng) * 2.0, pz + rand(rng) * 2.0
            add_segment!(tree, (px, py, pz), (dx, dy, dz), 0.5, Int32(-1))
        end

        tx, ty, tz = 5.0, 5.0, 5.0
        n = tree.segments.n

        # AK version
        costs_ak = zeros(200)
        t_ak = zeros(200)
        evaluate_all_connections!(costs_ak, t_ak, tree.segments, n, tx, ty, tz, gamma)

        # Scalar version
        costs_scalar = zeros(n)
        t_scalar = zeros(n)
        for i in 1:n
            t_opt, cost = optimize_bifurcation_point(
                tree.segments.proximal_x[i], tree.segments.proximal_y[i], tree.segments.proximal_z[i],
                tree.segments.distal_x[i], tree.segments.distal_y[i], tree.segments.distal_z[i],
                tree.segments.radius[i], tx, ty, tz, gamma,
            )
            costs_scalar[i] = cost
            t_scalar[i] = t_opt
        end

        for i in 1:n
            @test costs_ak[i] ≈ costs_scalar[i] atol = 1e-10
            @test t_ak[i] ≈ t_scalar[i] atol = 1e-10
        end
    end

    @testset "select_best_connection — matches brute force" begin
        rng = Random.MersenneTwister(99)
        tree = VascularTree("test", 200)
        for _ in 1:50
            px, py, pz = rand(rng, 3) .* 10.0
            dx, dy, dz = px + rand(rng) * 2.0, py + rand(rng) * 2.0, pz + rand(rng) * 2.0
            add_segment!(tree, (px, py, pz), (dx, dy, dz), 0.5, Int32(-1))
        end

        tx, ty, tz = 5.0, 5.0, 5.0
        n = tree.segments.n

        costs = zeros(200)
        bifurc_t = zeros(200)
        evaluate_all_connections!(costs, bifurc_t, tree.segments, n, tx, ty, tz, gamma)

        idx, t, cost = select_best_connection(costs, bifurc_t, n)

        # Brute force: find minimum cost
        min_cost = Inf
        min_idx = 0
        for i in 1:n
            if costs[i] < min_cost
                min_cost = costs[i]
                min_idx = i
            end
        end

        @test idx == min_idx
        @test cost ≈ min_cost atol = 1e-10
    end

end
