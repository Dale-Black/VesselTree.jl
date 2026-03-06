using Test
using VesselTree
using Random
using Distributions

@testset "Kassab Morphometry" begin

    params = kassab_coronary_params()
    gamma = params.gamma

    # --- Strahler ordering ---

    @testset "assign_strahler_orders! — hand-built tree" begin
        tree = VascularTree("test", 50)
        # Root: large diameter → high order
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 2.0, Int32(-1))
        # Two daughters: medium diameter
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.3, Int32(1))

        # Set radii to known Kassab diameters (in mm, params are in um)
        # Order 5: 184um diameter = 92um radius = 0.092mm
        tree.segments.radius[1] = 0.092
        # Order 3: 56um diameter = 28um radius = 0.028mm
        tree.segments.radius[2] = 0.028
        # Order 2: 30um diameter = 15um radius = 0.015mm
        tree.segments.radius[3] = 0.015

        assign_strahler_orders!(tree, params)

        topo = tree.topology
        @test topo.strahler_order[1] == Int32(5)
        @test topo.strahler_order[2] == Int32(3)
        @test topo.strahler_order[3] == Int32(2)
    end

    @testset "assign_strahler_orders! — AK kernel on many segments" begin
        n = 200
        tree = VascularTree("test", n + 10)
        rng = MersenneTwister(42)

        # Build a chain of segments with random radii
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.1 + rand(rng) * 0.5, Int32(-1))
        for i in 2:n
            r = 0.004 + rand(rng) * 0.5  # 4um to 504um radius → 8um to 1008um diameter
            add_segment!(tree, (Float64(i - 1), 0.0, 0.0), (Float64(i), rand(rng), 0.0), r, Int32(i - 1))
        end

        assign_strahler_orders!(tree, params)

        # Verify against scalar classify_order
        seg = tree.segments
        topo = tree.topology
        for i in 1:seg.n
            diameter_um = seg.radius[i] * 2.0 * 1000.0  # mm → um
            expected_order = classify_order(params, diameter_um)
            @test topo.strahler_order[i] == Int32(expected_order)
        end
    end

    @testset "assign_strahler_orders! — boundary diameters" begin
        tree = VascularTree("test", 50)

        # Test at exact boundary thresholds
        # diameter_bounds = [5, 11, 21, 43, 77, 140, 250, 450, 800, 1400, 2500, 4500, Inf]
        # Order 0: 5-11um, Order 1: 11-21um, etc.
        test_cases = [
            (5.0, 0),    # lower bound of order 0
            (11.0, 1),   # lower bound of order 1
            (21.0, 2),   # lower bound of order 2
            (140.0, 5),  # lower bound of order 5
            (4500.0, 11), # lower bound of order 11
        ]

        for (j, (diam_um, expected_order)) in enumerate(test_cases)
            tree2 = VascularTree("test", 10)
            radius_mm = diam_um / 2.0 / 1000.0
            add_segment!(tree2, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), radius_mm, Int32(-1))
            assign_strahler_orders!(tree2, params)
            @test tree2.topology.strahler_order[1] == Int32(expected_order)
        end
    end

    @testset "assign_strahler_orders! — below minimum" begin
        tree = VascularTree("test", 10)
        # Diameter = 4um (below minimum 5um boundary)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.002, Int32(-1))  # 4um diameter
        assign_strahler_orders!(tree, params)
        @test tree.topology.strahler_order[1] == Int32(-1)
    end

    # --- Asymmetry sampling ---

    @testset "sample_asymmetry — valid range" begin
        rng = MersenneTwister(42)
        for _ in 1:100
            s = sample_asymmetry(params, rng)
            @test 0.0 < s <= 1.0
        end
    end

    @testset "sample_asymmetry — KS test against Beta(2.5, 0.8)" begin
        rng = MersenneTwister(123)
        n_samples = 5000
        samples = [sample_asymmetry(params, rng) for _ in 1:n_samples]

        # Compare against Beta(2.5, 0.8) using one-sample KS test
        beta_dist = Beta(params.asymmetry_alpha, params.asymmetry_beta)
        sorted_samples = sort(samples)

        # Manual KS statistic
        ks_stat = 0.0
        for i in 1:n_samples
            ecdf_val = i / n_samples
            cdf_val = cdf(beta_dist, sorted_samples[i])
            ks_stat = max(ks_stat, abs(ecdf_val - cdf_val))
        end

        # Critical value at 0.05 significance: ~1.36/sqrt(n)
        critical = 1.36 / sqrt(n_samples)
        @test ks_stat < critical
    end

    @testset "sample_asymmetry — median near Beta theoretical" begin
        rng = MersenneTwister(456)
        samples = [sample_asymmetry(params, rng) for _ in 1:10000]
        m = sort(samples)[5000]
        # Beta(2.5, 0.8) theoretical median ≈ 0.812
        theoretical_median = quantile(Beta(params.asymmetry_alpha, params.asymmetry_beta), 0.5)
        @test m ≈ theoretical_median atol = 0.05
    end

    # --- Daughter radii ---

    @testset "compute_daughter_radii — Murray's law" begin
        r_parent = 1.0
        asymmetry = 0.5
        r_large, r_small = compute_daughter_radii(r_parent, asymmetry, gamma)

        # Murray's law: r_parent^gamma = r_large^gamma + r_small^gamma
        @test r_parent^gamma ≈ r_large^gamma + r_small^gamma rtol = 1e-10

        # r_small / r_large = asymmetry
        @test r_small / r_large ≈ asymmetry rtol = 1e-10
    end

    @testset "compute_daughter_radii — symmetric" begin
        r_parent = 2.0
        r_large, r_small = compute_daughter_radii(r_parent, 1.0, gamma)

        # Symmetric: r_large == r_small
        @test r_large ≈ r_small rtol = 1e-10

        # Murray's law
        @test r_parent^gamma ≈ 2 * r_large^gamma rtol = 1e-10
    end

    @testset "compute_daughter_radii — extreme asymmetry" begin
        r_parent = 1.0

        # Very asymmetric: small daughter is tiny
        r_large, r_small = compute_daughter_radii(r_parent, 0.1, gamma)
        @test r_parent^gamma ≈ r_large^gamma + r_small^gamma rtol = 1e-10
        @test r_small / r_large ≈ 0.1 rtol = 1e-10
        @test r_large > r_small

        # Nearly symmetric
        r_large2, r_small2 = compute_daughter_radii(r_parent, 0.99, gamma)
        @test r_parent^gamma ≈ r_large2^gamma + r_small2^gamma rtol = 1e-10
    end

    @testset "compute_daughter_radii — various parent radii" begin
        for r_parent in [0.01, 0.1, 1.0, 10.0, 100.0]
            for asymmetry in [0.2, 0.5, 0.76, 0.9, 1.0]
                r_large, r_small = compute_daughter_radii(r_parent, asymmetry, gamma)
                @test r_parent^gamma ≈ r_large^gamma + r_small^gamma rtol = 1e-8
                @test r_small <= r_large
                @test r_large > 0.0
                @test r_small > 0.0
            end
        end
    end

    # --- Apply Kassab radii ---

    @testset "apply_kassab_radii! — Murray's law preserved" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(42)
        grow_tree!(tree, domain, 20, params; rng=rng)

        apply_kassab_radii!(tree, params; rng=MersenneTwister(99))

        seg = tree.segments
        topo = tree.topology

        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                if c1 > 0 && c2 > 0
                    lhs = seg.radius[i]^gamma
                    rhs = seg.radius[c1]^gamma + seg.radius[c2]^gamma
                    @test lhs ≈ rhs rtol = 1e-6
                end
            end
        end
    end

    @testset "apply_kassab_radii! — all radii positive" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(55)
        grow_tree!(tree, domain, 30, params; rng=rng)
        apply_kassab_radii!(tree, params; rng=MersenneTwister(66))

        for i in 1:tree.segments.n
            @test tree.segments.radius[i] > 0.0
        end
    end

    @testset "apply_kassab_radii! — daughter asymmetry varies" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        tree = VascularTree("test", 500)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))

        rng = MersenneTwister(77)
        grow_tree!(tree, domain, 30, params; rng=rng)
        apply_kassab_radii!(tree, params; rng=MersenneTwister(88))

        seg = tree.segments
        topo = tree.topology

        # Collect asymmetry ratios at bifurcations
        ratios = Float64[]
        for i in 1:seg.n
            if topo.junction_type[i] == :bifurcation
                c1 = topo.child1_id[i]
                c2 = topo.child2_id[i]
                if c1 > 0 && c2 > 0
                    r1 = seg.radius[c1]
                    r2 = seg.radius[c2]
                    push!(ratios, min(r1, r2) / max(r1, r2))
                end
            end
        end

        # Should have some variation (not all symmetric)
        if length(ratios) > 2
            @test !all(r ≈ ratios[1] for r in ratios)
        end
    end

end
