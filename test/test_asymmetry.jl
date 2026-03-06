using Test
using VesselTree
using Random
using Statistics

@testset "CM-Implied Asymmetry (VESSEL-1039)" begin

    params = kassab_coronary_params()

    @testset "_cm_implied_asymmetry — mean tracks element diameter ratio" begin
        rng = MersenneTwister(42)
        for parent_ord in 1:5
            for d_ord in 0:(parent_ord - 1)
                expected = params.diameter_mean_elem[d_ord + 1] / params.diameter_mean_elem[parent_ord + 1]
                vals = [VesselTree._cm_implied_asymmetry(d_ord, parent_ord, params, rng) for _ in 1:500]
                m = mean(vals)
                # Mean should be within 15% of the diameter ratio
                @test abs(m - expected) / expected < 0.15
            end
        end
    end

    @testset "_cm_implied_asymmetry — always in (0, 1)" begin
        rng = MersenneTwister(123)
        for _ in 1:1000
            d_ord = rand(rng, 0:4)
            p_ord = d_ord + rand(rng, 1:3)
            p_ord = min(p_ord, params.n_orders - 1)
            asym = VesselTree._cm_implied_asymmetry(d_ord, p_ord, params, rng)
            @test 0.0 < asym < 1.0
        end
    end

    @testset "_cm_implied_asymmetry — has variance (not deterministic)" begin
        rng = MersenneTwister(42)
        vals = [VesselTree._cm_implied_asymmetry(1, 3, params, rng) for _ in 1:100]
        @test std(vals) > 0.01
    end

    @testset "_cm_implied_asymmetry — out-of-range fallback" begin
        rng = MersenneTwister(42)
        # Order beyond n_orders should return 0.5
        @test VesselTree._cm_implied_asymmetry(0, 99, params, rng) == 0.5
        @test VesselTree._cm_implied_asymmetry(99, 1, params, rng) == 0.5
    end

    @testset "_cm_implied_asymmetry — order-dependent (lower d_order → lower asymmetry)" begin
        rng = MersenneTwister(42)
        p_ord = 5
        # daughter order 1 vs daughter order 4: order-1 should have lower ratio
        vals_low = [VesselTree._cm_implied_asymmetry(1, p_ord, params, rng) for _ in 1:200]
        vals_high = [VesselTree._cm_implied_asymmetry(4, p_ord, params, rng) for _ in 1:200]
        @test mean(vals_low) < mean(vals_high)
    end

    @testset "subdivide_terminals! — element-level asymmetry matches Kassab" begin
        # Build tree with order-3 terminal and subdivide
        tree = VascularTree("test", 10000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        # Measure element-level asymmetry at bifurcations
        assign_strahler_orders!(tree, params)
        elements = group_into_elements(tree, params)

        # Build element connectivity
        conn = build_element_connectivity(elements, tree, params.n_orders)

        # Collect asymmetry ratios at element level
        asym_ratios = Float64[]
        elem_by_id = Dict(e.id => e for e in elements)
        seg_to_elem = Dict{Int, Int}()
        for e in elements
            for sid in e.segment_ids
                seg_to_elem[sid] = e.id
            end
        end

        # For each element with children, compute daughter/parent diameter ratio
        topo = tree.topology
        for e in elements
            # Find child elements: children of last segment in element
            last_seg = e.segment_ids[end]
            child_elems = Set{Int}()
            for cfield in [topo.child1_id[last_seg], topo.child2_id[last_seg]]
                cid = Int(cfield)
                if cid > 0 && haskey(seg_to_elem, cid)
                    push!(child_elems, seg_to_elem[cid])
                end
            end
            if length(child_elems) >= 2
                for ceid in child_elems
                    ce = elem_by_id[ceid]
                    if ce.mean_diameter_um < e.mean_diameter_um
                        push!(asym_ratios, ce.mean_diameter_um / e.mean_diameter_um)
                    end
                end
            end
        end

        if length(asym_ratios) >= 5
            med = median(asym_ratios)
            # CM-implied asymmetry should produce median < 0.8
            # (Beta(2.5,0.8) gave ~0.875 element-level)
            @test med < 0.80
            # And should be > 0.2 (not unreasonably low)
            @test med > 0.20
        end
    end

    @testset "subdivide_terminals! — Murray's law still exact with CM asymmetry" begin
        tree = VascularTree("test", 10000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))
        n_cco = tree.segments.n

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        gamma = params.gamma
        max_dev = 0.0
        n_checked = 0
        topo = tree.topology
        for i in (n_cco + 1):tree.segments.n
            c1 = Int(topo.child1_id[i]); c2 = Int(topo.child2_id[i])
            if c1 > 0 && c2 > 0
                rp = tree.segments.radius[i]
                r1 = tree.segments.radius[c1]
                r2 = tree.segments.radius[c2]
                lhs = rp^gamma
                rhs = r1^gamma + r2^gamma
                dev = abs(lhs - rhs) / max(lhs, 1e-30)
                dev > max_dev && (max_dev = dev)
                n_checked += 1
            end
        end
        @test n_checked > 5
        @test max_dev < 0.01  # < 1% deviation
    end

    @testset "subdivide_terminals! — r_daughter < r_parent with CM asymmetry" begin
        tree = VascularTree("test", 5000)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.028, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.0075, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.015, Int32(1))
        n_cco = tree.segments.n

        rng = MersenneTwister(42)
        subdivide_terminals!(tree, params; rng=rng)

        topo = tree.topology
        violations = 0
        for i in (n_cco + 1):tree.segments.n
            c1 = Int(topo.child1_id[i]); c2 = Int(topo.child2_id[i])
            if c1 > 0 && c2 > 0
                rp = tree.segments.radius[i]
                r1 = tree.segments.radius[c1]
                r2 = tree.segments.radius[c2]
                (r1 > rp * 1.001 || r2 > rp * 1.001) && (violations += 1)
            end
        end
        @test violations == 0
    end

    @testset "_cm_implied_asymmetry — works for all three artery types" begin
        rng = MersenneTwister(42)
        for pfn in [kassab_rca_params, kassab_lad_params, kassab_lcx_params]
            p = pfn()
            for _ in 1:100
                d_ord = rand(rng, 0:2)
                p_ord = d_ord + 1
                asym = VesselTree._cm_implied_asymmetry(d_ord, p_ord, p, rng)
                @test 0.0 < asym < 1.0
            end
        end
    end

    @testset "_cm_implied_asymmetry — CV propagation produces realistic spread" begin
        rng = MersenneTwister(42)
        # Parent=3, daughter=1: expected ratio ~0.525
        vals = [VesselTree._cm_implied_asymmetry(1, 3, params, rng) for _ in 1:500]
        cv = std(vals) / mean(vals)
        # CV should be meaningful (>5%) but not excessive (<50%)
        @test cv > 0.05
        @test cv < 0.50
    end

    @testset "subdivide_terminals! — CCO+subdivision still expands tree" begin
        domain = SphereDomain((0.0, 0.0, 0.0), 10.0)
        est_tree = VascularTree("est", 500)
        add_segment!(est_tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        grow_tree!(est_tree, domain, 30, params; rng=MersenneTwister(42))
        cap = estimate_subdivision_capacity(est_tree, params)

        tree = VascularTree("test", max(cap, 1000))
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 1.0, Int32(-1))
        grow_tree!(tree, domain, 30, params; rng=MersenneTwister(42))
        n_before = tree.segments.n

        subdivide_terminals!(tree, params; rng=MersenneTwister(99))
        @test tree.segments.n > n_before * 2
    end

end
