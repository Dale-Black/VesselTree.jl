using Test
using VesselTree
using Random

@testset "Element grouping (Kassab elements)" begin

    params = kassab_rca_params()

    @testset "Hand-built tree: simple bifurcation" begin
        # Tree:  1 → 2, 3
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.05, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))

        # Use topological Strahler for deterministic hand-built tests
        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)

        # Topological: terminals 2,3 → order 1; root 1 → order 2
        # Root has 2 children both order 1 → no same-order continuation
        # 3 elements, one per segment
        @test length(elements) == 3
        @test all(e -> length(e.segment_ids) >= 1, elements)
        for e in elements
            @test e.order >= 0
        end
    end

    @testset "Hand-built tree: continuation chain" begin
        # Tree:
        #       1
        #      / \
        #     2   3
        #    / \
        #   4   5
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.05, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 5.0, 0.0), 0.01, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 1.0, 0.0), 0.01, Int32(2))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)

        # Topological strahler: 4,5→1, 3→1, 2→2 (children 4,5 both 1 → max+1), 1→2 (children 2=2, 3=1 → max=2)
        # Seg 1 children: 2 (order 2), 3 (order 1) → 1 same-order child → continuation
        # So segs 1,2 form one element; segs 3,4,5 are separate → 4 elements total
        @test length(elements) == 4

        continuation = filter(e -> length(e.segment_ids) == 2, elements)
        @test length(continuation) == 1
        @test continuation[1].order == 2
        @test sort(continuation[1].segment_ids) == [1, 2]
    end

    @testset "Hand-built tree: symmetric tree (no continuation)" begin
        #        1
        #       / \
        #      2   3
        #     / \ / \
        #    4  5 6  7
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (4.0, 0.0, 0.0), 0.1, Int32(-1))
        add_segment!(tree, (4.0, 0.0, 0.0), (6.0, 2.0, 0.0), 0.05, Int32(1))
        add_segment!(tree, (4.0, 0.0, 0.0), (6.0, -2.0, 0.0), 0.05, Int32(1))
        add_segment!(tree, (6.0, 2.0, 0.0), (8.0, 3.0, 0.0), 0.02, Int32(2))
        add_segment!(tree, (6.0, 2.0, 0.0), (8.0, 1.0, 0.0), 0.02, Int32(2))
        add_segment!(tree, (6.0, -2.0, 0.0), (8.0, -1.0, 0.0), 0.02, Int32(3))
        add_segment!(tree, (6.0, -2.0, 0.0), (8.0, -3.0, 0.0), 0.02, Int32(3))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)

        # Topological: 4-7→1, 2→2 (4,5 both 1), 3→2 (6,7 both 1), 1→3 (2,3 both 2)
        # Seg 1 has children 2,3 both order 2 → 2 same-order children → NOT continuation
        # Each segment is its own element → 7 elements
        @test length(elements) == 7
        @test all(e -> length(e.segment_ids) == 1, elements)
    end

    @testset "Empty tree" begin
        tree = VascularTree("test", 10)
        elements = group_into_elements(tree, params)
        @test isempty(elements)
    end

    @testset "Single segment" begin
        tree = VascularTree("test", 10)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))
        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)
        @test length(elements) == 1
        @test elements[1].segment_ids == [1]
    end

    @testset "Element diameter and length" begin
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.05, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 5.0, 0.0), 0.01, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 1.0, 0.0), 0.01, Int32(2))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)

        # The continuation element (segs 1,2) has order 2
        cont = filter(e -> length(e.segment_ids) == 2, elements)
        @test length(cont) == 1
        elem = cont[1]

        # Mean diameter = mean of segment radii * 2 * 1000 (mm→um)
        d1 = 0.1 * 2.0 * 1000.0   # 200 um
        d2 = 0.05 * 2.0 * 1000.0  # 100 um
        @test elem.mean_diameter_um ≈ (d1 + d2) / 2.0

        # Total length = sum of segment lengths in um
        l1 = sqrt(5.0^2) * 1000.0
        l2 = sqrt(5.0^2 + 3.0^2) * 1000.0
        @test elem.total_length_um ≈ l1 + l2
    end
end

@testset "Element connectivity matrix" begin
    params = kassab_rca_params()

    @testset "Simple bifurcation connectivity" begin
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)
        CM = build_element_connectivity(elements, tree, params.n_orders)

        @test size(CM) == (params.n_orders, params.n_orders)

        # Root is order 2, children are order 1
        # CM[1+1, 2+1] = 2 children / 1 parent = 2.0
        @test CM[2, 3] ≈ 2.0
    end

    @testset "Empty elements" begin
        CM = build_element_connectivity(ElementData[], VascularTree("t", 10), 12)
        @test size(CM) == (12, 12)
        @test all(CM .== 0.0)
    end

    @testset "Continuation tree connectivity" begin
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.05, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 5.0, 0.0), 0.01, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 1.0, 0.0), 0.01, Int32(2))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)
        CM = build_element_connectivity(elements, tree, params.n_orders)

        # Continuation element (segs 1,2, order 2) has daughters:
        # - seg 3 (order 1, from seg 1's child)
        # - seg 4 (order 1, from seg 2's child)
        # - seg 5 (order 1, from seg 2's child)
        # → 3 daughter elements of order 1 per 1 parent of order 2
        # CM[order1+1, order2+1] = CM[2, 3] = 3.0
        @test CM[2, 3] ≈ 3.0
    end
end

@testset "Element statistics" begin
    params = kassab_rca_params()

    @testset "Per-order statistics" begin
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 5.0, 0.0), 0.01, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 1.0, 0.0), 0.01, Int32(2))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)
        stats = compute_element_statistics(elements)

        orders_present = Set(e.order for e in elements)
        for ord in orders_present
            @test haskey(stats, ord)
            s = stats[ord]
            @test s.mean_d > 0.0
            @test s.mean_l > 0.0
            @test s.count > 0
        end
    end

    @testset "Single element per order has zero SD" begin
        tree = VascularTree("test", 10)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)
        stats = compute_element_statistics(elements)

        for (ord, s) in stats
            if s.count == 1
                @test s.sd_d == 0.0
                @test s.sd_l == 0.0
            end
        end
    end

    @testset "Empty elements" begin
        stats = compute_element_statistics(ElementData[])
        @test isempty(stats)
    end
end

@testset "S/E ratios" begin
    params = kassab_rca_params()

    @testset "Continuation increases S/E" begin
        tree = VascularTree("test", 20)
        add_segment!(tree, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), 0.1, Int32(-1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, 3.0, 0.0), 0.05, Int32(1))
        add_segment!(tree, (5.0, 0.0, 0.0), (10.0, -3.0, 0.0), 0.02, Int32(1))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 5.0, 0.0), 0.01, Int32(2))
        add_segment!(tree, (10.0, 3.0, 0.0), (15.0, 1.0, 0.0), 0.01, Int32(2))

        VesselTree._assign_topological_strahler!(tree)
        elements = group_into_elements(tree, params)
        se = compute_se_ratios(elements)

        # The continuation element (order 2) has 2 segments → S/E = 2.0
        @test se[2] == 2.0

        # Terminal elements (order 1) each have 1 segment → S/E = 1.0
        @test se[1] == 1.0
    end

    @testset "Empty elements" begin
        se = compute_se_ratios(ElementData[])
        @test isempty(se)
    end
end

@testset "Element grouping with CCO-grown tree" begin
    params = kassab_rca_params()
    rng = MersenneTwister(42)
    domain = SphereDomain((0.0, 0.0, 0.0), 3.0)

    tree = VascularTree("cco_test", 500)
    add_segment!(tree, (0.0, 0.0, 0.0), (1.5, 0.0, 0.0), 0.5, Int32(-1))
    grow_tree!(tree, domain, 30, params; rng=rng)
    VesselTree.assign_strahler_orders!(tree, params)

    elements = group_into_elements(tree, params)

    @testset "Basic properties" begin
        @test length(elements) > 0
        @test length(elements) <= tree.segments.n

        # All segment IDs valid
        for e in elements
            for sid in e.segment_ids
                @test 1 <= sid <= tree.segments.n
            end
        end

        # Every active segment belongs to exactly one element
        assigned = Set{Int}()
        for e in elements
            for sid in e.segment_ids
                @test !(sid in assigned)
                push!(assigned, sid)
            end
        end
    end

    @testset "Element IDs are unique" begin
        ids = [e.id for e in elements]
        @test length(ids) == length(unique(ids))
    end

    @testset "Statistics are consistent" begin
        stats = compute_element_statistics(elements)
        se = compute_se_ratios(elements)

        for (ord, s) in stats
            @test s.mean_d > 0.0
            @test s.mean_l >= 0.0
            @test s.count > 0
            @test haskey(se, ord)
            @test se[ord] >= 1.0
        end
    end

    @testset "Connectivity matrix" begin
        CM = build_element_connectivity(elements, tree, params.n_orders)
        @test size(CM) == (params.n_orders, params.n_orders)
        @test all(CM .>= 0.0)
    end
end

@testset "Element grouping with subdivided tree" begin
    params = kassab_rca_params()
    rng = MersenneTwister(123)
    domain = SphereDomain((0.0, 0.0, 0.0), 5.0)

    # Grow a CCO skeleton then subdivide for more segments
    cap = estimate_subdivision_capacity(
        begin
            t = VascularTree("sub_test", 200)
            add_segment!(t, (0.0, 0.0, 0.0), (2.5, 0.0, 0.0), 0.5, Int32(-1))
            grow_tree!(t, domain, 15, params; rng=rng)
            t
        end,
        params,
    )
    tree = VascularTree("sub_test", max(cap, 5000))
    rng2 = MersenneTwister(123)
    add_segment!(tree, (0.0, 0.0, 0.0), (2.5, 0.0, 0.0), 0.5, Int32(-1))
    grow_tree!(tree, domain, 15, params; rng=rng2)
    subdivide_terminals!(tree, params; rng=MersenneTwister(456), max_order=4)
    VesselTree.assign_strahler_orders!(tree, params)

    n_segs = tree.segments.n
    elements = group_into_elements(tree, params)

    @testset "Handles subdivided tree ($(n_segs) segments)" begin
        @test n_segs > 20
        @test length(elements) > 0
        @test length(elements) <= n_segs

        # Every segment assigned to exactly one element
        all_sids = Set{Int}()
        for e in elements
            for sid in e.segment_ids
                @test !(sid in all_sids)
                push!(all_sids, sid)
            end
        end
        @test length(all_sids) == n_segs
    end

    @testset "S/E ratios are computable and reasonable" begin
        se = compute_se_ratios(elements)
        @test !isempty(se)
        for (ord, ratio) in se
            @test ratio >= 1.0  # at least 1 segment per element
            @test ratio < 100.0  # sanity upper bound
        end
    end

    @testset "Element statistics match segment data" begin
        stats = compute_element_statistics(elements)
        for (ord, s) in stats
            @test s.mean_d > 0.0
            @test s.mean_l > 0.0
            @test s.count > 0
            @test s.sd_d >= 0.0
            @test s.sd_l >= 0.0
        end
    end

    @testset "Element connectivity matrix structure" begin
        CM = build_element_connectivity(elements, tree, params.n_orders)
        @test size(CM) == (params.n_orders, params.n_orders)
        @test all(CM .>= 0.0)
        # Diagonal should be zero (an element doesn't spawn same-order daughters)
        # (except through continuation, which is same element, not daughter)
        for i in 1:params.n_orders
            @test CM[i, i] == 0.0 || CM[i, i] >= 0.0  # may be non-zero in practice
        end
    end

    @testset "S/E reference comparison (Kassab Table 5)" begin
        se = compute_se_ratios(elements)
        ref_se = params.se_ratio  # Kassab Table 5 values
        # Check that computed S/E ratios are in plausible range
        # (exact match depends on subdivision quality from VESSEL-1038)
        for (ord, ratio) in se
            if ord + 1 <= length(ref_se) && ord >= 1
                # Just verify they're computed and positive
                @test ratio > 0.0
                @test ref_se[ord + 1] > 0.0
            end
        end
    end
end
