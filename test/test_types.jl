using Test
using VesselTree

@testset "Core Types" begin

    @testset "SegmentData construction" begin
        # Default (CPU Vector) construction
        seg = SegmentData(100)
        @test seg.n == 0
        @test seg.capacity == 100
        @test length(seg.proximal_x) == 100
        @test length(seg.radius) == 100
        @test length(seg.flow) == 100

        # Parameterized construction
        seg2 = SegmentData{Vector{Float64}}(50)
        @test seg2.capacity == 50
        @test seg2.n == 0
    end

    @testset "TreeTopology construction" begin
        topo = TreeTopology(100)
        @test topo.n == 0
        @test topo.capacity == 100
        @test all(topo.parent_id .== Int32(-1))
        @test all(topo.child1_id .== Int32(-1))
        @test all(topo.is_terminal .== true)
        @test all(topo.junction_type .== :none)
        @test all(topo.strahler_order .== Int32(0))
    end

    @testset "VascularTree construction" begin
        tree = VascularTree("test_artery", 200)
        @test tree.name == "test_artery"
        @test n_segments(tree) == 0
        @test tree.root_segment_id == Int32(-1)
        @test tree.n_terminals == 0
        @test tree.n_bifurcations == 0
        @test tree.n_trifurcations == 0
        @test tree.segments.capacity == 200
        @test tree.topology.capacity == 200
    end

    @testset "add_segment! — root segment" begin
        tree = VascularTree("test", 100)
        proximal = (0.0, 0.0, 0.0)
        distal = (1.0, 0.0, 0.0)

        id = add_segment!(tree, proximal, distal, 0.5, Int32(-1))

        @test id == Int32(1)
        @test n_segments(tree) == 1
        @test tree.root_segment_id == Int32(1)
        @test tree.n_terminals == 1
        @test tree.n_bifurcations == 0

        # Check geometry
        @test tree.segments.proximal_x[1] ≈ 0.0
        @test tree.segments.distal_x[1] ≈ 1.0
        @test tree.segments.radius[1] ≈ 0.5
        @test tree.segments.seg_length[1] ≈ 1.0

        # Check topology
        @test tree.topology.parent_id[1] == Int32(-1)
        @test tree.topology.is_terminal[1] == true
        @test tree.topology.generation[1] == Int32(0)
    end

    @testset "add_segment! — bifurcation" begin
        tree = VascularTree("test", 100)

        # Root
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))

        # First child — creates bifurcation on parent
        id2 = add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.7, Int32(1))
        @test id2 == Int32(2)
        @test tree.topology.child1_id[1] == Int32(2)
        @test tree.topology.is_terminal[1] == false
        @test tree.topology.junction_type[1] == :bifurcation
        @test tree.n_bifurcations == 1
        @test tree.n_terminals == 1  # parent lost terminal, child gained it → net 0

        # Second child
        id3 = add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.7, Int32(1))
        @test id3 == Int32(3)
        @test tree.topology.child2_id[1] == Int32(3)
        @test tree.topology.junction_type[1] == :bifurcation
        @test tree.n_terminals == 2  # net +1 from second child
        @test tree.n_bifurcations == 1

        # Verify children
        children = get_children(tree.topology, Int32(1))
        @test length(children) == 2
        @test Int32(2) in children
        @test Int32(3) in children
    end

    @testset "add_segment! — trifurcation" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.5, Int32(1))

        # Third child — upgrade to trifurcation
        id4 = add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 0.0, 1.0), 0.3, Int32(1))
        @test id4 == Int32(4)
        @test tree.topology.child3_id[1] == Int32(4)
        @test tree.topology.junction_type[1] == :trifurcation
        @test tree.n_trifurcations == 1
        @test tree.n_bifurcations == 0  # upgraded from bifurcation

        children = get_children(tree.topology, Int32(1))
        @test length(children) == 3
    end

    @testset "add_segment! — generation tracking" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), 0.7, Int32(1))
        add_segment!(tree, (2.0, 0.0, 0.0), (3.0, 0.0, 0.0), 0.5, Int32(2))

        @test tree.topology.generation[1] == Int32(0)
        @test tree.topology.generation[2] == Int32(1)
        @test tree.topology.generation[3] == Int32(2)
    end

    @testset "add_segment! — segment length 3D" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (1.0, 2.0, 3.0), (4.0, 6.0, 3.0), 0.5, Int32(-1))

        # sqrt((4-1)^2 + (6-2)^2 + (3-3)^2) = sqrt(9+16) = 5.0
        @test tree.segments.seg_length[1] ≈ 5.0
    end

    @testset "add_segment! — capacity error" begin
        tree = VascularTree("test", 2)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), 0.5, Int32(1))

        @test_throws ErrorException add_segment!(tree, (2.0, 0.0, 0.0), (3.0, 0.0, 0.0), 0.3, Int32(2))
    end

    @testset "add_segment! — max children error" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.5, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 0.0, 1.0), 0.3, Int32(1))

        # Fourth child should error
        @test_throws ErrorException add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 0.0, -1.0), 0.2, Int32(1))
    end

    @testset "get_children — terminal node" begin
        tree = VascularTree("test", 100)
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))

        children = get_children(tree.topology, Int32(1))
        @test isempty(children)
    end

    @testset "SoA views for AK" begin
        tree = VascularTree("test", 1000)
        for i in 1:10
            add_segment!(tree, (0.0, 0.0, 0.0), (Float64(i), 0.0, 0.0), 0.1 * i, Int32(-1))
        end

        # Active view should have exactly n elements
        n = tree.segments.n
        @test n == 10
        radii_view = @view tree.segments.radius[1:n]
        @test length(radii_view) == 10
        @test radii_view[5] ≈ 0.5

        # Can use AK operations on views
        import AcceleratedKernels as AK
        max_r = AK.maximum(radii_view)
        @test max_r ≈ 1.0
    end

    @testset "Multi-level tree" begin
        tree = VascularTree("coronary", 100)

        # Build a small tree:
        #       1
        #      / \
        #     2   3
        #    / \
        #   4   5
        add_segment!(tree, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 1.0, Int32(-1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), 0.7, Int32(1))
        add_segment!(tree, (1.0, 0.0, 0.0), (2.0, -1.0, 0.0), 0.7, Int32(1))
        add_segment!(tree, (2.0, 1.0, 0.0), (3.0, 2.0, 0.0), 0.5, Int32(2))
        add_segment!(tree, (2.0, 1.0, 0.0), (3.0, 0.0, 0.0), 0.5, Int32(2))

        @test n_segments(tree) == 5
        @test tree.n_terminals == 3   # segments 3, 4, 5
        @test tree.n_bifurcations == 2  # segments 1 and 2
        @test tree.n_trifurcations == 0

        # Check terminal flags
        @test tree.topology.is_terminal[1] == false
        @test tree.topology.is_terminal[2] == false
        @test tree.topology.is_terminal[3] == true
        @test tree.topology.is_terminal[4] == true
        @test tree.topology.is_terminal[5] == true

        # Check parent chain
        @test tree.topology.parent_id[4] == Int32(2)
        @test tree.topology.parent_id[2] == Int32(1)
        @test tree.topology.parent_id[1] == Int32(-1)
    end

end
