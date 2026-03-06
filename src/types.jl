# Core data structures for VesselTree.jl
# SoA layout for AK compatibility, topology on CPU

"""
    SegmentData{V<:AbstractVector}

Structure-of-Arrays storage for vessel segment geometry and hemodynamics.
`V` is the vector type (e.g., `Vector{Float64}` for CPU, `CuVector{Float64}` for GPU).
Pre-allocated with `capacity`; only the first `n` entries are active.
"""
mutable struct SegmentData{V<:AbstractVector{Float64}}
    # Geometry — proximal point
    proximal_x::V
    proximal_y::V
    proximal_z::V
    # Geometry — distal point
    distal_x::V
    distal_y::V
    distal_z::V
    # Physical properties
    radius::V
    seg_length::V
    # Hemodynamics
    flow::V
    pressure_proximal::V
    pressure_distal::V
    resistance::V
    # Bookkeeping
    n::Int
    capacity::Int
end

function SegmentData{V}(capacity::Int) where {V<:AbstractVector{Float64}}
    z() = V(undef, capacity)
    return SegmentData{V}(
        z(), z(), z(),  # proximal
        z(), z(), z(),  # distal
        z(), z(),       # radius, length
        z(), z(), z(), z(),  # hemodynamics
        0, capacity,
    )
end

SegmentData(capacity::Int) = SegmentData{Vector{Float64}}(capacity)

"""
    TreeTopology

CPU-only topology data for the vascular tree.
Stores parent/child relationships, Strahler order, generation depth,
terminal status, and junction type.
"""
mutable struct TreeTopology
    parent_id::Vector{Int32}
    child1_id::Vector{Int32}
    child2_id::Vector{Int32}
    child3_id::Vector{Int32}      # -1 for non-trifurcation
    strahler_order::Vector{Int32}
    generation::Vector{Int32}
    is_terminal::Vector{Bool}
    junction_type::Vector{Symbol}  # :none, :bifurcation, :trifurcation
    n::Int
    capacity::Int
end

function TreeTopology(capacity::Int)
    return TreeTopology(
        fill(Int32(-1), capacity),   # parent_id
        fill(Int32(-1), capacity),   # child1_id
        fill(Int32(-1), capacity),   # child2_id
        fill(Int32(-1), capacity),   # child3_id
        zeros(Int32, capacity),      # strahler_order
        zeros(Int32, capacity),      # generation
        fill(true, capacity),        # is_terminal (default: leaf)
        fill(:none, capacity),       # junction_type
        0, capacity,
    )
end

"""
    VascularTree{V<:AbstractVector{Float64}}

Composite structure holding SoA geometry + CPU topology for a vascular tree.
"""
mutable struct VascularTree{V<:AbstractVector{Float64}}
    name::String
    segments::SegmentData{V}
    topology::TreeTopology
    root_segment_id::Int32
    n_terminals::Int
    n_bifurcations::Int
    n_trifurcations::Int
end

function VascularTree(name::String, capacity::Int)
    return VascularTree{Vector{Float64}}(
        name,
        SegmentData(capacity),
        TreeTopology(capacity),
        Int32(-1),
        0, 0, 0,
    )
end

function VascularTree{V}(name::String, capacity::Int) where {V<:AbstractVector{Float64}}
    return VascularTree{V}(
        name,
        SegmentData{V}(capacity),
        TreeTopology(capacity),
        Int32(-1),
        0, 0, 0,
    )
end

"""
    n_segments(tree::VascularTree) -> Int

Number of active segments in the tree.
"""
n_segments(tree::VascularTree) = tree.segments.n

"""
    add_segment!(tree, proximal, distal, radius, parent_id) -> Int32

Add a new segment to the tree. Returns the new segment's ID (1-based).
`proximal` and `distal` are 3-tuples (x, y, z).
`parent_id` should be Int32(-1) for the root segment.

If `parent_id` is valid, updates the parent's child pointers and junction type.
"""
function add_segment!(tree::VascularTree, proximal, distal, radius::Float64, parent_id::Int32)
    seg = tree.segments
    topo = tree.topology

    seg.n == seg.capacity && error("SegmentData at capacity ($(seg.capacity)). Cannot add more segments.")

    seg.n += 1
    topo.n += 1
    id = Int32(seg.n)

    # Geometry
    seg.proximal_x[id] = proximal[1]
    seg.proximal_y[id] = proximal[2]
    seg.proximal_z[id] = proximal[3]
    seg.distal_x[id] = distal[1]
    seg.distal_y[id] = distal[2]
    seg.distal_z[id] = distal[3]

    # Physical
    seg.radius[id] = radius
    dx = distal[1] - proximal[1]
    dy = distal[2] - proximal[2]
    dz = distal[3] - proximal[3]
    seg.seg_length[id] = sqrt(dx * dx + dy * dy + dz * dz)

    # Hemodynamics (initialized to zero)
    seg.flow[id] = 0.0
    seg.pressure_proximal[id] = 0.0
    seg.pressure_distal[id] = 0.0
    seg.resistance[id] = 0.0

    # Topology
    topo.parent_id[id] = parent_id
    topo.child1_id[id] = Int32(-1)
    topo.child2_id[id] = Int32(-1)
    topo.child3_id[id] = Int32(-1)
    topo.strahler_order[id] = Int32(0)
    topo.is_terminal[id] = true
    topo.junction_type[id] = :none

    # Update parent's children
    if parent_id > 0
        parent_gen = topo.generation[parent_id]
        topo.generation[id] = parent_gen + Int32(1)

        if topo.child1_id[parent_id] == Int32(-1)
            # First child — parent was terminal, now it's not; child is terminal → net 0
            topo.child1_id[parent_id] = id
            topo.is_terminal[parent_id] = false
            topo.junction_type[parent_id] = :bifurcation
            # n_terminals unchanged: parent loses terminal, child gains it
            tree.n_bifurcations += 1
        elseif topo.child2_id[parent_id] == Int32(-1)
            # Second child — net +1 terminal
            topo.child2_id[parent_id] = id
            tree.n_terminals += 1
        elseif topo.child3_id[parent_id] == Int32(-1)
            # Third child (upgrade to trifurcation)
            topo.child3_id[parent_id] = id
            topo.junction_type[parent_id] = :trifurcation
            tree.n_bifurcations -= 1
            tree.n_trifurcations += 1
            tree.n_terminals += 1
        else
            error("Segment $(parent_id) already has 3 children. Cannot add more.")
        end
    else
        # Root segment
        topo.generation[id] = Int32(0)
        tree.root_segment_id = id
        tree.n_terminals += 1
    end

    return id
end

"""
    get_children(topology, segment_id) -> Vector{Int32}

Return the child IDs of a segment. Empty vector for terminals.
"""
function get_children(topo::TreeTopology, id::Int32)
    children = Int32[]
    topo.child1_id[id] > 0 && push!(children, topo.child1_id[id])
    topo.child2_id[id] > 0 && push!(children, topo.child2_id[id])
    topo.child3_id[id] > 0 && push!(children, topo.child3_id[id])
    return children
end
