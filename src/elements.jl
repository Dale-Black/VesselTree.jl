# Element grouping: series of same-order segments form Kassab "elements"

"""
    ElementData

A Kassab element: a series of same-order segments connected in series.
Element diameter = mean of constituent segment diameters.
Element length = sum of constituent segment lengths.
"""
struct ElementData
    id::Int
    order::Int
    segment_ids::Vector{Int}
    mean_diameter_um::Float64
    total_length_um::Float64
end

"""
    group_into_elements(tree, params) -> Vector{ElementData}

Walk the tree and group consecutive same-order segments into elements.
Strahler orders must already be assigned.

An element boundary occurs when:
- A segment has a different order than its parent
- A segment has no parent (root)
- A segment is a child at a bifurcation where the parent has a different order
"""
function group_into_elements(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    n == 0 && return ElementData[]

    # Track which element each segment belongs to
    seg_to_elem = zeros(Int, n)
    elements = ElementData[]
    next_elem_id = 1

    # Process segments in topological order (root first)
    order = _topo_order(tree)

    for k in 1:length(order)  # sequential: tree traversal
        i = order[k]
        seg_to_elem[i] > 0 && continue  # already assigned

        my_order = Int(topo.strahler_order[i])

        # Check if this segment continues its parent's element
        parent_id = Int(topo.parent_id[i])
        if parent_id > 0
            parent_ord = Int(topo.strahler_order[parent_id])
            parent_elem = seg_to_elem[parent_id]

            # Continue the parent's element only if:
            # 1. Same order as parent
            # 2. Parent has exactly one child of the same order (continuation)
            if my_order == parent_ord && parent_elem > 0
                # Check if this is the sole same-order child (continuation)
                n_same_order_siblings = 0
                c1 = Int(topo.child1_id[parent_id])
                c2 = Int(topo.child2_id[parent_id])
                c3 = Int(topo.child3_id[parent_id])
                c1 > 0 && Int(topo.strahler_order[c1]) == parent_ord && (n_same_order_siblings += 1)
                c2 > 0 && Int(topo.strahler_order[c2]) == parent_ord && (n_same_order_siblings += 1)
                c3 > 0 && Int(topo.strahler_order[c3]) == parent_ord && (n_same_order_siblings += 1)

                if n_same_order_siblings == 1
                    # This is a continuation — extend the parent's element
                    seg_to_elem[i] = parent_elem
                    continue
                end
            end
        end

        # Start a new element
        seg_to_elem[i] = next_elem_id
        next_elem_id += 1
    end

    # Build element data from assignments
    n_elements = next_elem_id - 1
    elem_segments = [Int[] for _ in 1:n_elements]
    elem_orders = zeros(Int, n_elements)

    for i in 1:n  # sequential: collect segments per element
        eid = seg_to_elem[i]
        eid <= 0 && continue
        push!(elem_segments[eid], i)
        elem_orders[eid] = Int(topo.strahler_order[i])
    end

    for eid in 1:n_elements
        sids = elem_segments[eid]
        isempty(sids) && continue

        # Element diameter = mean of segment diameters (um)
        sum_d = 0.0
        sum_l = 0.0
        for sid in sids
            sum_d += seg.radius[sid] * 2.0 * 1000.0  # mm → um
            sum_l += seg.seg_length[sid] * 1000.0     # mm → um
        end
        mean_d = sum_d / length(sids)

        push!(elements, ElementData(eid, elem_orders[eid], sids, mean_d, sum_l))
    end

    return elements
end

"""
    build_element_connectivity(elements, n_orders) -> Matrix{Float64}

Build connectivity matrix from elements. CM[daughter_order+1, parent_order+1] =
average number of daughter elements of order (daughter_order) per parent element
of order (parent_order).
"""
function build_element_connectivity(
    elements::Vector{ElementData},
    tree::VascularTree,
    n_orders::Int,
)
    topo = tree.topology
    isempty(elements) && return zeros(Float64, n_orders, n_orders)

    # Build element ID lookup: segment_id → element
    seg_to_elem = Dict{Int, Int}()
    elem_by_id = Dict{Int, ElementData}()
    for elem in elements
        elem_by_id[elem.id] = elem
        for sid in elem.segment_ids
            seg_to_elem[sid] = elem.id
        end
    end

    # Count daughter elements per parent element
    daughter_counts = zeros(Float64, n_orders, n_orders)  # [daughter+1, parent+1]
    parent_elem_counts = zeros(Int, n_orders)

    for elem in elements  # sequential: element-level traversal
        p_order = elem.order
        p_order < 0 && continue
        p_order >= n_orders && continue

        # Find daughter elements: elements whose first segment's parent is in this element
        daughter_orders = Int[]
        for sid in elem.segment_ids
            # Check children of this segment
            c1 = Int(topo.child1_id[sid])
            c2 = Int(topo.child2_id[sid])
            c3 = Int(topo.child3_id[sid])

            for cid in (c1, c2, c3)
                cid <= 0 && continue
                child_eid = get(seg_to_elem, cid, 0)
                child_eid <= 0 && continue
                child_eid == elem.id && continue  # same element (continuation)
                child_elem = elem_by_id[child_eid]
                push!(daughter_orders, child_elem.order)
            end
        end

        if !isempty(daughter_orders)
            parent_elem_counts[p_order + 1] += 1
            for d_order in daughter_orders
                d_order < 0 && continue
                d_order >= n_orders && continue
                daughter_counts[d_order + 1, p_order + 1] += 1.0
            end
        end
    end

    # Normalize: divide by number of parent elements of each order
    CM = zeros(Float64, n_orders, n_orders)
    for j in 1:n_orders
        if parent_elem_counts[j] > 0
            for i in 1:n_orders
                CM[i, j] = daughter_counts[i, j] / parent_elem_counts[j]
            end
        end
    end

    return CM
end

"""
    compute_element_statistics(elements) -> Dict{Int, NamedTuple}

Compute per-order diameter and length statistics from elements.
Returns dict mapping order → (mean_d, sd_d, mean_l, sd_l, count).
"""
function compute_element_statistics(elements::Vector{ElementData})
    order_data = Dict{Int, Vector{Tuple{Float64, Float64}}}()

    for elem in elements  # sequential: collect data
        ord = elem.order
        ord < 0 && continue
        if !haskey(order_data, ord)
            order_data[ord] = Tuple{Float64, Float64}[]
        end
        push!(order_data[ord], (elem.mean_diameter_um, elem.total_length_um))
    end

    result = Dict{Int, NamedTuple{(:mean_d, :sd_d, :mean_l, :sd_l, :count), Tuple{Float64, Float64, Float64, Float64, Int}}}()
    for (ord, data) in order_data
        n = length(data)
        n == 0 && continue

        sum_d = sum(d for (d, _) in data)
        sum_l = sum(l for (_, l) in data)
        mean_d = sum_d / n
        mean_l = sum_l / n

        sd_d = 0.0
        sd_l = 0.0
        if n > 1
            var_d = sum((d - mean_d)^2 for (d, _) in data) / (n - 1)
            var_l = sum((l - mean_l)^2 for (_, l) in data) / (n - 1)
            sd_d = sqrt(max(var_d, 0.0))
            sd_l = sqrt(max(var_l, 0.0))
        end

        result[ord] = (mean_d=mean_d, sd_d=sd_d, mean_l=mean_l, sd_l=sd_l, count=n)
    end

    return result
end

"""
    compute_se_ratios(elements) -> Dict{Int, Float64}

Compute segments-per-element (S/E) ratio for each Strahler order.
"""
function compute_se_ratios(elements::Vector{ElementData})
    order_seg_counts = Dict{Int, Vector{Int}}()
    for elem in elements  # sequential: collect
        ord = elem.order
        ord < 0 && continue
        if !haskey(order_seg_counts, ord)
            order_seg_counts[ord] = Int[]
        end
        push!(order_seg_counts[ord], length(elem.segment_ids))
    end

    result = Dict{Int, Float64}()
    for (ord, counts) in order_seg_counts
        result[ord] = sum(counts) / length(counts)
    end
    return result
end
