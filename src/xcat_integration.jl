function _xcat_recompute_tree_topology!(tree::VascularTree)
    topo = tree.topology
    tree.n_terminals = 0
    tree.n_bifurcations = 0
    tree.n_trifurcations = 0

    for id in 1:topo.n
        children = get_children(topo, Int32(id))
        n_children = length(children)
        if n_children == 0
            topo.is_terminal[id] = true
            topo.junction_type[id] = :none
            tree.n_terminals += 1
        elseif n_children == 1
            topo.is_terminal[id] = false
            topo.junction_type[id] = :none
        elseif n_children == 2
            topo.is_terminal[id] = false
            topo.junction_type[id] = :bifurcation
            tree.n_bifurcations += 1
        elseif n_children == 3
            topo.is_terminal[id] = false
            topo.junction_type[id] = :trifurcation
            tree.n_trifurcations += 1
        else
            error("Unexpected child count $n_children at segment $id")
        end
    end
    return tree
end

function _xcat_segment_radius(centerline::XCATCenterline, i::Int)
    i2 = min(i + 1, length(centerline.radii))
    return 0.5 * (centerline.radii[i] + centerline.radii[i2])
end

function _xcat_total_segments(tree::XCATCenterlineTree)
    total = 0
    for seg in values(tree.segments)
        total += max(length(seg.centers) - 1, 0)
    end
    return total
end

function xcat_tree_to_vascular_tree(
    tree::XCATCenterlineTree;
    capacity::Int=max(64, _xcat_total_segments(tree) * 2),
)
    vtree = VascularTree(tree.name, capacity)
    segment_point_to_id = Dict{Tuple{String, Int}, Int32}()

    root = tree.segments[tree.root_segment]
    prev_id = Int32(-1)
    for i in 1:(length(root.centers) - 1)
        seg_id = add_segment!(
            vtree,
            Tuple(root.centers[i]),
            Tuple(root.centers[i + 1]),
            _xcat_segment_radius(root, i),
            prev_id,
        )
        segment_point_to_id[(root.name, i + 1)] = seg_id
        prev_id = seg_id
    end

    remaining = copy(tree.connections)
    progress = true
    while !isempty(remaining) && progress
        progress = false
        next_remaining = XCATTreeConnection[]
        for conn in remaining
            parent_key = (conn.parent_segment, conn.parent_index)
            if !haskey(segment_point_to_id, parent_key)
                push!(next_remaining, conn)
                continue
            end

            parent_id = segment_point_to_id[parent_key]
            child = tree.segments[conn.child_segment]
            branch_start = child.centers[conn.child_index]
            prev_child_id = Int32(parent_id)
            for i in conn.child_index:(length(child.centers) - 1)
                prox = i == conn.child_index ? tree.segments[conn.parent_segment].centers[conn.parent_index] : child.centers[i]
                dist = child.centers[i + 1]
                seg_id = add_segment!(
                    vtree,
                    Tuple(prox),
                    Tuple(dist),
                    _xcat_segment_radius(child, i),
                    prev_child_id,
                )
                if i == conn.child_index
                    segment_point_to_id[(conn.child_segment, conn.child_index)] = parent_id
                end
                segment_point_to_id[(conn.child_segment, i + 1)] = seg_id
                prev_child_id = seg_id
            end
            progress = true
        end
        remaining = next_remaining
    end

    isempty(remaining) || error("Unresolved XCAT tree connections for $(tree.name): $(remaining)")
    _xcat_recompute_tree_topology!(vtree)
    return vtree
end

function xcat_trees_to_vascular_trees(
    trees::Dict{String, XCATCenterlineTree};
    include_aorta::Bool=false,
    capacity::Union{Nothing, Int}=nothing,
    capacity_map::Dict{String, Int}=Dict{String, Int}(),
)
    converted = Dict{String, VascularTree}()
    for (name, tree) in trees
        if !include_aorta && name == "AORTA"
            continue
        end
        tree_capacity = haskey(capacity_map, name) ? capacity_map[name] : something(capacity, max(64, _xcat_total_segments(tree) * 2))
        converted[name] = xcat_tree_to_vascular_tree(tree; capacity=tree_capacity)
    end
    return converted
end

function xcat_default_coronary_surface_names(; phase::String="dias")
    prefix = string(phase, "_")
    return [
        prefix * "aorta",
        prefix * "lad1",
        prefix * "lad2",
        prefix * "lad3",
        prefix * "lcx",
        prefix * "rca1",
        prefix * "rca2",
    ]
end

function xcat_import_coronary_trees(
    nrb_path::AbstractString;
    phase::String="dias",
    surfaces::Union{Nothing, Vector{XCATNurbsSurface}}=nothing,
    include_aorta::Bool=false,
    capacity::Union{Nothing, Int}=nothing,
    capacity_map::Dict{String, Int}=Dict{String, Int}(),
    centerline_resample::Bool=true,
    centerline_min_spacing_mm::Float64=0.75,
    centerline_max_spacing_mm::Float64=1.5,
    centerline_radius_spacing_factor::Float64=0.75,
    trim_root_cap::Bool=true,
    root_cap_radius_fraction::Float64=0.75,
    root_cap_max_trim_fraction::Float64=0.25,
)
    spec = xcat_heart_organ_spec(; phase=phase)
    imported = nrb_import_fixed_trees(
        nrb_path,
        spec;
        surfaces=surfaces,
        capacity=capacity,
        capacity_map=capacity_map,
        centerline_resample=centerline_resample,
        centerline_min_spacing_mm=centerline_min_spacing_mm,
        centerline_max_spacing_mm=centerline_max_spacing_mm,
        centerline_radius_spacing_factor=centerline_radius_spacing_factor,
        trim_root_cap=trim_root_cap,
        root_cap_radius_fraction=root_cap_radius_fraction,
        root_cap_max_trim_fraction=root_cap_max_trim_fraction,
    )
    if include_aorta && !isnothing(spec.root_reference_surface_name)
        centerline_map = Dict(line.name => line for line in imported.centerlines)
        if haskey(centerline_map, spec.root_reference_surface_name)
            imported.xcat_trees["AORTA"] = XCATCenterlineTree(
                "AORTA",
                Dict(spec.root_reference_surface_name => centerline_map[spec.root_reference_surface_name]),
                spec.root_reference_surface_name,
                XCATTreeConnection[],
            )
        end
    end
    return imported
end

function generate_xcat_coronary_forest(
    nrb_path::AbstractString,
    params::MorphometricParams;
    phase::String="dias",
    rng::AbstractRNG=Random.default_rng(),
    verbose::Bool=false,
    kassab::Bool=true,
    target_interior::Int=80_000,
    max_candidates::Int=500_000,
    batch_size::Int=8192,
    tree_capacity::Int=5000,
    enforce_phase1_morphometry::Bool=false,
    phase1_length_sigma::Float64=5.0,
    phase1_morphometry_by_tree::Dict{String, Bool}=Dict{String, Bool}(),
    phase1_length_sigma_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    phase1_use_fixed_handoff_terminal_radius::Bool=true,
    territory_tree_distance_weight::Float64=0.0,
    centerline_resample::Bool=true,
    centerline_min_spacing_mm::Float64=0.75,
    centerline_max_spacing_mm::Float64=1.5,
    centerline_radius_spacing_factor::Float64=0.75,
    trim_root_cap::Bool=true,
    root_cap_radius_fraction::Float64=0.75,
    root_cap_max_trim_fraction::Float64=0.25,
    target_terminals::Dict{String, Int}=Dict(
        "LAD" => 10,
        "LCX" => 8,
        "RCA" => 8,
    ),
    territory_fractions::Dict{String, Float64}=Dict(
        "LAD" => 0.40,
        "LCX" => 0.25,
        "RCA" => 0.35,
    ),
)
    spec = xcat_heart_organ_spec(
        ;
        phase=phase,
        target_terminals=target_terminals,
        territory_fractions=territory_fractions,
    )
    return generate_nrb_continuation_forest(
        nrb_path,
        spec,
        params;
        rng=rng,
        verbose=verbose,
        kassab=kassab,
        target_interior=target_interior,
        max_candidates=max_candidates,
        batch_size=batch_size,
        tree_capacity=tree_capacity,
        enforce_phase1_morphometry=enforce_phase1_morphometry,
        phase1_length_sigma=phase1_length_sigma,
        phase1_morphometry_by_tree=phase1_morphometry_by_tree,
        phase1_length_sigma_by_tree=phase1_length_sigma_by_tree,
        phase1_use_fixed_handoff_terminal_radius=phase1_use_fixed_handoff_terminal_radius,
        territory_tree_distance_weight=territory_tree_distance_weight,
        centerline_resample=centerline_resample,
        centerline_min_spacing_mm=centerline_min_spacing_mm,
        centerline_max_spacing_mm=centerline_max_spacing_mm,
        centerline_radius_spacing_factor=centerline_radius_spacing_factor,
        trim_root_cap=trim_root_cap,
        root_cap_radius_fraction=root_cap_radius_fraction,
        root_cap_max_trim_fraction=root_cap_max_trim_fraction,
    )
end

function generate_xcat_kassab_coronary(
    nrb_path::AbstractString,
    params::MorphometricParams;
    phase::String="dias",
    rng::AbstractRNG=Random.default_rng(),
    verbose::Bool=false,
    handoff_order::Int=6,
    subdivision_max_order::Union{Nothing, Int}=nothing,
    apply_geometry::Bool=true,
    apply_kassab_radius_refinement::Bool=false,
    enforce_full_cutoff::Bool=false,
    max_tree_capacity::Int=2_000_000,
    target_interior::Int=80_000,
    max_candidates::Int=500_000,
    batch_size::Int=8192,
    tree_capacity::Int=10_000,
    enforce_phase1_morphometry::Bool=false,
    phase1_length_sigma::Float64=5.0,
    phase1_morphometry_by_tree::Dict{String, Bool}=Dict{String, Bool}(),
    phase1_length_sigma_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    phase1_use_fixed_handoff_terminal_radius::Bool=true,
    territory_tree_distance_weight::Float64=0.0,
    centerline_resample::Bool=true,
    centerline_min_spacing_mm::Float64=0.75,
    centerline_max_spacing_mm::Float64=1.5,
    centerline_radius_spacing_factor::Float64=0.75,
    trim_root_cap::Bool=true,
    root_cap_radius_fraction::Float64=0.75,
    root_cap_max_trim_fraction::Float64=0.25,
    target_terminals::Dict{String, Int}=Dict(
        "LAD" => 300,
        "LCX" => 220,
        "RCA" => 260,
    ),
    territory_fractions::Dict{String, Float64}=Dict(
        "LAD" => 0.40,
        "LCX" => 0.25,
        "RCA" => 0.35,
    ),
)
    spec = xcat_heart_organ_spec(
        ;
        phase=phase,
        target_terminals=target_terminals,
        territory_fractions=territory_fractions,
    )
    return generate_nrb_kassab_forest(
        nrb_path,
        spec,
        params;
        rng=rng,
        verbose=verbose,
        handoff_order=handoff_order,
        subdivision_max_order=subdivision_max_order,
        apply_geometry=apply_geometry,
        apply_kassab_radius_refinement=apply_kassab_radius_refinement,
        enforce_full_cutoff=enforce_full_cutoff,
        max_tree_capacity=max_tree_capacity,
        target_interior=target_interior,
        max_candidates=max_candidates,
        batch_size=batch_size,
        tree_capacity=tree_capacity,
        enforce_phase1_morphometry=enforce_phase1_morphometry,
        phase1_length_sigma=phase1_length_sigma,
        phase1_morphometry_by_tree=phase1_morphometry_by_tree,
        phase1_length_sigma_by_tree=phase1_length_sigma_by_tree,
        phase1_use_fixed_handoff_terminal_radius=phase1_use_fixed_handoff_terminal_radius,
        territory_tree_distance_weight=territory_tree_distance_weight,
        centerline_resample=centerline_resample,
        centerline_min_spacing_mm=centerline_min_spacing_mm,
        centerline_max_spacing_mm=centerline_max_spacing_mm,
        centerline_radius_spacing_factor=centerline_radius_spacing_factor,
        trim_root_cap=trim_root_cap,
        root_cap_radius_fraction=root_cap_radius_fraction,
        root_cap_max_trim_fraction=root_cap_max_trim_fraction,
    )
end
