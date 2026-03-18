struct NRBTreeSpec
    name::String
    surface_names::Vector{String}
    connections::Vector{Tuple{String, String}}
    target_terminals::Int
    territory_fraction::Float64
end

function NRBTreeSpec(
    name::AbstractString,
    surface_names::AbstractVector{<:AbstractString};
    connections::Union{Nothing, AbstractVector{<:Tuple{<:AbstractString, <:AbstractString}}}=nothing,
    target_terminals::Int=0,
    territory_fraction::Float64=1.0,
)
    seg_names = String.(surface_names)
    conns =
        if isnothing(connections)
            [(seg_names[i], seg_names[i + 1]) for i in 1:(length(seg_names) - 1)]
        else
            [(String(parent), String(child)) for (parent, child) in connections]
        end
    return NRBTreeSpec(String(name), seg_names, conns, target_terminals, territory_fraction)
end

struct NRBOrganSpec
    name::String
    outer_surface_name::String
    cavity_surface_names::Vector{String}
    tree_specs::Vector{NRBTreeSpec}
    root_reference_surface_name::Union{Nothing, String}
end

function xcat_heart_organ_spec(
    ;
    phase::String="dias",
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
    prefix = string(phase, "_")
    tree_specs = [
        NRBTreeSpec(
            "LAD",
            [prefix * "lad1", prefix * "lad2", prefix * "lad3"];
            target_terminals=get(target_terminals, "LAD", 300),
            territory_fraction=get(territory_fractions, "LAD", 0.40),
        ),
        NRBTreeSpec(
            "LCX",
            [prefix * "lcx"];
            target_terminals=get(target_terminals, "LCX", 220),
            territory_fraction=get(territory_fractions, "LCX", 0.25),
        ),
        NRBTreeSpec(
            "RCA",
            [prefix * "rca1", prefix * "rca2"];
            target_terminals=get(target_terminals, "RCA", 260),
            territory_fraction=get(territory_fractions, "RCA", 0.35),
        ),
    ]
    return NRBOrganSpec(
        "XCAT Heart",
        prefix * "pericardium",
        xcat_default_cavity_names(; phase=phase),
        tree_specs,
        prefix * "aorta",
    )
end

function nrb_vessel_surface_names(spec::NRBOrganSpec; include_reference::Bool=true)
    names = String[]
    for tree_spec in spec.tree_specs
        append!(names, tree_spec.surface_names)
    end
    if include_reference && !isnothing(spec.root_reference_surface_name)
        push!(names, spec.root_reference_surface_name)
    end
    return unique(names)
end

function nrb_tree_target_terminals(spec::NRBOrganSpec)
    return Dict(tree_spec.name => tree_spec.target_terminals for tree_spec in spec.tree_specs)
end

function nrb_tree_territory_fractions(spec::NRBOrganSpec)
    return Dict(tree_spec.name => tree_spec.territory_fraction for tree_spec in spec.tree_specs)
end

function _nrb_surface_dict(centerlines::AbstractVector{XCATCenterline})
    return Dict(line.name => line for line in centerlines)
end

function _nrb_build_tree(
    centerline_map::Dict{String, XCATCenterline},
    tree_spec::NRBTreeSpec;
    reference_line::Union{Nothing, XCATCenterline}=nothing,
    trim_root_cap::Bool=true,
    root_cap_radius_fraction::Float64=0.75,
    root_cap_max_trim_fraction::Float64=0.25,
)
    isempty(tree_spec.surface_names) && error("Tree spec `$(tree_spec.name)` has no surface names.")
    root_name = first(tree_spec.surface_names)
    haskey(centerline_map, root_name) || error("Missing root surface `$root_name` for tree `$(tree_spec.name)`.")

    root = centerline_map[root_name]
    if reference_line !== nothing
        root = _orient_root_to_aorta(root, reference_line)
    end
    if trim_root_cap
        root = _trim_root_cap_by_radius(
            root;
            min_radius_fraction=root_cap_radius_fraction,
            max_trim_fraction=root_cap_max_trim_fraction,
        )
    end

    segments = Dict(root.name => root)
    pending = copy(tree_spec.connections)
    connections = XCATTreeConnection[]

    while !isempty(pending)
        progress = false
        next_pending = Tuple{String, String}[]
        for (parent_name, child_name) in pending
            if !haskey(segments, parent_name)
                push!(next_pending, (parent_name, child_name))
                continue
            end
            haskey(centerline_map, child_name) || error("Missing child surface `$child_name` for tree `$(tree_spec.name)`.")
            child_oriented, conn = _make_tree_connection(segments[parent_name], centerline_map[child_name])
            segments[child_oriented.name] = child_oriented
            push!(connections, conn)
            progress = true
        end
        progress || error("Could not resolve NRB tree connections for `$(tree_spec.name)`: $next_pending")
        pending = next_pending
    end

    return XCATCenterlineTree(tree_spec.name, segments, root.name, connections)
end

function nrb_shell_domain(
    nrb_path::AbstractString,
    spec::NRBOrganSpec;
    surfaces::Union{Nothing, Vector{XCATNurbsSurface}}=nothing,
    outer_n_u::Int=160,
    outer_n_v::Int=96,
    cavity_n_u::Int=120,
    cavity_n_v::Int=72,
    rng::AbstractRNG=Random.default_rng(),
    target_interior::Int=120_000,
    max_candidates::Int=400_000,
    batch_size::Int=8192,
)
    surface_set = isnothing(surfaces) ? parse_xcat_nrb(nrb_path) : surfaces
    object_map = xcat_object_dict(surface_set)

    haskey(object_map, spec.outer_surface_name) || error("Outer surface `$(spec.outer_surface_name)` not found in $nrb_path")
    outer_surface = object_map[spec.outer_surface_name]
    outer_points, outer_normals = xcat_sampled_surface_rows(
        outer_surface;
        n_u=outer_n_u,
        n_v=outer_n_v,
        orient_outward=true,
    )

    cavity_points = Matrix{Float64}[]
    cavity_normals = Matrix{Float64}[]
    actual_names = String[]
    for name in spec.cavity_surface_names
        haskey(object_map, name) || error("Cavity surface `$name` not found in $nrb_path")
        pts, nrms = xcat_sampled_surface_rows(
            object_map[name];
            n_u=cavity_n_u,
            n_v=cavity_n_v,
            orient_outward=true,
        )
        push!(cavity_points, pts)
        push!(cavity_normals, nrms)
        push!(actual_names, name)
    end

    domain = _shell_domain_from_matrices(
        outer_points,
        outer_normals,
        cavity_points,
        cavity_normals;
        cavity_names=actual_names,
        rng=rng,
        target_interior=target_interior,
        max_candidates=max_candidates,
        batch_size=batch_size,
        coordinate_scale=1.0,
    )

    return domain, surface_set
end

function nrb_import_fixed_trees(
    nrb_path::AbstractString,
    spec::NRBOrganSpec;
    surfaces::Union{Nothing, Vector{XCATNurbsSurface}}=nothing,
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
    surface_set = isnothing(surfaces) ? parse_xcat_nrb(nrb_path) : surfaces
    vessel_surfaces = xcat_select_objects(
        surface_set;
        names=nrb_vessel_surface_names(spec; include_reference=true),
    )
    centerlines = XCATCenterline[]
    for surface in vessel_surfaces
        line = xcat_centerline_from_surface(surface)
        if centerline_resample
            line = xcat_resample_centerline(
                line;
                min_spacing_mm=centerline_min_spacing_mm,
                max_spacing_mm=centerline_max_spacing_mm,
                radius_spacing_factor=centerline_radius_spacing_factor,
            )
        end
        push!(centerlines, line)
    end
    centerline_map = _nrb_surface_dict(centerlines)
    reference_line =
        isnothing(spec.root_reference_surface_name) ? nothing :
        (haskey(centerline_map, spec.root_reference_surface_name) ? centerline_map[spec.root_reference_surface_name] : nothing)

    xcat_trees = Dict{String, XCATCenterlineTree}()
    for tree_spec in spec.tree_specs
        xcat_trees[tree_spec.name] = _nrb_build_tree(
            centerline_map,
            tree_spec;
            reference_line=reference_line,
            trim_root_cap=trim_root_cap,
            root_cap_radius_fraction=root_cap_radius_fraction,
            root_cap_max_trim_fraction=root_cap_max_trim_fraction,
        )
    end

    fixed_trees = xcat_trees_to_vascular_trees(
        xcat_trees;
        include_aorta=false,
        capacity=capacity,
        capacity_map=capacity_map,
    )

    return (
        trees=fixed_trees,
        vessel_surfaces=vessel_surfaces,
        centerlines=centerlines,
        xcat_trees=xcat_trees,
        surfaces=surface_set,
        spec=spec,
    )
end

function _resolve_tree_params(
    spec::NRBOrganSpec,
    default_params::MorphometricParams,
    tree_params::Dict{String, MorphometricParams},
)
    resolved = Dict{String, MorphometricParams}()
    for tree_spec in spec.tree_specs
        resolved[tree_spec.name] = get(tree_params, tree_spec.name, _params_for_tree(tree_spec.name, default_params))
    end
    return resolved
end

function _scaled_fixed_prefix_radii(
    fixed_trees::Dict{String, VascularTree},
    root_target_diameter_um_by_tree::Dict{String, Float64},
)
    scaled = Dict{String, Vector{Float64}}()
    for (name, tree) in fixed_trees
        radii = copy(tree.segments.radius[1:tree.segments.n])
        if haskey(root_target_diameter_um_by_tree, name)
            root_id = Int(tree.root_segment_id)
            root_diameter_um = tree.segments.radius[root_id] * 2000.0
            target_diameter_um = root_target_diameter_um_by_tree[name]
            if root_diameter_um > 0.0 && target_diameter_um > 0.0
                radii .*= target_diameter_um / root_diameter_um
            end
        end
        scaled[name] = radii
    end
    return scaled
end

function generate_nrb_continuation_forest(
    nrb_path::AbstractString,
    spec::NRBOrganSpec,
    params::MorphometricParams;
    rng::AbstractRNG=Random.default_rng(),
    verbose::Bool=false,
    kassab::Bool=true,
    handoff_order::Union{Nothing, Int}=nothing,
    target_interior::Int=80_000,
    max_candidates::Int=500_000,
    batch_size::Int=8192,
    tree_capacity::Int=5000,
    tree_params::Dict{String, MorphometricParams}=Dict{String, MorphometricParams}(),
    enforce_phase1_morphometry::Bool=false,
    phase1_length_sigma::Float64=5.0,
    phase1_morphometry_by_tree::Dict{String, Bool}=Dict{String, Bool}(),
    phase1_length_sigma_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    phase1_use_fixed_handoff_terminal_radius::Bool=true,
    fixed_prefix_radius_blend_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    fixed_prefix_root_target_diameter_um_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    territory_tree_distance_weight::Float64=0.0,
    centerline_resample::Bool=true,
    centerline_min_spacing_mm::Float64=0.75,
    centerline_max_spacing_mm::Float64=1.5,
    centerline_radius_spacing_factor::Float64=0.75,
    trim_root_cap::Bool=true,
    root_cap_radius_fraction::Float64=0.75,
    root_cap_max_trim_fraction::Float64=0.25,
)
    domain, surfaces = nrb_shell_domain(
        nrb_path,
        spec;
        rng=rng,
        target_interior=target_interior,
        max_candidates=max_candidates,
        batch_size=batch_size,
    )

    imported = nrb_import_fixed_trees(
        nrb_path,
        spec;
        surfaces=surfaces,
        capacity=tree_capacity,
        centerline_resample=centerline_resample,
        centerline_min_spacing_mm=centerline_min_spacing_mm,
        centerline_max_spacing_mm=centerline_max_spacing_mm,
        centerline_radius_spacing_factor=centerline_radius_spacing_factor,
        trim_root_cap=trim_root_cap,
        root_cap_radius_fraction=root_cap_radius_fraction,
        root_cap_max_trim_fraction=root_cap_max_trim_fraction,
    )
    fixed_trees = imported.trees
    fixed_prefix_radii = _scaled_fixed_prefix_radii(fixed_trees, fixed_prefix_root_target_diameter_um_by_tree)
    working_trees = deepcopy(fixed_trees)
    per_tree_params = _resolve_tree_params(spec, params, tree_params)
    tree_configs = fixed_tree_configs(
        fixed_trees;
        target_terminals=nrb_tree_target_terminals(spec),
        territory_fractions=nrb_tree_territory_fractions(spec),
    )
    forest = continue_coronary_forest(
        domain,
        working_trees,
        params;
        tree_configs=tree_configs,
        rng=rng,
        verbose=verbose,
        kassab=kassab,
        tree_params=per_tree_params,
        handoff_order=handoff_order,
        enforce_phase1_morphometry=enforce_phase1_morphometry,
        phase1_length_sigma=phase1_length_sigma,
        phase1_morphometry_by_tree=phase1_morphometry_by_tree,
        phase1_length_sigma_by_tree=phase1_length_sigma_by_tree,
        phase1_use_fixed_handoff_terminal_radius=phase1_use_fixed_handoff_terminal_radius,
        fixed_prefix_radii=fixed_prefix_radii,
        fixed_prefix_radius_blend_by_tree=fixed_prefix_radius_blend_by_tree,
        territory_tree_distance_weight=territory_tree_distance_weight,
    )

    return (
        forest=forest,
        domain=domain,
        fixed_trees=fixed_trees,
        tree_configs=tree_configs,
        vessel_surfaces=imported.vessel_surfaces,
        centerlines=imported.centerlines,
        xcat_trees=imported.xcat_trees,
        surfaces=surfaces,
        spec=spec,
    )
end

function generate_nrb_kassab_forest(
    nrb_path::AbstractString,
    spec::NRBOrganSpec,
    params::MorphometricParams;
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
    tree_params::Dict{String, MorphometricParams}=Dict{String, MorphometricParams}(),
    enforce_phase1_morphometry::Bool=false,
    phase1_length_sigma::Float64=5.0,
    phase1_morphometry_by_tree::Dict{String, Bool}=Dict{String, Bool}(),
    phase1_length_sigma_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    phase1_use_fixed_handoff_terminal_radius::Bool=true,
    fixed_prefix_radius_blend_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    fixed_prefix_root_target_diameter_um_by_tree::Dict{String, Float64}=Dict{String, Float64}(),
    territory_tree_distance_weight::Float64=0.0,
    centerline_resample::Bool=true,
    centerline_min_spacing_mm::Float64=0.75,
    centerline_max_spacing_mm::Float64=1.5,
    centerline_radius_spacing_factor::Float64=0.75,
    trim_root_cap::Bool=true,
    root_cap_radius_fraction::Float64=0.75,
    root_cap_max_trim_fraction::Float64=0.25,
)
    t_start = time()

    domain, surfaces = nrb_shell_domain(
        nrb_path,
        spec;
        rng=rng,
        target_interior=target_interior,
        max_candidates=max_candidates,
        batch_size=batch_size,
    )

    raw_import = nrb_import_fixed_trees(
        nrb_path,
        spec;
        surfaces=surfaces,
        capacity=tree_capacity,
        centerline_resample=centerline_resample,
        centerline_min_spacing_mm=centerline_min_spacing_mm,
        centerline_max_spacing_mm=centerline_max_spacing_mm,
        centerline_radius_spacing_factor=centerline_radius_spacing_factor,
        trim_root_cap=trim_root_cap,
        root_cap_radius_fraction=root_cap_radius_fraction,
        root_cap_max_trim_fraction=root_cap_max_trim_fraction,
    )
    per_tree_params = _resolve_tree_params(spec, params, tree_params)

    capacity_map = Dict{String, Int}()
    for tree_spec in spec.tree_specs
        name = tree_spec.name
        xtree = raw_import.xcat_trees[name]
        tp = per_tree_params[name]
        ho = isnothing(subdivision_max_order) ? min(handoff_order, tp.n_orders - 1) : min(subdivision_max_order, tp.n_orders - 1)
        imported_segments = _xcat_total_segments(xtree)
        continuation_capacity = imported_segments + tree_spec.target_terminals * 3
        subdivision_expansion = max(estimate_total_segments(ho, tp), 1)
        estimated_total = max(continuation_capacity * subdivision_expansion * 3, tree_capacity)
        capped_total = min(estimated_total, max_tree_capacity)
        capacity_map[name] = capped_total
        if verbose && capped_total < estimated_total
            println("  $(name): estimated capacity $(estimated_total) capped at $(max_tree_capacity)")
        end
    end

    imported = nrb_import_fixed_trees(
        nrb_path,
        spec;
        surfaces=surfaces,
        capacity_map=capacity_map,
        centerline_resample=centerline_resample,
        centerline_min_spacing_mm=centerline_min_spacing_mm,
        centerline_max_spacing_mm=centerline_max_spacing_mm,
        centerline_radius_spacing_factor=centerline_radius_spacing_factor,
        trim_root_cap=trim_root_cap,
        root_cap_radius_fraction=root_cap_radius_fraction,
        root_cap_max_trim_fraction=root_cap_max_trim_fraction,
    )
    fixed_trees = imported.trees
    fixed_prefix_radii = _scaled_fixed_prefix_radii(fixed_trees, fixed_prefix_root_target_diameter_um_by_tree)
    working_trees = deepcopy(fixed_trees)
    tree_configs = fixed_tree_configs(
        fixed_trees;
        target_terminals=nrb_tree_target_terminals(spec),
        territory_fractions=nrb_tree_territory_fractions(spec),
    )

    if verbose
        println("Phase 1: NRB fixed-trunk continuation...")
    end
    t_phase1 = time()
    forest = continue_coronary_forest(
        domain,
        working_trees,
        params;
        tree_configs=tree_configs,
        rng=rng,
        verbose=verbose,
        kassab=true,
        tree_params=per_tree_params,
        handoff_order=handoff_order,
        enforce_phase1_morphometry=enforce_phase1_morphometry,
        phase1_length_sigma=phase1_length_sigma,
        phase1_morphometry_by_tree=phase1_morphometry_by_tree,
        phase1_length_sigma_by_tree=phase1_length_sigma_by_tree,
        phase1_use_fixed_handoff_terminal_radius=phase1_use_fixed_handoff_terminal_radius,
        fixed_prefix_radii=fixed_prefix_radii,
        fixed_prefix_radius_blend_by_tree=fixed_prefix_radius_blend_by_tree,
        territory_tree_distance_weight=territory_tree_distance_weight,
    )
    t_phase1_end = time()
    if verbose
        println("  Phase 1 time: $(round(t_phase1_end - t_phase1, digits=1))s")
        println("Phase 2: Statistical subdivision...")
    end

    t_phase2 = time()
    for tree_spec in spec.tree_specs
        name = tree_spec.name
        tree = forest.trees[name]
        tp = per_tree_params[name]
        ho = isnothing(subdivision_max_order) ? min(handoff_order, tp.n_orders - 1) : min(subdivision_max_order, tp.n_orders - 1)
        subdivide_terminals!(tree, tp; rng=rng, max_order=ho, domain=domain, enforce_full_cutoff=enforce_full_cutoff)
        if verbose
            println("  $(name): $(tree.segments.n) segments after subdivision")
        end
    end

    for tree_spec in spec.tree_specs
        name = tree_spec.name
        tree = forest.trees[name]
        tp = per_tree_params[name]
        seg = tree.segments
        topo = tree.topology
        for i in 1:seg.n
            if topo.is_terminal[i]
                ord = Int(topo.strahler_order[i])
                if ord <= 1 && ord + 1 <= length(tp.diameter_mean_elem)
                    seg.radius[i] = tp.diameter_mean_elem[ord + 1] / 2.0 / 1000.0
                end
            end
        end
        update_radii!(tree, tp.gamma)
        if apply_kassab_radius_refinement
            apply_full_kassab_radii!(tree, tp; rng=rng)
        end
        if haskey(fixed_prefix_radii, name)
            alpha = get(fixed_prefix_radius_blend_by_tree, name, 0.0)
            alpha > 0.0 && _apply_fixed_prefix_radii_blend!(tree, fixed_prefix_radii[name], alpha)
        end
    end
    t_phase2_end = time()

    if verbose && apply_geometry
        println("  Phase 2 time: $(round(t_phase2_end - t_phase2, digits=1))s")
        println("Phase 3: Barabasi geometry...")
    elseif verbose
        println("  Phase 2 time: $(round(t_phase2_end - t_phase2, digits=1))s")
        println("Phase 3: skipped geometry")
    end

    t_phase3 = time()
    if apply_geometry
        for tree_spec in spec.tree_specs
            name = tree_spec.name
            tree = forest.trees[name]
            tp = per_tree_params[name]
            topo_ord = _topo_order(tree)
            for i in topo_ord
                if tree.topology.junction_type[i] == :bifurcation
                    apply_junction_geometry!(tree, Int32(i), tp)
                elseif tree.topology.junction_type[i] == :trifurcation
                    apply_trifurcation_geometry!(tree, Int32(i), tp)
                end
            end
            _project_tree_to_domain!(tree, domain)
            if verbose
                println("  $(name): geometry applied + projected to domain")
            end
        end
    end
    t_phase3_end = time()

    result = (
        forest=CoronaryForest(forest.trees, forest.territory_map, params, per_tree_params),
        domain=domain,
        fixed_trees=fixed_trees,
        tree_configs=tree_configs,
        vessel_surfaces=imported.vessel_surfaces,
        centerlines=imported.centerlines,
        xcat_trees=imported.xcat_trees,
        surfaces=surfaces,
        spec=spec,
    )

    if verbose
        println("  Phase 3 time: $(round(t_phase3_end - t_phase3, digits=1))s")
        println("Phase 4: Summary")
        for (name, tree) in sort(collect(result.forest.trees), by=x -> x[1])
            println("  $(name): $(tree.segments.n) segments, $(tree.n_terminals) terminals")
        end
        println("  Total time: $(round(time() - t_start, digits=1))s")
    end

    return result
end

function apply_vasodilation!(
    tree::VascularTree;
    factor::Float64=1.6,
    min_diameter_um::Float64=0.0,
    max_diameter_um::Float64=400.0,
)
    seg = tree.segments
    for i in 1:seg.n
        diameter_um = seg.radius[i] * 2000.0
        if min_diameter_um <= diameter_um <= max_diameter_um
            seg.radius[i] *= factor
        end
    end
    return tree
end

function compute_root_flow_mLmin(tree::VascularTree, params::MorphometricParams)
    compute_resistances!(tree, params.blood_viscosity)
    compute_flows!(tree, params)
    compute_pressures!(tree, params)
    return tree.segments.flow[Int(tree.root_segment_id)] * (1e6 * 60.0)
end
