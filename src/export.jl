# JLD2 persistence and VTP export

"""
    save_tree(filename, tree::VascularTree)

Save a VascularTree to a JLD2 file. Only active segments (1:n) are stored.
"""
function save_tree(filename::AbstractString, tree::VascularTree)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    JLD2.jldsave(filename;
        name = tree.name,
        n = n,
        root_segment_id = tree.root_segment_id,
        n_terminals = tree.n_terminals,
        n_bifurcations = tree.n_bifurcations,
        n_trifurcations = tree.n_trifurcations,
        proximal_x = seg.proximal_x[1:n],
        proximal_y = seg.proximal_y[1:n],
        proximal_z = seg.proximal_z[1:n],
        distal_x = seg.distal_x[1:n],
        distal_y = seg.distal_y[1:n],
        distal_z = seg.distal_z[1:n],
        radius = seg.radius[1:n],
        seg_length = seg.seg_length[1:n],
        flow = seg.flow[1:n],
        pressure_proximal = seg.pressure_proximal[1:n],
        pressure_distal = seg.pressure_distal[1:n],
        resistance = seg.resistance[1:n],
        parent_id = topo.parent_id[1:n],
        child1_id = topo.child1_id[1:n],
        child2_id = topo.child2_id[1:n],
        child3_id = topo.child3_id[1:n],
        strahler_order = topo.strahler_order[1:n],
        generation = topo.generation[1:n],
        is_terminal = topo.is_terminal[1:n],
        junction_type = String.(topo.junction_type[1:n]),
    )
    return filename
end

"""
    load_tree(filename) -> VascularTree

Load a VascularTree from a JLD2 file.
"""
function load_tree(filename::AbstractString)
    d = JLD2.load(filename)
    n = d["n"]
    capacity = n  # exact fit

    tree = VascularTree(d["name"], capacity)
    seg = tree.segments
    topo = tree.topology

    seg.n = n
    topo.n = n
    tree.root_segment_id = d["root_segment_id"]
    tree.n_terminals = d["n_terminals"]
    tree.n_bifurcations = d["n_bifurcations"]
    tree.n_trifurcations = d["n_trifurcations"]

    seg.proximal_x[1:n] .= d["proximal_x"]
    seg.proximal_y[1:n] .= d["proximal_y"]
    seg.proximal_z[1:n] .= d["proximal_z"]
    seg.distal_x[1:n] .= d["distal_x"]
    seg.distal_y[1:n] .= d["distal_y"]
    seg.distal_z[1:n] .= d["distal_z"]
    seg.radius[1:n] .= d["radius"]
    seg.seg_length[1:n] .= d["seg_length"]
    seg.flow[1:n] .= d["flow"]
    seg.pressure_proximal[1:n] .= d["pressure_proximal"]
    seg.pressure_distal[1:n] .= d["pressure_distal"]
    seg.resistance[1:n] .= d["resistance"]

    topo.parent_id[1:n] .= d["parent_id"]
    topo.child1_id[1:n] .= d["child1_id"]
    topo.child2_id[1:n] .= d["child2_id"]
    topo.child3_id[1:n] .= d["child3_id"]
    topo.strahler_order[1:n] .= d["strahler_order"]
    topo.generation[1:n] .= d["generation"]
    topo.is_terminal[1:n] .= d["is_terminal"]
    topo.junction_type[1:n] .= Symbol.(d["junction_type"])

    return tree
end

"""
    export_centerlines_vtp(tree::VascularTree, filename) -> String

Export tree centerlines as VTK PolyData (.vtp). Each segment is a line
with radius, Strahler order, and flow as cell data. Returns the output path.
"""
function export_centerlines_vtp(tree::VascularTree, filename::AbstractString)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    # 2 points per segment (proximal + distal)
    points = zeros(Float64, 3, 2 * n)
    for i in 1:n
        points[1, 2i - 1] = seg.proximal_x[i]
        points[2, 2i - 1] = seg.proximal_y[i]
        points[3, 2i - 1] = seg.proximal_z[i]
        points[1, 2i] = seg.distal_x[i]
        points[2, 2i] = seg.distal_y[i]
        points[3, 2i] = seg.distal_z[i]
    end

    # Each segment is a line cell connecting its two points
    cells = [MeshCell(PolyData.Lines(), [2i - 1, 2i]) for i in 1:n]

    outpath = vtk_grid(filename, points, cells) do vtk
        vtk["radius"] = seg.radius[1:n]
        vtk["strahler_order"] = Int.(topo.strahler_order[1:n])
        vtk["flow"] = seg.flow[1:n]
        vtk["is_terminal"] = Int.(topo.is_terminal[1:n])
    end

    return first(outpath)
end

"""
    export_forest_vtp(forest::CoronaryForest, directory)

Export each tree as a separate VTP file in the given directory.
"""
function export_forest_vtp(forest::CoronaryForest, directory::AbstractString)
    isdir(directory) || mkpath(directory)
    paths = String[]
    for (name, tree) in sort(collect(forest.trees), by=x -> x[1])
        filepath = joinpath(directory, name)
        push!(paths, export_centerlines_vtp(tree, filepath))
    end
    return paths
end
