# JLD2 persistence, VTP export, STL mesh, graph JSON, CSV

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

"""
    export_stl(tree::VascularTree, filename; circumferential_resolution=16)

Export tree as binary STL mesh. Each segment becomes a cylinder approximated
by `circumferential_resolution` facets. Vertex generation uses AK kernel.
"""
function export_stl(tree::VascularTree, filename::AbstractString; circumferential_resolution::Int=16)
    n = tree.segments.n
    seg = tree.segments
    nf = circumferential_resolution  # facets around circumference
    # Each segment produces 2*nf triangles (nf quads = 2*nf tris)
    total_triangles = 2 * nf * n

    # Pre-allocate vertex ring buffers: for each segment, compute nf vertices at each end
    # Using AK for bulk angle computation
    angles = zeros(nf)
    AK.foreachindex(angles) do j
        angles[j] = 2.0 * Float64(π) * (j - 1) / nf
    end
    cos_angles = cos.(angles)
    sin_angles = sin.(angles)

    # Write binary STL
    open(filename, "w") do io
        # 80-byte header
        header = zeros(UInt8, 80)
        hdr = "VesselTree STL"
        copyto!(header, 1, Vector{UInt8}(hdr), 1, length(hdr))
        write(io, header)
        # Number of triangles
        write(io, UInt32(total_triangles))

        for i in 1:n
            px, py, pz = seg.proximal_x[i], seg.proximal_y[i], seg.proximal_z[i]
            dx, dy, dz = seg.distal_x[i], seg.distal_y[i], seg.distal_z[i]
            r = seg.radius[i]

            # Segment axis
            ax = dx - px; ay = dy - py; az = dz - pz
            alen = sqrt(ax^2 + ay^2 + az^2)
            if alen < 1e-15
                ax, ay, az = 1.0, 0.0, 0.0
                alen = 1.0
            end
            ax /= alen; ay /= alen; az /= alen

            # Two orthonormal vectors perpendicular to axis
            u1, u2, u3 = _find_perpendicular(ax, ay, az)
            v1 = ay * u3 - az * u2
            v2 = az * u1 - ax * u3
            v3 = ax * u2 - ay * u1

            # Generate ring vertices at proximal and distal ends
            for j in 1:nf
                j2 = j % nf + 1
                c1, s1 = cos_angles[j], sin_angles[j]
                c2, s2 = cos_angles[j2], sin_angles[j2]

                # 4 vertices of the quad
                p1x = px + r * (c1 * u1 + s1 * v1)
                p1y = py + r * (c1 * u2 + s1 * v2)
                p1z = pz + r * (c1 * u3 + s1 * v3)

                p2x = px + r * (c2 * u1 + s2 * v1)
                p2y = py + r * (c2 * u2 + s2 * v2)
                p2z = pz + r * (c2 * u3 + s2 * v3)

                d1x = dx + r * (c1 * u1 + s1 * v1)
                d1y = dy + r * (c1 * u2 + s1 * v2)
                d1z = dz + r * (c1 * u3 + s1 * v3)

                d2x = dx + r * (c2 * u1 + s2 * v1)
                d2y = dy + r * (c2 * u2 + s2 * v2)
                d2z = dz + r * (c2 * u3 + s2 * v3)

                # Triangle 1: p1, d1, p2
                _write_stl_triangle(io, p1x, p1y, p1z, d1x, d1y, d1z, p2x, p2y, p2z)
                # Triangle 2: p2, d1, d2
                _write_stl_triangle(io, p2x, p2y, p2z, d1x, d1y, d1z, d2x, d2y, d2z)
            end
        end
    end
    return filename
end

function _write_stl_triangle(io::IO,
    ax, ay, az, bx, by, bz, cx, cy, cz)
    # Compute normal
    e1x = bx - ax; e1y = by - ay; e1z = bz - az
    e2x = cx - ax; e2y = cy - ay; e2z = cz - az
    nx = e1y * e2z - e1z * e2y
    ny = e1z * e2x - e1x * e2z
    nz = e1x * e2y - e1y * e2x
    nlen = sqrt(nx^2 + ny^2 + nz^2)
    if nlen > 0
        nx /= nlen; ny /= nlen; nz /= nlen
    end
    write(io, Float32(nx), Float32(ny), Float32(nz))
    write(io, Float32(ax), Float32(ay), Float32(az))
    write(io, Float32(bx), Float32(by), Float32(bz))
    write(io, Float32(cx), Float32(cy), Float32(cz))
    write(io, UInt16(0))  # attribute byte count
end

"""
    export_graph_json(tree::VascularTree, filename)

Export tree topology as JSON graph with nodes (bifurcations + terminals) and
edges (segments with length, radius, flow, order).
"""
function export_graph_json(tree::VascularTree, filename::AbstractString)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    # Build node and edge lists
    nodes = []
    edges = []

    for i in 1:n
        # Each segment's distal point is a node
        push!(nodes, Dict(
            "id" => i,
            "x" => seg.distal_x[i],
            "y" => seg.distal_y[i],
            "z" => seg.distal_z[i],
            "radius" => seg.radius[i],
            "is_terminal" => topo.is_terminal[i],
            "strahler_order" => Int(topo.strahler_order[i]),
            "junction_type" => String(topo.junction_type[i]),
        ))

        push!(edges, Dict(
            "id" => i,
            "source" => Int(topo.parent_id[i]),
            "target" => i,
            "length" => seg.seg_length[i],
            "radius" => seg.radius[i],
            "flow" => seg.flow[i],
            "strahler_order" => Int(topo.strahler_order[i]),
        ))
    end

    graph = Dict("nodes" => nodes, "edges" => edges, "name" => tree.name, "n_segments" => n)

    open(filename, "w") do io
        # Simple JSON serialization (no external JSON dependency)
        _write_json(io, graph)
    end
    return filename
end

function _write_json(io::IO, obj::Dict)
    write(io, "{")
    pairs = collect(obj)
    for (idx, (k, v)) in enumerate(pairs)
        write(io, "\"", string(k), "\":")
        _write_json(io, v)
        idx < length(pairs) && write(io, ",")
    end
    write(io, "}")
end

function _write_json(io::IO, arr::AbstractVector)
    write(io, "[")
    for (idx, v) in enumerate(arr)
        _write_json(io, v)
        idx < length(arr) && write(io, ",")
    end
    write(io, "]")
end

function _write_json(io::IO, s::AbstractString)
    write(io, "\"", s, "\"")
end

function _write_json(io::IO, b::Bool)
    write(io, b ? "true" : "false")
end

function _write_json(io::IO, x::Number)
    write(io, string(x))
end

"""
    export_wenbo_txt(tree::VascularTree, filename)

Export tree in Wenbo's text format for flow simulation. Each line has 11
whitespace-separated columns:

    node_id  parent_id  direction  diameter(μm)  length(cm)  0  0  0  x(cm)  y(cm)  z(cm)

- IDs are 0-indexed (node 0 = root with parent -1)
- Direction: "r" (right/first child) or "l" (left/second child); root = "r"
- Diameter in microns (radius_mm × 2000)
- Length in centimeters (mm × 0.1) — Wenbo's physics does `L_m = len * 0.01`
- Coordinates in centimeters (mm × 0.1)
"""
function export_wenbo_txt(tree::VascularTree, filename::AbstractString)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    # Build a direction map: is this segment child1 ("r") or child2 ("l") of its parent?
    direction = fill("r", n)
    for i in 1:n
        pid = Int(topo.parent_id[i])
        pid <= 0 && continue
        if Int(topo.child2_id[pid]) == i
            direction[i] = "l"
        end
        # child3 (trifurcation) — rare, mark as "l"
        if Int(topo.child3_id[pid]) == i
            direction[i] = "l"
        end
    end

    open(filename, "w") do io
        for i in 1:n
            node_id = i - 1  # 0-indexed
            parent_id = Int(topo.parent_id[i]) <= 0 ? -1 : Int(topo.parent_id[i]) - 1
            diam_um = seg.radius[i] * 2000.0      # radius (mm) → diameter (μm)
            len_cm = seg.seg_length[i] * 0.1       # mm → cm
            dir = direction[i]
            # Distal endpoint as node position (mm → cm)
            x = seg.distal_x[i] * 0.1
            y = seg.distal_y[i] * 0.1
            z = seg.distal_z[i] * 0.1

            println(io,
                lpad(node_id, 6), "  ", lpad(parent_id, 6), "  ", dir, "  ",
                lpad(string(round(diam_um, digits=6)), 12), "  ",
                lpad(string(round(len_cm, digits=9)), 12), "  ",
                "    0  ", "    0  ", "    0  ",
                lpad(string(round(x, digits=9)), 12), "  ",
                lpad(string(round(y, digits=9)), 12), "  ",
                lpad(string(round(z, digits=9)), 12))
        end
    end
    return filename
end

"""
    export_forest_wenbo_txt(forest::CoronaryForest, directory)

Export each tree as a separate Wenbo-format text file.
"""
function export_forest_wenbo_txt(forest::CoronaryForest, directory::AbstractString)
    isdir(directory) || mkpath(directory)
    paths = String[]
    for (name, tree) in sort(collect(forest.trees), by=x -> x[1])
        filepath = joinpath(directory, "$(name).txt")
        push!(paths, export_wenbo_txt(tree, filepath))
    end
    return paths
end

"""
    export_csv(tree::VascularTree, filename)

Export all segments as a flat CSV for spreadsheet analysis.
"""
function export_csv(tree::VascularTree, filename::AbstractString)
    n = tree.segments.n
    seg = tree.segments
    topo = tree.topology

    open(filename, "w") do io
        println(io, "id,proximal_x,proximal_y,proximal_z,distal_x,distal_y,distal_z,radius,length,flow,pressure_proximal,pressure_distal,resistance,parent_id,child1_id,child2_id,child3_id,strahler_order,generation,is_terminal,junction_type")
        for i in 1:n
            println(io,
                i, ",",
                seg.proximal_x[i], ",", seg.proximal_y[i], ",", seg.proximal_z[i], ",",
                seg.distal_x[i], ",", seg.distal_y[i], ",", seg.distal_z[i], ",",
                seg.radius[i], ",", seg.seg_length[i], ",",
                seg.flow[i], ",", seg.pressure_proximal[i], ",", seg.pressure_distal[i], ",",
                seg.resistance[i], ",",
                topo.parent_id[i], ",", topo.child1_id[i], ",", topo.child2_id[i], ",", topo.child3_id[i], ",",
                topo.strahler_order[i], ",", topo.generation[i], ",",
                topo.is_terminal[i], ",", topo.junction_type[i],
            )
        end
    end
    return filename
end
