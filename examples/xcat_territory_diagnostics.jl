using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function build_domain_and_trees(nrb_path::AbstractString)
    spec = xcat_heart_organ_spec(; target_terminals=Dict("LAD"=>18, "LCX"=>14, "RCA"=>18), territory_fractions=Dict("LAD"=>0.40, "LCX"=>0.25, "RCA"=>0.35))
    domain, surfaces = nrb_shell_domain(nrb_path, spec; rng=MersenneTwister(321), target_interior=60_000, max_candidates=300_000, batch_size=8192)
    imported = nrb_import_fixed_trees(
        nrb_path,
        spec;
        surfaces=surfaces,
        capacity=10_000,
        centerline_resample=true,
        centerline_min_spacing_mm=0.75,
        centerline_max_spacing_mm=1.5,
        centerline_radius_spacing_factor=0.75,
    )
    trees = imported.trees
    configs = fixed_tree_configs(trees; target_terminals=Dict("LAD"=>18, "LCX"=>14, "RCA"=>18), territory_fractions=Dict("LAD"=>0.40, "LCX"=>0.25, "RCA"=>0.35))
    return domain, trees, configs
end

function inside_voxels(domain, tmap)
    pts = Tuple{Int,Int,Int,Int,NTuple{3,Float64}}[]
    nx, ny, nz = tmap.dims
    for iz in 1:nz, iy in 1:ny, ix in 1:nx
        x = tmap.origin[1] + (ix - 0.5) * tmap.cell_size
        y = tmap.origin[2] + (iy - 0.5) * tmap.cell_size
        z = tmap.origin[3] + (iz - 0.5) * tmap.cell_size
        p = (x, y, z)
        in_domain(domain, p) || continue
        flat = ix + (iy - 1) * nx + (iz - 1) * nx * ny
        push!(pts, (ix, iy, iz, flat, p))
    end
    return pts
end

function actual_fractions(tmap, inside)
    counts = Dict(name => 0 for name in tmap.tree_names)
    for (_, _, _, flat, _) in inside
        owner = tmap.tree_names[tmap.cells[flat]]
        counts[owner] += 1
    end
    total = max(length(inside), 1)
    return Dict(name => counts[name] / total for name in keys(counts))
end

function nearest_tree_name(trees, point)
    best_name = first(keys(trees))
    best_dist = Inf
    for (name, tree) in trees
        n = tree.segments.n
        buf = zeros(n)
        compute_all_distances!(buf, tree.segments, point[1], point[2], point[3], n)
        d = minimum(@view buf[1:n])
        if d < best_dist
            best_dist = d
            best_name = name
        end
    end
    return best_name, best_dist
end

function mismatch_stats(tmap, inside, trees)
    mismatches = Dict(name => 0 for name in tmap.tree_names)
    counts = Dict(name => 0 for name in tmap.tree_names)
    for (_, _, _, flat, p) in inside
        assigned = tmap.tree_names[tmap.cells[flat]]
        nearest, _ = nearest_tree_name(trees, p)
        counts[assigned] += 1
        if nearest != assigned
            mismatches[assigned] += 1
        end
    end
    return Dict(name => mismatches[name] / max(counts[name], 1) for name in keys(counts))
end

function component_stats(tmap, inside)
    nx, ny, nz = tmap.dims
    inside_set = Set((ix,iy,iz) for (ix,iy,iz,_,_) in inside)
    visited = Set{Tuple{Int,Int,Int}}()
    results = Dict(name => Int[] for name in tmap.tree_names)
    neigh = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
    for (ix, iy, iz, flat, _) in inside
        cell = (ix,iy,iz)
        cell in visited && continue
        assigned = tmap.tree_names[tmap.cells[flat]]
        stack = [cell]
        push!(visited, cell)
        size = 0
        while !isempty(stack)
            cx, cy, cz = pop!(stack)
            size += 1
            for (dx,dy,dz) in neigh
                nb = (cx+dx, cy+dy, cz+dz)
                nb in inside_set || continue
                nb in visited && continue
                nix, niy, niz = nb
                nflat = nix + (niy - 1) * nx + (niz - 1) * nx * ny
                tmap.tree_names[tmap.cells[nflat]] == assigned || continue
                push!(visited, nb)
                push!(stack, nb)
            end
        end
        push!(results[assigned], size)
    end
    return results
end

function summarize(label, tmap, inside, trees)
    println("
=== $(label) ===")
    fracs = actual_fractions(tmap, inside)
    mism = mismatch_stats(tmap, inside, trees)
    comps = component_stats(tmap, inside)
    for name in sort(collect(keys(fracs)))
        parts = sort(comps[name]; rev=true)
        largest = isempty(parts) ? 0 : first(parts)
        second = length(parts) >= 2 ? parts[2] : 0
        println(name,
            ": frac=", round(fracs[name], digits=4),
            ", mismatch=", round(mism[name], digits=4),
            ", components=", length(parts),
            ", largest=", largest,
            ", second=", second)
    end
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    domain, trees, configs = build_domain_and_trees(nrb_path)
    tmap_root = initialize_territories(domain, configs)
    inside = inside_voxels(domain, tmap_root)
    summarize("root-anchor territory", tmap_root, inside, trees)
    tmap_tree = initialize_territories(domain, configs, trees; tree_distance_weight=1.0)
    summarize("tree-distance territory", tmap_tree, inside, trees)
    for w in (0.05, 0.1, 0.2)
        tmap_blend = initialize_territories(domain, configs, trees; tree_distance_weight=w)
        summarize("blended territory weight=$(w)", tmap_blend, inside, trees)
    end
end

main()
