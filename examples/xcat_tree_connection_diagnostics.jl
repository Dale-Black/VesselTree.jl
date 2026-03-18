using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function parent_child_lengths(tree::XCATCenterlineTree)
    rows = []
    for conn in tree.connections
        parent = tree.segments[conn.parent_segment]
        child = tree.segments[conn.child_segment]
        push!(rows, (
            parent=conn.parent_segment,
            child=conn.child_segment,
            gap=conn.gap_mm,
            parent_index=conn.parent_index,
            child_index=conn.child_index,
            parent_n=length(parent.centers),
            child_n=length(child.centers),
            parent_tail=length(parent.centers) - conn.parent_index,
            child_length=xcat_centerline_length_mm(child),
        ))
    end
    return rows
end

function report(label, imported)
    println("\n=== ", label, " ===")
    for name in ["LAD", "LCX", "RCA"]
        xtree = imported.xcat_trees[name]
        println(name, ": total imported seg=", sum(max(length(seg.centers)-1, 0) for seg in values(xtree.segments)))
        for row in parent_child_lengths(xtree)
            println("  ", row.parent, " -> ", row.child,
                ", gap=", round(row.gap, digits=3),
                ", parent_index=", row.parent_index, "/", row.parent_n,
                ", parent_tail_pts=", row.parent_tail,
                ", child_pts=", row.child_n,
                ", child_len_mm=", round(row.child_length, digits=3))
        end
    end
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    report("no resampling", xcat_import_coronary_trees(nrb_path; centerline_resample=false))
    report("with resampling", xcat_import_coronary_trees(nrb_path; centerline_resample=true, centerline_min_spacing_mm=0.75, centerline_max_spacing_mm=1.5, centerline_radius_spacing_factor=0.75))
end

main()
