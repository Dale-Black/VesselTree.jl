using Random
using Statistics
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function segment_lengths_mm(tree::VascularTree)
    return tree.segments.seg_length[1:tree.segments.n]
end

function summarize_tree(name, tree)
    lengths = segment_lengths_mm(tree)
    println(name,
        ": seg=", tree.segments.n,
        ", terminals=", tree.n_terminals,
        ", min_len=", round(minimum(lengths), digits=4),
        " mm, med_len=", round(median(lengths), digits=4),
        " mm")
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    spec = xcat_heart_organ_spec()

    println("\nImported trees without resampling")
    raw = nrb_import_fixed_trees(nrb_path, spec; centerline_resample=false)
    for name in ["LAD", "LCX", "RCA"]
        summarize_tree(name, raw.trees[name])
    end

    println("\nImported trees with resampling")
    smooth = nrb_import_fixed_trees(
        nrb_path, spec;
        centerline_resample=true,
        centerline_min_spacing_mm=0.75,
        centerline_max_spacing_mm=1.5,
        centerline_radius_spacing_factor=0.75,
    )
    for name in ["LAD", "LCX", "RCA"]
        summarize_tree(name, smooth.trees[name])
    end
end

main()
