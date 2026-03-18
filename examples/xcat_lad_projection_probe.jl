using Random
using VesselTree

const DEFAULT_XCAT_NRB = "/home/molloi-lab/smb_mount/shared_drive/XCAT Phantom/xcat_adult_nrb_files/vmale50_heart.nrb"

function build_result(nrb_path::AbstractString)
    rng = MersenneTwister(321)
    return generate_xcat_kassab_coronary(
        nrb_path,
        kassab_coronary_params();
        phase="dias",
        rng=rng,
        verbose=false,
        handoff_order=6,
        subdivision_max_order=6,
        apply_geometry=false,
        target_interior=60_000,
        max_candidates=300_000,
        batch_size=8192,
        tree_capacity=10_000,
        target_terminals=Dict("LAD" => 18, "LCX" => 14, "RCA" => 18),
        territory_fractions=Dict("LAD" => 0.40, "LCX" => 0.25, "RCA" => 0.35),
    )
end

function main()
    nrb_path = isempty(ARGS) ? DEFAULT_XCAT_NRB : ARGS[1]
    result = build_result(nrb_path)
    tree = result.forest.trees["LAD"]
    seg = tree.segments
    topo = tree.topology
    for seg_id in (169, 175, 181, 187)
        px,py,pz = seg.proximal_x[seg_id], seg.proximal_y[seg_id], seg.proximal_z[seg_id]
        dx,dy,dz = seg.distal_x[seg_id], seg.distal_y[seg_id], seg.distal_z[seg_id]
        println("seg=", seg_id,
            " in_prox=", in_domain(result.domain, (px,py,pz)),
            " in_dist=", in_domain(result.domain, (dx,dy,dz)),
            " sd_prox=", signed_distance(result.domain, (px,py,pz)),
            " sd_dist=", signed_distance(result.domain, (dx,dy,dz)),
            " len=", seg.seg_length[seg_id],
            " order=", topo.strahler_order[seg_id])
    end
end

main()
