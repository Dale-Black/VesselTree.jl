# Kassab morphometry: Strahler ordering, asymmetry sampling, daughter radii

"""
    assign_strahler_orders!(tree::VascularTree, params::MorphometricParams)

Assign Strahler order to each segment based on its diameter using AK kernel.
Diameter is converted from mm (internal) to um (Kassab scale) and classified
using params.diameter_bounds.
"""
function assign_strahler_orders!(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    n == 0 && return nothing

    rad_view = @view seg.radius[1:n]
    order_view = @view topo.strahler_order[1:n]

    # Diameter bounds for classification (in um)
    bounds = params.diameter_bounds
    n_orders = params.n_orders

    AK.foreachindex(order_view) do i
        diameter_um = rad_view[i] * 2.0 * 1000.0  # mm → um
        ord = Int32(-1)
        # Scan from highest order downward
        for k in n_orders:-1:1
            if diameter_um >= bounds[k]
                ord = Int32(k - 1)  # 0-indexed
                break
            end
        end
        order_view[i] = ord
    end
    return nothing
end

"""
    sample_asymmetry(params::MorphometricParams, rng::AbstractRNG) -> Float64

Sample asymmetry ratio (r_small/r_large) from Beta distribution.
Clamped to (0, 1].
"""
function sample_asymmetry(params::MorphometricParams, rng::AbstractRNG)
    d = Beta(params.asymmetry_alpha, params.asymmetry_beta)
    s = rand(rng, d)
    return clamp(s, 1e-6, 1.0)
end

"""
    compute_daughter_radii(r_parent, asymmetry, gamma) -> (r_large, r_small)

Compute daughter radii from parent radius and asymmetry ratio.
Murray's law: r_parent^gamma = r_large^gamma + r_small^gamma
Asymmetry: r_small = asymmetry * r_large

Returns (r_large, r_small) where r_large >= r_small.
"""
function compute_daughter_radii(r_parent::Float64, asymmetry::Float64, gamma::Float64)
    # r_parent^gamma = r_large^gamma + (asymmetry * r_large)^gamma
    # r_parent^gamma = r_large^gamma * (1 + asymmetry^gamma)
    # r_large = r_parent / (1 + asymmetry^gamma)^(1/gamma)
    r_large = r_parent / (1.0 + asymmetry^gamma)^(1.0 / gamma)
    r_small = asymmetry * r_large
    return (r_large, r_small)
end

"""
    apply_kassab_radii!(tree::VascularTree, params::MorphometricParams; rng=Random.default_rng())

For each bifurcation, sample asymmetry from Kassab's Beta distribution and
assign daughter radii. Then propagate via Murray's law bottom-up.
"""
function apply_kassab_radii!(
    tree::VascularTree,
    params::MorphometricParams;
    rng::AbstractRNG=Random.default_rng(),
)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    gamma = params.gamma

    # sequential: tree traversal — process top-down to set asymmetry at each bifurcation
    order = VesselTree._topo_order(tree)

    # First pass: assign asymmetry-based daughter radii top-down
    for k in 1:length(order)
        i = order[k]
        if topo.junction_type[i] == :bifurcation
            c1 = Int(topo.child1_id[i])
            c2 = Int(topo.child2_id[i])
            if c1 > 0 && c2 > 0
                asymmetry = sample_asymmetry(params, rng)
                r_large, r_small = compute_daughter_radii(seg.radius[i], asymmetry, gamma)
                # Assign larger radius to larger existing daughter
                if seg.radius[c1] >= seg.radius[c2]
                    seg.radius[c1] = r_large
                    seg.radius[c2] = r_small
                else
                    seg.radius[c1] = r_small
                    seg.radius[c2] = r_large
                end
            end
        end
    end

    # Second pass: propagate bottom-up to ensure Murray's law consistency
    update_radii!(tree, gamma)

    return nothing
end
