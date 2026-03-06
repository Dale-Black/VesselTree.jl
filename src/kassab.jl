# Kassab morphometry: Strahler ordering, asymmetry sampling, daughter radii

"""
    assign_strahler_orders!(tree::VascularTree, params::MorphometricParams)

Assign Strahler order to each segment based on its diameter using AK kernel.
Diameter is converted from mm (internal) to um (Kassab scale) and classified
using params.diameter_bounds. This is the simple fixed-bound method.
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

"""Alias for assign_strahler_orders! — the simple fixed-bound method."""
const assign_strahler_orders_simple! = assign_strahler_orders!

"""
    assign_diameter_defined_strahler!(tree, params; max_iter=20) -> Int

Iterative diameter-defined Strahler ordering (Jiang et al. 1994, Eq 3A/3B).

Algorithm:
1. Initialize with topological Strahler ordering (bottom-up from terminals)
2. Compute per-order mean diameter D_n and SD_n
3. Compute non-overlapping diameter bounds from Eq 3A/3B
4. Reassign each vessel to the order whose range contains its diameter
5. Repeat until < 1% of vessels change order

Returns number of iterations taken to converge.
"""
function assign_diameter_defined_strahler!(
    tree::VascularTree,
    params::MorphometricParams;
    max_iter::Int=20,
)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    n == 0 && return 0

    # Step 1: Topological Strahler ordering (bottom-up)
    _assign_topological_strahler!(tree)

    # Get diameters in um (used throughout)
    diameters = Vector{Float64}(undef, n)
    for i in 1:n  # sequential: simple array fill
        diameters[i] = seg.radius[i] * 2.0 * 1000.0
    end

    # Iterative refinement
    for iter in 1:max_iter
        # Step 2: Compute per-order statistics
        max_order = 0
        for i in 1:n
            ord = Int(topo.strahler_order[i])
            ord > max_order && (max_order = ord)
        end
        max_order < 1 && break

        n_per_order = zeros(Int, max_order + 1)     # 0-indexed orders
        sum_d = zeros(Float64, max_order + 1)
        sum_d2 = zeros(Float64, max_order + 1)

        for i in 1:n  # sequential: accumulate statistics
            ord = Int(topo.strahler_order[i])
            ord < 0 && continue
            idx = ord + 1
            idx > length(n_per_order) && continue
            n_per_order[idx] += 1
            sum_d[idx] += diameters[i]
            sum_d2[idx] += diameters[i]^2
        end

        mean_d = zeros(Float64, max_order + 1)
        sd_d = zeros(Float64, max_order + 1)
        for k in 1:(max_order + 1)
            if n_per_order[k] > 0
                mean_d[k] = sum_d[k] / n_per_order[k]
                if n_per_order[k] > 1
                    var = (sum_d2[k] / n_per_order[k]) - mean_d[k]^2
                    sd_d[k] = sqrt(max(var, 0.0))
                end
            end
        end

        # Step 3: Compute diameter bounds (Eq 3A/3B from Jiang 1994)
        n_bounds = max_order + 2
        bounds = Vector{Float64}(undef, n_bounds)
        bounds[1] = 0.0       # D'1(0) = 0 (lowest order has no lower bound)

        for k in 2:(max_order + 1)
            # Lower bound of order k-1: D'1(k-1) = [(D_{k-2}+SD_{k-2}) + (D_{k-1}-SD_{k-1})] / 2
            d_lo = mean_d[k - 1] + sd_d[k - 1]  # D_{k-2} + SD_{k-2} (prev order)
            d_hi = mean_d[k] - sd_d[k]           # D_{k-1} - SD_{k-1} (current order)
            bounds[k] = max((d_lo + d_hi) / 2.0, bounds[k - 1] + 0.001)  # enforce monotonicity
        end
        bounds[n_bounds] = Inf  # D'2(max_order) = Inf

        # Step 4: Reassign orders by diameter range
        n_changed = 0
        for i in 1:n  # sequential: reassignment
            d = diameters[i]
            old_ord = Int(topo.strahler_order[i])
            new_ord = Int32(-1)
            for k in (max_order + 1):-1:1
                if d >= bounds[k]
                    new_ord = Int32(k - 1)
                    break
                end
            end
            if new_ord != old_ord
                topo.strahler_order[i] = new_ord
                n_changed += 1
            end
        end

        # Step 5: Check convergence
        frac_changed = n_changed / n
        frac_changed < 0.01 && return iter
    end

    return max_iter
end

"""
    _assign_topological_strahler!(tree)

Assign topological Strahler orders bottom-up: terminals get order 1,
at bifurcations parent = max(children) + (1 if equal).
"""
function _assign_topological_strahler!(tree::VascularTree)
    topo = tree.topology
    n = tree.segments.n
    n == 0 && return nothing

    # BFS order (top-down), then process bottom-up
    order = _topo_order(tree)

    # Initialize all to 0 (will be set bottom-up)
    for i in 1:n  # sequential: initialization
        topo.strahler_order[i] = Int32(0)
    end

    # Bottom-up pass
    for k in length(order):-1:1
        i = order[k]
        c1 = Int(topo.child1_id[i])
        c2 = Int(topo.child2_id[i])
        c3 = Int(topo.child3_id[i])

        if c1 <= 0
            # Terminal → order 1
            topo.strahler_order[i] = Int32(1)
        elseif c3 > 0
            # Trifurcation
            o1 = Int(topo.strahler_order[c1])
            o2 = Int(topo.strahler_order[c2])
            o3 = Int(topo.strahler_order[c3])
            m = max(o1, o2, o3)
            n_max = (o1 == m ? 1 : 0) + (o2 == m ? 1 : 0) + (o3 == m ? 1 : 0)
            topo.strahler_order[i] = Int32(n_max >= 2 ? m + 1 : m)
        elseif c2 > 0
            # Bifurcation
            o1 = Int(topo.strahler_order[c1])
            o2 = Int(topo.strahler_order[c2])
            if o1 == o2
                topo.strahler_order[i] = Int32(o1 + 1)
            else
                topo.strahler_order[i] = Int32(max(o1, o2))
            end
        else
            # Single child — inherit order
            topo.strahler_order[i] = topo.strahler_order[c1]
        end
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

"""
    build_empirical_connectivity(tree::VascularTree, params::MorphometricParams) -> Matrix{Float64}

Build empirical connectivity matrix from tree.
CM[daughter_order+1, parent_order+1] = count of order-(daughter_order) daughters
per order-(parent_order) parent. Strahler orders must already be assigned.
"""
function build_empirical_connectivity(tree::VascularTree, params::MorphometricParams)
    topo = tree.topology
    seg = tree.segments
    n = seg.n
    n_orders = params.n_orders

    CM = zeros(Float64, n_orders, n_orders)

    # sequential: scan all bifurcations and record parent-daughter order pairs
    for i in 1:n
        if topo.junction_type[i] == :bifurcation || topo.junction_type[i] == :trifurcation
            p_order = Int(topo.strahler_order[i])
            p_order < 0 && continue

            children = get_children(topo, Int32(i))
            for cid in children
                d_order = Int(topo.strahler_order[cid])
                d_order < 0 && continue
                if d_order < n_orders && p_order < n_orders
                    CM[d_order + 1, p_order + 1] += 1.0
                end
            end
        end
    end

    return CM
end

"""
    validate_connectivity(empirical::Matrix{Float64}, reference::Matrix{Float64}) -> (chi2, p_value)

Chi-squared test comparing empirical connectivity matrix to Kassab reference.
Returns (chi2_statistic, p_value). Only counts cells where reference > 0.
"""
function validate_connectivity(empirical::Matrix{Float64}, reference::Matrix{Float64})
    chi2 = 0.0
    dof = 0

    for j in axes(reference, 2)
        # Normalize column to get probabilities
        ref_col_sum = sum(reference[:, j])
        emp_col_sum = sum(empirical[:, j])
        ref_col_sum <= 0 && continue
        emp_col_sum <= 0 && continue

        for i in axes(reference, 1)
            expected = reference[i, j] / ref_col_sum * emp_col_sum
            expected < 0.5 && continue  # skip sparse cells

            observed = empirical[i, j]
            chi2 += (observed - expected)^2 / expected
            dof += 1
        end
    end

    dof = max(dof - 1, 1)

    # Approximate p-value using chi-squared distribution
    d = Chisq(dof)
    p_val = 1.0 - cdf(d, chi2)

    return (chi2, p_val)
end

"""
    constrain_connectivity!(tree::VascularTree, params::MorphometricParams)

Post-hoc adjustment: for each bifurcation, check if the parent-daughter order
combination is valid (has nonzero entry in Kassab connectivity matrix).
If invalid, adjust daughter radius to the nearest valid order.
Then propagate via Murray's law.
"""
function constrain_connectivity!(tree::VascularTree, params::MorphometricParams)
    topo = tree.topology
    seg = tree.segments
    n = seg.n
    gamma = params.gamma
    n_orders = params.n_orders

    # First ensure orders are assigned
    assign_strahler_orders!(tree, params)

    CM = params.connectivity_matrix
    changed = false

    # sequential: scan all bifurcations
    for i in 1:n
        if topo.junction_type[i] != :bifurcation
            continue
        end

        p_order = Int(topo.strahler_order[i])
        p_order < 1 && continue  # order 0 parents have no children in connectivity matrix

        children = get_children(topo, Int32(i))
        for cid in children
            d_order = Int(topo.strahler_order[cid])
            d_order < 0 && continue

            # Check if this combination has nonzero probability in Kassab matrix
            if p_order < n_orders && d_order < n_orders
                if CM[d_order + 1, p_order + 1] <= 0.0
                    # Invalid combination — find nearest valid daughter order
                    best_order = _find_nearest_valid_order(CM, p_order, d_order, n_orders)
                    if best_order >= 0 && best_order != d_order
                        # Set daughter radius to mean diameter of target order
                        target_diam_um = params.diameter_mean[best_order + 1]
                        seg.radius[cid] = target_diam_um / 2.0 / 1000.0  # um → mm
                        changed = true
                    end
                end
            end
        end
    end

    # Re-propagate Murray's law if any changes were made
    if changed
        update_radii!(tree, gamma)
    end

    return nothing
end

"""
    _find_nearest_valid_order(CM, parent_order, current_order, n_orders) -> Int

Find the nearest valid daughter order for a given parent order.
"""
function _find_nearest_valid_order(CM::Matrix{Float64}, parent_order::Int, current_order::Int, n_orders::Int)
    p_col = parent_order + 1
    p_col > n_orders && return current_order

    # Search outward from current order
    for delta in 0:(n_orders-1)
        for candidate in [current_order - delta, current_order + delta]
            if 0 <= candidate < n_orders
                if CM[candidate + 1, p_col] > 0.0
                    return candidate
                end
            end
        end
    end
    return current_order
end

"""
    apply_full_kassab_radii!(tree, params; rng)

Post-hoc Kassab radius refinement. Top-down asymmetry assignment with floor
enforcement, followed by bottom-up Murray's law propagation.

Unlike apply_kassab_radii!, this:
- Handles trifurcations (3 daughters)
- Enforces minimum radius floor (vessel_cutoff_um / 2 / 1000)
- Uses subtree size to assign larger radius to continuation
"""
function apply_full_kassab_radii!(
    tree::VascularTree,
    params::MorphometricParams;
    rng::AbstractRNG=Random.default_rng(),
)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    gamma = params.gamma
    min_radius = params.vessel_cutoff_um / 2.0 / 1000.0  # 8um diameter → 0.004mm radius

    # Compute subtree sizes for each segment (sequential: BFS bottom-up)
    order = _topo_order(tree)
    subtree_size = ones(Int, n)
    for k in length(order):-1:1
        i = order[k]
        c1 = Int(topo.child1_id[i])
        c2 = Int(topo.child2_id[i])
        c3 = Int(topo.child3_id[i])
        c1 > 0 && (subtree_size[i] += subtree_size[c1])
        c2 > 0 && (subtree_size[i] += subtree_size[c2])
        c3 > 0 && (subtree_size[i] += subtree_size[c3])
    end

    # Top-down pass: assign asymmetric daughter radii
    for k in 1:length(order)
        i = order[k]

        if topo.junction_type[i] == :bifurcation
            c1 = Int(topo.child1_id[i])
            c2 = Int(topo.child2_id[i])
            (c1 <= 0 || c2 <= 0) && continue

            asymmetry = sample_asymmetry(params, rng)
            r_large, r_small = compute_daughter_radii(seg.radius[i], asymmetry, gamma)

            # Floor enforcement
            r_large = max(r_large, min_radius)
            r_small = max(r_small, min_radius)

            # Assign larger radius to larger subtree (continuation)
            if subtree_size[c1] >= subtree_size[c2]
                seg.radius[c1] = r_large
                seg.radius[c2] = r_small
            else
                seg.radius[c1] = r_small
                seg.radius[c2] = r_large
            end

        elseif topo.junction_type[i] == :trifurcation
            c1 = Int(topo.child1_id[i])
            c2 = Int(topo.child2_id[i])
            c3 = Int(topo.child3_id[i])
            (c1 <= 0 || c2 <= 0 || c3 <= 0) && continue

            # Sample two asymmetry ratios for 3 daughters
            # Split: r_parent^gamma = r1^gamma + r2^gamma + r3^gamma
            # First split: r_parent into (r_large_pair, r3)
            asym1 = sample_asymmetry(params, rng)
            r_large_pair, r3 = compute_daughter_radii(seg.radius[i], asym1, gamma)

            # Second split: r_large_pair into (r1, r2)
            asym2 = sample_asymmetry(params, rng)
            r1, r2 = compute_daughter_radii(r_large_pair, asym2, gamma)

            # Floor enforcement
            r1 = max(r1, min_radius)
            r2 = max(r2, min_radius)
            r3 = max(r3, min_radius)

            # Sort by subtree size: largest subtree gets largest radius
            children = [(c1, subtree_size[c1]), (c2, subtree_size[c2]), (c3, subtree_size[c3])]
            sort!(children, by=x -> x[2], rev=true)
            radii = sort([r1, r2, r3], rev=true)

            seg.radius[children[1][1]] = radii[1]
            seg.radius[children[2][1]] = radii[2]
            seg.radius[children[3][1]] = radii[3]
        end
    end

    # Bottom-up: enforce Murray's law from terminals up
    update_radii!(tree, gamma)

    # Final floor enforcement pass (in case Murray propagation pushed some below floor)
    for i in 1:n
        if seg.radius[i] < min_radius
            seg.radius[i] = min_radius
        end
    end

    return nothing
end
