# Validation framework: statistical comparison against Kassab data

"""
    ValidationReport

Holds all validation metrics for a generated vascular tree.
"""
Base.@kwdef mutable struct ValidationReport
    n_segments::Int = 0
    n_bifurcations::Int = 0
    n_trifurcations::Int = 0

    # Diameter KS test p-values by Strahler order
    diameter_ks_pvalues::Dict{Int, Float64} = Dict{Int, Float64}()

    # Asymmetry ratios
    asymmetry_ratios::Vector{Float64} = Float64[]
    asymmetry_median::Float64 = 0.0

    # Length/Diameter ratio
    ld_ratios::Vector{Float64} = Float64[]
    ld_ratio_mean::Float64 = 0.0
    ld_ratio_std::Float64 = 0.0

    # Branching angles
    branching_angles::Vector{Float64} = Float64[]
    angle_mean::Float64 = 0.0
    angle_std::Float64 = 0.0

    # Extra metrics (connectivity chi2, trifurcation %, Murray deviation, lambda)
    extra::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

"""
    compute_diameters!(buf, segments)

Compute diameters (2*radius) for all active segments using AK.
"""
function compute_diameters!(buf::Vector{Float64}, seg::SegmentData)
    n = seg.n
    n == 0 && return nothing
    rad_view = @view seg.radius[1:n]
    dia_view = @view buf[1:n]
    AK.foreachindex(dia_view) do i
        dia_view[i] = 2.0 * rad_view[i]
    end
    return nothing
end

"""
    compute_ld_ratios!(buf, segments)

Compute length/diameter ratios for all active segments using AK.
"""
function compute_ld_ratios!(buf::Vector{Float64}, seg::SegmentData)
    n = seg.n
    n == 0 && return nothing
    len_view = @view seg.seg_length[1:n]
    rad_view = @view seg.radius[1:n]
    ld_view = @view buf[1:n]
    AK.foreachindex(ld_view) do i
        d = 2.0 * rad_view[i]
        ld_view[i] = d > 0.0 ? len_view[i] / d : 0.0
    end
    return nothing
end

"""
    compute_asymmetry_ratios(tree) -> Vector{Float64}

Compute asymmetry ratios (r_small / r_large) at all bifurcations.
"""
function compute_asymmetry_ratios(tree::VascularTree)
    topo = tree.topology
    seg = tree.segments
    n = seg.n
    ratios = Float64[]

    for i in 1:n
        topo.junction_type[i] != :bifurcation && continue
        c1 = Int(topo.child1_id[i])
        c2 = Int(topo.child2_id[i])
        (c1 <= 0 || c2 <= 0) && continue

        r1 = seg.radius[c1]
        r2 = seg.radius[c2]
        r_large = max(r1, r2)
        r_small = min(r1, r2)
        r_large > 0.0 && push!(ratios, r_small / r_large)
    end

    return ratios
end

"""
    compute_branching_angles(tree) -> Vector{Float64}

Compute the angle between daughter segments at each bifurcation.
"""
function compute_branching_angles(tree::VascularTree)
    topo = tree.topology
    seg = tree.segments
    n = seg.n
    angles = Float64[]

    for i in 1:n
        topo.junction_type[i] != :bifurcation && continue
        c1 = Int(topo.child1_id[i])
        c2 = Int(topo.child2_id[i])
        (c1 <= 0 || c2 <= 0) && continue

        # Direction vectors of daughters
        d1x = seg.distal_x[c1] - seg.proximal_x[c1]
        d1y = seg.distal_y[c1] - seg.proximal_y[c1]
        d1z = seg.distal_z[c1] - seg.proximal_z[c1]
        d2x = seg.distal_x[c2] - seg.proximal_x[c2]
        d2y = seg.distal_y[c2] - seg.proximal_y[c2]
        d2z = seg.distal_z[c2] - seg.proximal_z[c2]

        len1 = sqrt(d1x^2 + d1y^2 + d1z^2)
        len2 = sqrt(d2x^2 + d2y^2 + d2z^2)
        (len1 <= 0.0 || len2 <= 0.0) && continue

        dot_val = (d1x * d2x + d1y * d2y + d1z * d2z) / (len1 * len2)
        dot_val = clamp(dot_val, -1.0, 1.0)
        push!(angles, acos(dot_val))
    end

    return angles
end

"""
    compute_murray_deviation(tree, gamma) -> (max_dev, mean_dev)

Compute Murray's law deviation at all bifurcations and trifurcations.
Deviation = |r_parent^gamma - sum(r_child^gamma)| / r_parent^gamma.
"""
function compute_murray_deviation(tree::VascularTree, gamma::Float64)
    topo = tree.topology
    seg = tree.segments
    n = seg.n
    deviations = Float64[]

    for i in 1:n
        (topo.junction_type[i] != :bifurcation && topo.junction_type[i] != :trifurcation) && continue

        r_parent_g = seg.radius[i]^gamma
        r_parent_g <= 0.0 && continue

        r_sum = 0.0
        c1 = Int(topo.child1_id[i])
        c2 = Int(topo.child2_id[i])
        c3 = Int(topo.child3_id[i])
        c1 > 0 && (r_sum += seg.radius[c1]^gamma)
        c2 > 0 && (r_sum += seg.radius[c2]^gamma)
        c3 > 0 && (r_sum += seg.radius[c3]^gamma)

        dev = abs(r_parent_g - r_sum) / r_parent_g
        push!(deviations, dev)
    end

    isempty(deviations) && return (0.0, 0.0)
    return (maximum(deviations), sum(deviations) / length(deviations))
end

"""
    compute_trifurcation_pct(tree) -> Float64

Compute trifurcation percentage: n_trifurcations / (n_bifurcations + n_trifurcations) * 100.
"""
function compute_trifurcation_pct(tree::VascularTree)
    total = tree.n_bifurcations + tree.n_trifurcations
    total == 0 && return 0.0
    return tree.n_trifurcations / total * 100.0
end

"""
    compute_lambda_distribution(tree) -> Vector{Float64}

Compute lambda = l / w for each junction, where:
- l = minimum daughter segment length (distance between junction and nearest downstream node)
- w = 2 * pi * r_parent (parent circumference)

P(lambda -> 0) > 0 indicates stable trifurcations exist.
"""
function compute_lambda_distribution(tree::VascularTree)
    topo = tree.topology
    seg = tree.segments
    n = seg.n
    lambdas = Float64[]

    for i in 1:n
        (topo.junction_type[i] != :bifurcation && topo.junction_type[i] != :trifurcation) && continue

        w = 2π * seg.radius[i]
        w <= 0.0 && continue

        # Minimum daughter length
        min_len = Inf
        c1 = Int(topo.child1_id[i])
        c2 = Int(topo.child2_id[i])
        c3 = Int(topo.child3_id[i])
        c1 > 0 && (min_len = min(min_len, seg.seg_length[c1]))
        c2 > 0 && (min_len = min(min_len, seg.seg_length[c2]))
        c3 > 0 && (min_len = min(min_len, seg.seg_length[c3]))

        min_len < Inf && push!(lambdas, min_len / w)
    end

    return lambdas
end

"""
    validate_tree(tree, params) -> ValidationReport

Run all validation metrics on a generated tree.
"""
function validate_tree(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    report = ValidationReport()

    report.n_segments = n
    report.n_bifurcations = tree.n_bifurcations
    report.n_trifurcations = tree.n_trifurcations

    n == 0 && return report

    # --- Diameters by order + KS tests ---
    diameters_buf = zeros(n)
    compute_diameters!(diameters_buf, seg)

    # Assign Strahler orders
    assign_strahler_orders!(tree, params)

    for ord in 0:(params.n_orders - 1)
        diams_in_order = Float64[]
        for i in 1:n
            if topo.strahler_order[i] == ord
                push!(diams_in_order, diameters_buf[i] * 1000.0)  # convert to um
            end
        end
        length(diams_in_order) < 5 && continue

        # Simple KS-like statistic: compare mean to Kassab mean
        mean_d = sum(diams_in_order) / length(diams_in_order)
        expected_mean = params.diameter_mean[ord + 1]
        expected_sd = params.diameter_sd[ord + 1]

        # Approximate p-value using z-test of mean
        if expected_sd > 0.0 && length(diams_in_order) > 1
            se = expected_sd / sqrt(length(diams_in_order))
            z = abs(mean_d - expected_mean) / se
            # Two-tailed p from normal: p ≈ 2*exp(-z^2/2) (rough approximation)
            p_approx = min(1.0, 2.0 * exp(-z^2 / 2.0))
            report.diameter_ks_pvalues[ord] = p_approx
        end
    end

    # --- Asymmetry ratios ---
    report.asymmetry_ratios = compute_asymmetry_ratios(tree)
    if !isempty(report.asymmetry_ratios)
        sorted_ratios = sort(report.asymmetry_ratios)
        report.asymmetry_median = sorted_ratios[div(length(sorted_ratios) + 1, 2)]
    end

    # --- L/D ratio ---
    ld_buf = zeros(n)
    compute_ld_ratios!(ld_buf, seg)
    report.ld_ratios = ld_buf[1:n]
    if n > 0
        report.ld_ratio_mean = sum(report.ld_ratios) / n
        if n > 1
            sq_sum = sum((x - report.ld_ratio_mean)^2 for x in report.ld_ratios)
            report.ld_ratio_std = sqrt(sq_sum / (n - 1))
        end
    end

    # --- Branching angles ---
    report.branching_angles = compute_branching_angles(tree)
    if !isempty(report.branching_angles)
        na = length(report.branching_angles)
        report.angle_mean = sum(report.branching_angles) / na
        if na > 1
            sq_sum = sum((a - report.angle_mean)^2 for a in report.branching_angles)
            report.angle_std = sqrt(sq_sum / (na - 1))
        end
    end

    # --- Murray's law deviation ---
    max_dev, mean_dev = compute_murray_deviation(tree, params.gamma)
    report.extra[:murray_max_deviation] = max_dev
    report.extra[:murray_mean_deviation] = mean_dev

    # --- Trifurcation percentage ---
    report.extra[:trifurcation_pct] = compute_trifurcation_pct(tree)

    # --- P(lambda) distribution ---
    lambdas = compute_lambda_distribution(tree)
    report.extra[:lambda_values] = lambdas

    # --- Connectivity chi-squared ---
    # (strahler orders already assigned above)
    cm_empirical = build_empirical_connectivity(tree, params)
    chi2, p_val = validate_connectivity(cm_empirical, params.connectivity_matrix)
    report.extra[:connectivity_chi2] = chi2
    report.extra[:connectivity_pval] = p_val

    return report
end

"""
    print_report([io], report)

Print a formatted summary of the validation report.
"""
function print_report(io::IO, report::ValidationReport)
    println(io, "=== VesselTree Validation Report ===")
    println(io, "  Segments: $(report.n_segments)")
    println(io, "  Bifurcations: $(report.n_bifurcations)")
    println(io, "  Trifurcations: $(report.n_trifurcations)")
    println(io, "")

    if !isempty(report.diameter_ks_pvalues)
        println(io, "  Diameter KS p-values by order:")
        for ord in sort(collect(keys(report.diameter_ks_pvalues)))
            p = report.diameter_ks_pvalues[ord]
            println(io, "    Order $ord: p = $(round(p, digits=4))")
        end
        println(io, "")
    end

    if !isempty(report.asymmetry_ratios)
        println(io, "  Asymmetry ratio:")
        println(io, "    N = $(length(report.asymmetry_ratios))")
        println(io, "    Median = $(round(report.asymmetry_median, digits=4))")
        println(io, "")
    end

    println(io, "  L/D ratio:")
    println(io, "    Mean = $(round(report.ld_ratio_mean, digits=2))")
    println(io, "    Std  = $(round(report.ld_ratio_std, digits=2))")
    println(io, "")

    if !isempty(report.branching_angles)
        println(io, "  Branching angles (rad):")
        println(io, "    N = $(length(report.branching_angles))")
        println(io, "    Mean = $(round(report.angle_mean, digits=4))")
        println(io, "    Std  = $(round(report.angle_std, digits=4))")
        println(io, "")
    end

    if haskey(report.extra, :murray_max_deviation)
        println(io, "  Murray's law deviation:")
        println(io, "    Max  = $(report.extra[:murray_max_deviation])")
        println(io, "    Mean = $(report.extra[:murray_mean_deviation])")
        println(io, "")
    end

    if haskey(report.extra, :trifurcation_pct)
        println(io, "  Trifurcation percentage: $(round(report.extra[:trifurcation_pct], digits=2))%")
        println(io, "")
    end

    if haskey(report.extra, :lambda_values) && !isempty(report.extra[:lambda_values])
        lambdas = report.extra[:lambda_values]
        println(io, "  P(lambda) distribution:")
        println(io, "    N = $(length(lambdas))")
        mean_l = sum(lambdas) / length(lambdas)
        println(io, "    Mean lambda = $(round(mean_l, digits=4))")
        println(io, "    Min lambda  = $(round(minimum(lambdas), digits=4))")
        println(io, "")
    end

    if haskey(report.extra, :connectivity_chi2)
        println(io, "  Connectivity matrix:")
        println(io, "    Chi-squared = $(round(report.extra[:connectivity_chi2], digits=2))")
        println(io, "    p-value     = $(round(report.extra[:connectivity_pval], digits=4))")
    end
end

print_report(report::ValidationReport) = print_report(stdout, report)
