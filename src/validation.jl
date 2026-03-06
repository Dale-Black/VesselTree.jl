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

        # Proper one-sample KS test against N(mean, sd) from Kassab tables
        mu = params.diameter_mean[ord + 1]
        sd = params.diameter_sd[ord + 1]
        if sd > 0.0
            d_normal = Normal(mu, sd)
            _, p = ks_test_onesample(diams_in_order, x -> cdf(d_normal, x))
            report.diameter_ks_pvalues[ord] = p
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

# --- Proper KS test + comprehensive Kassab validation ---

"""
    ks_test_onesample(samples, cdf_func) -> (D, p_value)

One-sample Kolmogorov-Smirnov test. Compares `samples` to a theoretical
distribution given by `cdf_func(x)`. Returns the KS statistic D and an
asymptotic p-value (accurate for n > 35).
"""
function ks_test_onesample(samples::Vector{Float64}, cdf_func)
    n = length(samples)
    n == 0 && return (0.0, 1.0)
    sorted = sort(samples)
    D = 0.0
    for i in 1:n
        F_emp = i / n
        F_emp_prev = (i - 1) / n
        F_theo = cdf_func(sorted[i])
        D = max(D, abs(F_emp - F_theo), abs(F_emp_prev - F_theo))
    end
    # Kolmogorov's asymptotic formula: p ≈ 2*sum_{k=1}^inf (-1)^(k+1) * exp(-2*k^2*n*D^2)
    # First term approximation (accurate for moderate D):
    lambda = (sqrt(n) + 0.12 + 0.11 / sqrt(n)) * D
    if lambda <= 0.0
        p = 1.0
    elseif lambda > 3.1
        p = 0.0
    else
        # Smirnov's series (first 4 terms)
        p = 0.0
        for k in 1:4
            p += (-1)^(k + 1) * exp(-2.0 * k^2 * lambda^2)
        end
        p = clamp(2.0 * p, 0.0, 1.0)
    end
    return (D, p)
end

"""
    validate_diameters_per_order(tree, params) -> Dict{Int, NamedTuple}

Proper KS test of diameter distribution per Strahler order against
N(mean_um, sd_um) from Kassab tables. Returns dict mapping order -> (n, D, p, mean_um).
"""
function validate_diameters_per_order(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    assign_strahler_orders!(tree, params)

    results = Dict{Int, NamedTuple{(:n, :D, :p, :mean_um), Tuple{Int, Float64, Float64, Float64}}}()
    for ord in 0:(params.n_orders - 1)
        diams_um = Float64[]
        for i in 1:n
            if topo.strahler_order[i] == ord
                push!(diams_um, seg.radius[i] * 2.0 * 1000.0)
            end
        end
        length(diams_um) < 5 && continue
        mu = params.diameter_mean[ord + 1]
        sd = params.diameter_sd[ord + 1]
        sd <= 0.0 && continue
        d_normal = Normal(mu, sd)
        D, p = ks_test_onesample(diams_um, x -> cdf(d_normal, x))
        mean_um = sum(diams_um) / length(diams_um)
        results[ord] = (n=length(diams_um), D=D, p=p, mean_um=mean_um)
    end
    return results
end

"""
    validate_lengths_per_order(tree, params) -> Dict{Int, NamedTuple}

Proper KS test of length distribution per Strahler order against
N(mean_um, sd_um) from Kassab tables.
"""
function validate_lengths_per_order(tree::VascularTree, params::MorphometricParams)
    seg = tree.segments
    topo = tree.topology
    n = seg.n
    assign_strahler_orders!(tree, params)

    results = Dict{Int, NamedTuple{(:n, :D, :p, :mean_um), Tuple{Int, Float64, Float64, Float64}}}()
    for ord in 0:(params.n_orders - 1)
        lengths_um = Float64[]
        for i in 1:n
            if topo.strahler_order[i] == ord
                push!(lengths_um, seg.seg_length[i] * 1000.0)  # mm -> um
            end
        end
        length(lengths_um) < 5 && continue
        mu = params.length_mean[ord + 1]
        sd = params.length_sd[ord + 1]
        sd <= 0.0 && continue
        d_normal = Normal(mu, sd)
        D, p = ks_test_onesample(lengths_um, x -> cdf(d_normal, x))
        mean_um = sum(lengths_um) / length(lengths_um)
        results[ord] = (n=length(lengths_um), D=D, p=p, mean_um=mean_um)
    end
    return results
end

"""
    validate_segment_counts(tree, params) -> Dict{Int, Int}

Count segments per Strahler order.
"""
function validate_segment_counts(tree::VascularTree, params::MorphometricParams)
    topo = tree.topology
    n = tree.segments.n
    assign_strahler_orders!(tree, params)
    counts = Dict{Int, Int}()
    for i in 1:n
        ord = Int(topo.strahler_order[i])
        ord < 0 && continue
        counts[ord] = get(counts, ord, 0) + 1
    end
    return counts
end

"""
    validate_asymmetry_ks(tree, params) -> (n, D, p, median)

KS test of asymmetry ratios against Beta(alpha, beta) from params.
"""
function validate_asymmetry_ks(tree::VascularTree, params::MorphometricParams)
    ratios = compute_asymmetry_ratios(tree)
    isempty(ratios) && return (n=0, D=0.0, p=1.0, median=0.0)
    d_beta = Beta(params.asymmetry_alpha, params.asymmetry_beta)
    D, p = ks_test_onesample(ratios, x -> cdf(d_beta, x))
    sorted = sort(ratios)
    med = sorted[div(length(sorted) + 1, 2)]
    return (n=length(ratios), D=D, p=p, median=med)
end

"""
    validate_connectivity_matrix(params) -> Dict{Symbol, Any}

Audit the connectivity matrix: column sums, diagonal dominance, monotonicity.
"""
function validate_connectivity_matrix(params::MorphometricParams)
    CM = params.connectivity_matrix
    n_orders = params.n_orders
    result = Dict{Symbol, Any}()

    col_sums = Float64[]
    for j in 2:n_orders  # parent orders 1-11
        s = sum(CM[:, j])
        push!(col_sums, s)
    end
    result[:column_sums] = col_sums
    result[:all_sums_ge_2] = all(s >= 1.9 for s in col_sums)
    result[:all_sums_le_3] = all(s <= 3.1 for s in col_sums)

    # Diagonal dominance: CM[n, n+1] (order n-1 daughter for order n parent)
    # should be the largest entry in each column
    diag_dominant = true
    for j in 2:n_orders
        diag_val = CM[j - 1, j]  # CM[n-1+1, n+1] = CM[j-1, j]
        # Actually: parent order j-1, daughter order j-2 → CM[j-1, j]
        col_max = maximum(CM[:, j])
        if diag_val < col_max - 1e-10
            diag_dominant = false
        end
    end
    result[:diagonal_dominant] = diag_dominant

    return result
end

"""
    generate_report_card(tree, params) -> Dict{Symbol, Any}

Comprehensive pass/fail report card for a generated tree.
Returns a dict with individual metric results and an overall grade.
"""
function generate_report_card(tree::VascularTree, params::MorphometricParams)
    card = Dict{Symbol, Any}()
    passed = 0
    total = 0

    # 1. Per-order diameter KS tests
    diam_results = validate_diameters_per_order(tree, params)
    diam_pass = 0
    diam_total = length(diam_results)
    for (ord, r) in diam_results
        r.p > 0.05 && (diam_pass += 1)
    end
    card[:diameter_ks] = diam_results
    card[:diameter_ks_pass] = diam_pass
    card[:diameter_ks_total] = diam_total
    if diam_total > 0
        total += 1
        diam_pass >= max(1, ceil(Int, diam_total * 0.5)) && (passed += 1)
    end

    # 2. Per-order length KS tests
    len_results = validate_lengths_per_order(tree, params)
    len_pass = 0
    len_total = length(len_results)
    for (ord, r) in len_results
        r.p > 0.05 && (len_pass += 1)
    end
    card[:length_ks] = len_results
    card[:length_ks_pass] = len_pass
    card[:length_ks_total] = len_total
    if len_total > 0
        total += 1
        len_pass >= max(1, ceil(Int, len_total * 0.5)) && (passed += 1)
    end

    # 3. Connectivity chi-squared
    assign_strahler_orders!(tree, params)
    cm_emp = build_empirical_connectivity(tree, params)
    chi2, p_conn = validate_connectivity(cm_emp, params.connectivity_matrix)
    card[:connectivity_chi2] = chi2
    card[:connectivity_pval] = p_conn
    total += 1
    p_conn > 0.01 && (passed += 1)

    # 4. Asymmetry KS
    asym = validate_asymmetry_ks(tree, params)
    card[:asymmetry] = asym
    if asym.n >= 5
        total += 1
        asym.p > 0.05 && (passed += 1)
    end

    # 5. Murray's law
    max_dev, mean_dev = compute_murray_deviation(tree, params.gamma)
    card[:murray_max_dev] = max_dev
    card[:murray_mean_dev] = mean_dev
    total += 1
    mean_dev < 0.01 && (passed += 1)

    # 6. Segment count per order
    counts = validate_segment_counts(tree, params)
    card[:segment_counts] = counts
    card[:n_orders_populated] = length(counts)
    total += 1
    length(counts) >= 2 && (passed += 1)

    # 7. Trifurcation percentage
    tri_pct = compute_trifurcation_pct(tree)
    card[:trifurcation_pct] = tri_pct

    # Overall grade
    card[:passed] = passed
    card[:total] = total
    card[:grade] = total > 0 ? passed / total : 0.0

    return card
end

"""
    print_report_card([io], card)

Print a formatted pass/fail report card.
"""
function print_report_card(io::IO, card::Dict{Symbol, Any})
    println(io, "=== Kassab Validation Report Card ===")
    println(io, "")

    # Diameter KS
    if haskey(card, :diameter_ks)
        dp = card[:diameter_ks_pass]; dt = card[:diameter_ks_total]
        status = dp >= max(1, ceil(Int, dt * 0.5)) ? "PASS" : "FAIL"
        println(io, "  [$status] Diameter KS: $dp/$dt orders pass (p > 0.05)")
        for ord in sort(collect(keys(card[:diameter_ks])))
            r = card[:diameter_ks][ord]
            mark = r.p > 0.05 ? "v" : "x"
            println(io, "    [$mark] Order $ord: n=$(r.n), D=$(round(r.D, digits=4)), p=$(round(r.p, digits=4)), mean=$(round(r.mean_um, digits=1))um")
        end
        println(io, "")
    end

    # Length KS
    if haskey(card, :length_ks)
        lp = card[:length_ks_pass]; lt = card[:length_ks_total]
        status = lp >= max(1, ceil(Int, lt * 0.5)) ? "PASS" : "FAIL"
        println(io, "  [$status] Length KS: $lp/$lt orders pass (p > 0.05)")
        println(io, "")
    end

    # Connectivity
    if haskey(card, :connectivity_pval)
        p = card[:connectivity_pval]
        status = p > 0.01 ? "PASS" : "FAIL"
        println(io, "  [$status] Connectivity chi2: $(round(card[:connectivity_chi2], digits=2)), p=$(round(p, digits=4))")
        println(io, "")
    end

    # Asymmetry
    if haskey(card, :asymmetry)
        a = card[:asymmetry]
        if a.n >= 5
            status = a.p > 0.05 ? "PASS" : "FAIL"
            println(io, "  [$status] Asymmetry KS: n=$(a.n), D=$(round(a.D, digits=4)), p=$(round(a.p, digits=4)), median=$(round(a.median, digits=4))")
        else
            println(io, "  [SKIP] Asymmetry KS: n=$(a.n) (need >= 5)")
        end
        println(io, "")
    end

    # Murray's law
    if haskey(card, :murray_mean_dev)
        md = card[:murray_mean_dev]
        status = md < 0.01 ? "PASS" : "FAIL"
        println(io, "  [$status] Murray deviation: mean=$(round(md, digits=6)), max=$(round(card[:murray_max_dev], digits=6))")
        println(io, "")
    end

    # Segment counts
    if haskey(card, :segment_counts)
        n_ord = card[:n_orders_populated]
        status = n_ord >= 2 ? "PASS" : "FAIL"
        println(io, "  [$status] Orders populated: $n_ord")
        for ord in sort(collect(keys(card[:segment_counts])))
            println(io, "    Order $ord: $(card[:segment_counts][ord]) segments")
        end
        println(io, "")
    end

    # Trifurcation
    if haskey(card, :trifurcation_pct)
        println(io, "  [INFO] Trifurcation: $(round(card[:trifurcation_pct], digits=2))%")
        println(io, "")
    end

    # Grade
    println(io, "  Overall: $(card[:passed])/$(card[:total]) metrics pass ($(round(card[:grade] * 100, digits=0))%)")
end

print_report_card(card::Dict{Symbol, Any}) = print_report_card(stdout, card)
