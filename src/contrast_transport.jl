struct ContrastTransportResult{T<:AbstractFloat}
    times::Vector{Float64}
    concentration::Matrix{T}        # segment-average concentration [mg/mL]
    outlet_concentration::Union{Nothing, Matrix{T}}
    transit_time_s::Vector{Float64}
    segment_volume_mL::Vector{Float64}
end

"""
    gamma_variate_input(times; amplitude, t0, tmax, alpha)

Gamma-variate input curve used in the earlier brain-phantom workflow.
Returns concentration in the same units as `amplitude`.
"""
function gamma_variate_input(
    times::AbstractVector{<:Real};
    amplitude::Real=1.0,
    t0::Real=0.0,
    tmax::Real,
    alpha::Real,
)
    t = Float64.(times)
    y = zeros(Float64, length(t))
    denom = Float64(tmax - t0)
    denom > 0.0 || error("Expected tmax > t0, got t0=$(t0), tmax=$(tmax).")

    for i in eachindex(t)
        if t[i] > t0
            tprime = (t[i] - t0) / denom
            y[i] = Float64(amplitude) * (tprime^Float64(alpha)) * exp(Float64(alpha) * (1.0 - tprime))
        end
    end
    return y
end

"""
    root_pulse_input(times, concentration; duration=dt)

Construct a short inlet pulse. This is the closest direct analogue to
"root concentration at t=0" when a full inlet time-course is not supplied.
"""
function root_pulse_input(
    times::AbstractVector{<:Real},
    concentration::Real;
    duration::Union{Nothing, Real}=nothing,
)
    t = Float64.(times)
    length(t) >= 2 || error("Need at least two time points to build a pulse input.")
    dt = t[2] - t[1]
    dt > 0.0 || error("Times must be strictly increasing.")
    pulse_duration = isnothing(duration) ? dt : Float64(duration)

    curve = zeros(Float64, length(t))
    for i in eachindex(t)
        if t[i] <= t[1] + pulse_duration + 1e-12
            curve[i] = Float64(concentration)
        end
    end
    return curve
end

function _uniform_dt(times::Vector{Float64}; rtol::Float64=1e-6, atol::Float64=1e-9)
    length(times) >= 2 || error("Need at least two time points.")
    dt = times[2] - times[1]
    dt > 0.0 || error("Times must be strictly increasing.")
    for i in 3:length(times)
        next_dt = times[i] - times[i - 1]
        isapprox(next_dt, dt; rtol=rtol, atol=atol) || error("Contrast transport currently requires a uniform time step.")
    end
    return dt
end

function _resolve_root_input_curve(times::Vector{Float64}, root_input)
    if root_input isa Number
        return root_pulse_input(times, Float64(root_input))
    elseif root_input isa AbstractVector
        length(root_input) == length(times) || error("Root input length $(length(root_input)) does not match times length $(length(times)).")
        return Float64.(root_input)
    else
        error("Unsupported root_input type $(typeof(root_input)). Use a scalar pulse concentration or a vector time-course.")
    end
end

@inline function _sample_uniform_history(
    history::AbstractVector{T},
    t_query::Float64,
    t0::Float64,
    dt::Float64,
    current_idx::Int,
) where {T<:AbstractFloat}
    if t_query < t0
        return 0.0
    end

    pos = (t_query - t0) / dt + 1.0
    lo = floor(Int, pos)
    hi = ceil(Int, pos)
    hi <= current_idx || return Float64(history[current_idx])
    lo = max(lo, 1)
    hi = max(hi, 1)

    if lo == hi
        return Float64(history[lo])
    end

    w = pos - lo
    return (1.0 - w) * Float64(history[lo]) + w * Float64(history[hi])
end

function _segment_volume_mL(tree::VascularTree, seg_id::Int)
    seg = tree.segments
    r_mm = seg.radius[seg_id]
    l_mm = seg.seg_length[seg_id]
    volume_mm3 = π * r_mm^2 * l_mm
    return volume_mm3 / 1000.0
end

function _segment_transit_times_s(
    tree::VascularTree;
    min_flow_m3_s::Float64=1e-18,
    min_transit_time_s::Float64=1e-9,
)
    seg = tree.segments
    n = seg.n
    τ = zeros(Float64, n)
    volumes = zeros(Float64, n)
    for i in 1:n
        vol_mL = _segment_volume_mL(tree, i)
        volumes[i] = vol_mL
        flow_m3_s = seg.flow[i]
        abs(flow_m3_s) > min_flow_m3_s || error("Segment $i has near-zero flow $(flow_m3_s). Compute hemodynamics before contrast transport.")
        vol_m3 = vol_mL * 1e-6
        τ[i] = max(vol_m3 / abs(flow_m3_s), min_transit_time_s)
    end
    return τ, volumes
end

function _prepare_contrast_hemodynamics!(
    tree::VascularTree,
    params::Union{Nothing, MorphometricParams},
    recompute_hemodynamics::Bool,
)
    if recompute_hemodynamics
        isnothing(params) && error("MorphometricParams are required when recompute_hemodynamics=true.")
        compute_resistances!(tree, params.blood_viscosity)
        compute_flows!(tree, params)
        compute_pressures!(tree, params)
    end
    return tree
end

"""
    simulate_contrast_transport(tree, times; root_input, params, recompute_hemodynamics=true,
                                storage_type=Float32, return_outlet=false)

Simulate iodine transport through a directed arterial tree.

Model:
- each segment is a well-mixed control volume,
- inlet concentration is inherited from the parent outlet concentration,
- outlet concentration is the inlet concentration delayed by the segment transit time,
- segment-average concentration follows a mass-balance ODE:
  `dC/dt = (Q/V) * (Cin - Cout)`.

Concentrations are reported in `mg/mL` if the root input is given in `mg/mL`.
"""
function simulate_contrast_transport(
    tree::VascularTree,
    times::AbstractVector{<:Real};
    root_input,
    params::Union{Nothing, MorphometricParams}=nothing,
    recompute_hemodynamics::Bool=params !== nothing,
    storage_type::Type{T}=Float32,
    return_outlet::Bool=false,
    min_flow_m3_s::Float64=1e-18,
    min_transit_time_s::Float64=1e-9,
) where {T<:AbstractFloat}
    work_tree = deepcopy(tree)
    _prepare_contrast_hemodynamics!(work_tree, params, recompute_hemodynamics)

    t = Float64.(times)
    dt = _uniform_dt(t)
    root_curve = _resolve_root_input_curve(t, root_input)

    seg = work_tree.segments
    topo = work_tree.topology
    n = seg.n
    nt = length(t)
    order = _topo_order(work_tree)
    τ, volumes = _segment_transit_times_s(
        work_tree;
        min_flow_m3_s=min_flow_m3_s,
        min_transit_time_s=min_transit_time_s,
    )

    conc = zeros(T, n, nt)
    out = zeros(T, n, nt)
    t0 = t[1]

    for ti in 1:nt
        t_now = t[ti]
        for seg_id in order
            parent = Int(topo.parent_id[seg_id])

            cin_now =
                if parent > 0
                    Float64(out[parent, ti])
                else
                    root_curve[ti]
                end

            delayed_t = t_now - τ[seg_id]
            cout_now =
                if parent > 0
                    _sample_uniform_history(view(out, parent, :), delayed_t, t0, dt, ti)
                else
                    _sample_uniform_history(root_curve, delayed_t, t0, dt, ti)
                end

            out[seg_id, ti] = T(cout_now)

            if ti == 1
                conc[seg_id, ti] = zero(T)
            else
                cin_prev =
                    if parent > 0
                        Float64(out[parent, ti - 1])
                    else
                        root_curve[ti - 1]
                    end
                cout_prev = Float64(out[seg_id, ti - 1])
                dcdt = (cin_prev - cout_prev) / τ[seg_id]
                conc[seg_id, ti] = T(max(Float64(conc[seg_id, ti - 1]) + dt * dcdt, 0.0))
            end
        end
    end

    return ContrastTransportResult(
        t,
        conc,
        return_outlet ? out : nothing,
        τ,
        volumes,
    )
end

function segment_mass_mg(result::ContrastTransportResult)
    n, nt = size(result.concentration)
    mass = similar(result.concentration)
    for i in 1:n
        mass[i, :] .= result.concentration[i, :] .* result.segment_volume_mL[i]
    end
    return mass
end

function _forest_root_input(root_inputs, name::String)
    haskey(root_inputs, name) || error("Missing root input for tree `$name`.")
    return root_inputs[name]
end

"""
    simulate_forest_contrast(forest, times; root_inputs, recompute_hemodynamics=true, ...)

Run contrast transport independently on each tree in a `CoronaryForest`.
`root_inputs` maps tree name to either:
- a scalar concentration, interpreted as a short inlet pulse at `t=0`, or
- a full concentration time-course vector aligned with `times`.
"""
function simulate_forest_contrast(
    forest::CoronaryForest,
    times::AbstractVector{<:Real};
    root_inputs,
    recompute_hemodynamics::Bool=true,
    storage_type::Type{T}=Float32,
    return_outlet::Bool=false,
    threaded::Bool=false,
    min_flow_m3_s::Float64=1e-18,
    min_transit_time_s::Float64=1e-9,
) where {T<:AbstractFloat}
    names = sort!(collect(keys(forest.trees)))
    per_tree_results = Vector{ContrastTransportResult{T}}(undef, length(names))

    if threaded && length(names) > 1
        Threads.@threads for idx in eachindex(names)
            name = names[idx]
            tree = forest.trees[name]
            params = get(forest.tree_params, name, forest.params)
            per_tree_results[idx] = simulate_contrast_transport(
                tree,
                times;
                root_input=_forest_root_input(root_inputs, name),
                params=params,
                recompute_hemodynamics=recompute_hemodynamics,
                storage_type=storage_type,
                return_outlet=return_outlet,
                min_flow_m3_s=min_flow_m3_s,
                min_transit_time_s=min_transit_time_s,
            )
        end
    else
        for idx in eachindex(names)
            name = names[idx]
            tree = forest.trees[name]
            params = get(forest.tree_params, name, forest.params)
            per_tree_results[idx] = simulate_contrast_transport(
                tree,
                times;
                root_input=_forest_root_input(root_inputs, name),
                params=params,
                recompute_hemodynamics=recompute_hemodynamics,
                storage_type=storage_type,
                return_outlet=return_outlet,
                min_flow_m3_s=min_flow_m3_s,
                min_transit_time_s=min_transit_time_s,
            )
        end
    end

    return Dict(names[idx] => per_tree_results[idx] for idx in eachindex(names))
end
