using Printf

struct UnifiedViewerPointLayerData
    key::String
    label::String
    x::Vector{Float32}
    y::Vector{Float32}
    z::Vector{Float32}
    size::Float32
    color::String
    opacity::Float32
    visible::Bool
end

struct UnifiedViewerLineLayerData
    key::String
    label::String
    x::Vector{Float32}
    y::Vector{Float32}
    z::Vector{Float32}
    width::Float32
    color::String
    opacity::Float32
    visible::Bool
end

struct UnifiedViewerData
    title::String
    times::Vector{Float32}
    iodine_trees::Vector{ContrastViewerTreeData}
    point_layers::Vector{UnifiedViewerPointLayerData}
    line_layers::Vector{UnifiedViewerLineLayerData}
    cmax::Float32
    blue_reference::Float32
end

function _subsample_rows_evenly(mat::Matrix{Float64}, max_points::Int)
    max_points <= 0 && return Matrix{Float64}(undef, 0, size(mat, 2))
    n = size(mat, 1)
    n <= max_points && return copy(mat)
    idx = unique(round.(Int, range(1, n, length=max_points)))
    return mat[idx, :]
end

function _f32_xyz(mat::Matrix{Float64})
    return Float32.(mat[:, 1]), Float32.(mat[:, 2]), Float32.(mat[:, 3])
end

function _prepare_domain_layers(domain::AbstractDomain; max_outer_points::Int=12000, max_cavity_points::Int=6000)
    point_layers = UnifiedViewerPointLayerData[]
    if domain isa CSVShellDomain
        outer = _subsample_rows_evenly(domain.outer_surface_points, max_outer_points)
        x, y, z = _f32_xyz(outer)
        push!(point_layers, UnifiedViewerPointLayerData(
            "domain",
            "Domain Surface",
            x, y, z,
            1.6f0,
            "#7aa6a1",
            0.18f0,
            false,
        ))
        if !isempty(domain.cavity_surface_points)
            mats = Matrix{Float64}[]
            for pts in domain.cavity_surface_points
                push!(mats, _subsample_rows_evenly(pts, max(1, ceil(Int, max_cavity_points / max(length(domain.cavity_surface_points), 1)))))
            end
            total_rows = sum(size(m, 1) for m in mats)
            cavity = Matrix{Float64}(undef, total_rows, 3)
            offset = 0
            for m in mats
                cavity[(offset + 1):(offset + size(m, 1)), :] .= m
                offset += size(m, 1)
            end
            x, y, z = _f32_xyz(cavity)
            push!(point_layers, UnifiedViewerPointLayerData(
                "cavities",
                "Cardiac Cavities",
                x, y, z,
                1.4f0,
                "#c9d7df",
                0.12f0,
                false,
            ))
        end
    elseif domain isa CSVVolumeDomain
        surface = _subsample_rows_evenly(domain.surface_points, max_outer_points)
        x, y, z = _f32_xyz(surface)
        push!(point_layers, UnifiedViewerPointLayerData(
            "domain",
            "Domain Surface",
            x, y, z,
            1.6f0,
            "#7aa6a1",
            0.18f0,
            false,
        ))
    end
    return point_layers
end

function _sample_surface_points(
    surfaces::AbstractVector{XCATNurbsSurface};
    n_u::Int=18,
    n_v::Int=10,
    max_points::Int=6000,
)
    mats = Matrix{Float64}[]
    for surface in surfaces
        pts, _ = xcat_sampled_surface_rows(surface; n_u=n_u, n_v=n_v, orient_outward=false)
        push!(mats, pts)
    end
    total_rows = sum(size(m, 1) for m in mats)
    total_rows == 0 && return Matrix{Float64}(undef, 0, 3)
    out = Matrix{Float64}(undef, total_rows, 3)
    offset = 0
    for m in mats
        out[(offset + 1):(offset + size(m, 1)), :] .= m
        offset += size(m, 1)
    end
    return _subsample_rows_evenly(out, max_points)
end

function _segment_ids_for_static_layer(tree::VascularTree, ids::Vector{Int}; max_segments::Int)
    isempty(ids) && return Int[]
    length(ids) <= max_segments && return sort(ids)
    radii = tree.segments.radius[ids]
    n_large = clamp(round(Int, 0.35 * max_segments), 1, max_segments)
    order = sortperm(radii; rev=true)
    chosen = ids[order[1:min(n_large, length(order))]]
    chosen_set = Set(chosen)
    remaining = [id for id in ids if !(id in chosen_set)]
    append!(chosen, _sample_evenly(remaining, max_segments - length(chosen)))
    chosen = unique(chosen)
    if length(chosen) > max_segments
        chosen = _sample_evenly(sort(chosen), max_segments)
    end
    sort!(chosen)
    return chosen
end

function _line_arrays_from_segments(tree::VascularTree, ids::AbstractVector{<:Integer})
    seg = tree.segments
    x = Float32[]
    y = Float32[]
    z = Float32[]
    for id in ids
        push!(x, Float32(seg.proximal_x[id]), Float32(seg.distal_x[id]), NaN32)
        push!(y, Float32(seg.proximal_y[id]), Float32(seg.distal_y[id]), NaN32)
        push!(z, Float32(seg.proximal_z[id]), Float32(seg.distal_z[id]), NaN32)
    end
    return x, y, z
end

function _combine_line_layers(entries::Vector{Tuple{VascularTree, Vector{Int}}})
    x = Float32[]
    y = Float32[]
    z = Float32[]
    for (tree, ids) in entries
        tx, ty, tz = _line_arrays_from_segments(tree, ids)
        append!(x, tx)
        append!(y, ty)
        append!(z, tz)
    end
    return x, y, z
end

function _prepare_fixed_trunk_layer(fixed_trees::Dict{String, VascularTree})
    isempty(fixed_trees) && return nothing
    entries = Tuple{VascularTree, Vector{Int}}[]
    for name in sort(collect(keys(fixed_trees)))
        tree = fixed_trees[name]
        push!(entries, (tree, collect(1:tree.segments.n)))
    end
    x, y, z = _combine_line_layers(entries)
    return UnifiedViewerLineLayerData(
        "fixed_trunks",
        "Imported Fixed Trunks",
        x, y, z,
        3.0f0,
        "#1f4e79",
        0.9f0,
        true,
    )
end

function _prepare_grown_tree_layer(
    forest::CoronaryForest,
    fixed_trees::Dict{String, VascularTree};
    max_segments_per_tree::Int=3500,
)
    entries = Tuple{VascularTree, Vector{Int}}[]
    for name in sort(collect(keys(forest.trees)))
        tree = forest.trees[name]
        n_fixed = haskey(fixed_trees, name) ? fixed_trees[name].segments.n : 0
        grown_ids = collect((n_fixed + 1):tree.segments.n)
        isempty(grown_ids) && continue
        chosen = _segment_ids_for_static_layer(tree, grown_ids; max_segments=max_segments_per_tree)
        push!(entries, (tree, chosen))
    end
    isempty(entries) && return nothing
    x, y, z = _combine_line_layers(entries)
    return UnifiedViewerLineLayerData(
        "grown_tree",
        "Grown Tree",
        x, y, z,
        1.2f0,
        "#6f6f73",
        0.45f0,
        true,
    )
end

function prepare_unified_viewer_data(
    forest::CoronaryForest,
    results::Dict{String, <:ContrastTransportResult};
    title::AbstractString="Unified Vascular Viewer",
    domain::Union{Nothing, AbstractDomain}=nothing,
    original_surfaces::AbstractVector{XCATNurbsSurface}=XCATNurbsSurface[],
    fixed_trees::Dict{String, VascularTree}=Dict{String, VascularTree}(),
    domain_max_points::Int=12000,
    cavity_max_points::Int=6000,
    surface_n_u::Int=18,
    surface_n_v::Int=10,
    original_surface_max_points::Int=6000,
    grown_max_segments_per_tree::Int=3500,
    iodine_min_radius_um::Float64=4.0,
    iodine_max_segments_per_tree::Int=3000,
    time_stride::Int=1,
    radius_scale::Float64=10.0,
    blue_reference_mg_mL::Float64=3.0,
)
    iodine = prepare_contrast_viewer_data(
        forest,
        results;
        min_radius_um=iodine_min_radius_um,
        max_segments_per_tree=iodine_max_segments_per_tree,
        time_stride=time_stride,
        radius_scale=radius_scale,
        blue_reference_mg_mL=blue_reference_mg_mL,
    )

    point_layers = UnifiedViewerPointLayerData[]
    if domain !== nothing
        append!(point_layers, _prepare_domain_layers(domain; max_outer_points=domain_max_points, max_cavity_points=cavity_max_points))
    end
    if !isempty(original_surfaces)
        pts = _sample_surface_points(original_surfaces; n_u=surface_n_u, n_v=surface_n_v, max_points=original_surface_max_points)
        if size(pts, 1) > 0
            x, y, z = _f32_xyz(pts)
            push!(point_layers, UnifiedViewerPointLayerData(
                "xcat_surfaces",
                "XCAT Vessel Surfaces",
                x, y, z,
                1.8f0,
                "#b36b00",
                0.28f0,
                false,
            ))
        end
    end

    line_layers = UnifiedViewerLineLayerData[]
    fixed_layer = _prepare_fixed_trunk_layer(fixed_trees)
    fixed_layer !== nothing && push!(line_layers, fixed_layer)
    grown_layer = _prepare_grown_tree_layer(forest, fixed_trees; max_segments_per_tree=grown_max_segments_per_tree)
    grown_layer !== nothing && push!(line_layers, grown_layer)

    return UnifiedViewerData(String(title), iodine.times, iodine.trees, point_layers, line_layers, iodine.cmax, iodine.blue_reference)
end

function export_unified_viewer_html(
    path::AbstractString,
    viewer::UnifiedViewerData;
    title::AbstractString=viewer.title,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "<!DOCTYPE html>")
        println(io, "<html lang=\"en\">")
        println(io, "<head>")
        println(io, "  <meta charset=\"utf-8\" />")
        println(io, "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />")
        println(io, "  <title>", title, "</title>")
        println(io, "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>")
        println(io, "  <style>")
        println(io, "    body { margin: 0; font-family: Arial, sans-serif; background: #f7f7f3; color: #222; }")
        println(io, "    #wrap { display: grid; grid-template-rows: auto auto auto 1fr; height: 100vh; }")
        println(io, "    #header { padding: 10px 16px 6px 16px; border-bottom: 1px solid #ddd; background: #fbfbf8; }")
        println(io, "    #controls { display: grid; grid-template-columns: auto 1fr auto; gap: 12px; align-items: center; padding: 10px 16px 4px 16px; border-bottom: 1px solid #e4e4df; background: #fdfdf9; }")
        println(io, "    #layerRow { padding: 0 16px 10px 16px; border-bottom: 1px solid #e4e4df; background: #fdfdf9; display: flex; flex-wrap: wrap; gap: 14px; }")
        println(io, "    #plot { width: 100%; height: calc(100vh - 168px); }")
        println(io, "    .meta { font-size: 13px; color: #555; }")
        println(io, "    .time-readout { font-size: 14px; font-variant-numeric: tabular-nums; min-width: 86px; text-align: right; }")
        println(io, "    .button-row { display: flex; gap: 8px; }")
        println(io, "    .layer-toggle { display: inline-flex; gap: 6px; align-items: center; font-size: 13px; color: #444; }")
        println(io, "    button { border: 1px solid #cfcfc8; background: white; color: #222; border-radius: 6px; padding: 6px 12px; cursor: pointer; }")
        println(io, "    button:hover { background: #f2f2eb; }")
        println(io, "    input[type=\"range\"] { width: 100%; }")
        println(io, "  </style>")
        println(io, "</head>")
        println(io, "<body>")
        println(io, "<div id=\"wrap\">")
        println(io, "  <div id=\"header\">")
        println(io, "    <div><strong>", title, "</strong></div>")
        println(io, "    <div class=\"meta\">Single interactive scene with toggleable geometry layers. Iodine colors: near 0 = gray, around 3 mg/mL = blue, near peak = red.</div>")
        println(io, "  </div>")
        println(io, "  <div id=\"controls\">")
        println(io, "    <div class=\"button-row\">")
        println(io, "      <button id=\"playButton\" type=\"button\">Play</button>")
        println(io, "      <button id=\"pauseButton\" type=\"button\">Pause</button>")
        println(io, "    </div>")
        println(io, "    <input id=\"timeSlider\" type=\"range\" min=\"0\" step=\"1\" value=\"0\" />")
        println(io, "    <div class=\"time-readout\" id=\"timeValue\">0.00 s</div>")
        println(io, "  </div>")
        println(io, "  <div id=\"layerRow\"></div>")
        println(io, "  <div id=\"plot\"></div>")
        println(io, "</div>")
        println(io, "<script>")
        print(io, "const viewerData = {\n  times: ")
        _write_js_vector(io, viewer.times)
        print(io, ",\n  cmax: ")
        _write_js_number(io, viewer.cmax <= 0 ? 1.0 : viewer.cmax)
        print(io, ",\n  blueReference: ")
        _write_js_number(io, viewer.blue_reference)
        print(io, ",\n  pointLayers: [\n")
        for (i, layer) in enumerate(viewer.point_layers)
            i > 1 && print(io, ",\n")
            print(io, "    { key: "); _write_js_string(io, layer.key)
            print(io, ", label: "); _write_js_string(io, layer.label)
            print(io, ", x: "); _write_js_vector(io, layer.x)
            print(io, ", y: "); _write_js_vector(io, layer.y)
            print(io, ", z: "); _write_js_vector(io, layer.z)
            print(io, ", size: "); _write_js_number(io, layer.size)
            print(io, ", color: "); _write_js_string(io, layer.color)
            print(io, ", opacity: "); _write_js_number(io, layer.opacity)
            print(io, ", visible: ", layer.visible ? "true" : "false")
            print(io, " }")
        end
        print(io, "\n  ],\n  lineLayers: [\n")
        for (i, layer) in enumerate(viewer.line_layers)
            i > 1 && print(io, ",\n")
            print(io, "    { key: "); _write_js_string(io, layer.key)
            print(io, ", label: "); _write_js_string(io, layer.label)
            print(io, ", x: "); _write_js_vector(io, layer.x)
            print(io, ", y: "); _write_js_vector(io, layer.y)
            print(io, ", z: "); _write_js_vector(io, layer.z)
            print(io, ", width: "); _write_js_number(io, layer.width)
            print(io, ", color: "); _write_js_string(io, layer.color)
            print(io, ", opacity: "); _write_js_number(io, layer.opacity)
            print(io, ", visible: ", layer.visible ? "true" : "false")
            print(io, " }")
        end
        print(io, "\n  ],\n  iodineTrees: [\n")
        for (i, tree) in enumerate(viewer.iodine_trees)
            i > 1 && print(io, ",\n")
            print(io, "    { name: "); _write_js_string(io, tree.name)
            print(io, ", pointX: "); _write_js_vector(io, tree.point_x)
            print(io, ", pointY: "); _write_js_vector(io, tree.point_y)
            print(io, ", pointZ: "); _write_js_vector(io, tree.point_z)
            print(io, ", pointSize: "); _write_js_vector(io, tree.point_size)
            print(io, ", lengthMm: "); _write_js_vector(io, tree.segment_length_mm)
            print(io, ", diameterUm: "); _write_js_vector(io, tree.segment_diameter_um)
            print(io, ", conc: "); _write_js_matrix_rows(io, tree.concentration_by_time)
            print(io, " }")
        end
        println(io, "\n  ]\n};")
        println(io, "const plotDiv = document.getElementById('plot');")
        println(io, "const slider = document.getElementById('timeSlider');")
        println(io, "const timeValue = document.getElementById('timeValue');")
        println(io, "const playButton = document.getElementById('playButton');")
        println(io, "const pauseButton = document.getElementById('pauseButton');")
        println(io, "const layerRow = document.getElementById('layerRow');")
        println(io, "const viewerColorscale = [[0.00, '#8c8c8c'],[0.15, '#b2b2b2'],[0.55, '#6a8dff'],[0.78, '#1048ff'],[1.00, '#d7191c']];")
        println(io, "function concentrationToColorValue(concentration) { const cmax = Math.max(viewerData.cmax, 1e-6); const blueReference = Math.max(Math.min(viewerData.blueReference, cmax), 1e-6); const lowFraction = 0.78; if (cmax <= blueReference + 1e-6) { return concentration; } if (concentration <= blueReference) { return (concentration / blueReference) * (lowFraction * cmax); } return lowFraction * cmax + ((concentration - blueReference) / (cmax - blueReference)) * ((1.0 - lowFraction) * cmax); }")
        println(io, "function displayColors(concentrations) { return concentrations.map(concentrationToColorValue); }")
        println(io, "function pointCustomData(tree, concentrations) { return concentrations.map((c, j) => [c, tree.lengthMm[j], tree.diameterUm[j]]); }")
        println(io, "const traces = []; const layerTraceMap = {};")
        println(io, "viewerData.lineLayers.forEach((layer) => { const traceIndex = traces.length; traces.push({ type: 'scatter3d', mode: 'lines', name: layer.label, x: layer.x, y: layer.y, z: layer.z, line: {color: layer.color, width: layer.width}, opacity: layer.opacity, hoverinfo: 'skip', visible: layer.visible, showlegend: false }); layerTraceMap[layer.key] = [traceIndex]; });")
        println(io, "viewerData.pointLayers.forEach((layer) => { const traceIndex = traces.length; traces.push({ type: 'scatter3d', mode: 'markers', name: layer.label, x: layer.x, y: layer.y, z: layer.z, marker: { size: layer.size, color: layer.color, opacity: layer.opacity }, hoverinfo: 'skip', visible: layer.visible, showlegend: false }); layerTraceMap[layer.key] = [traceIndex]; });")
        println(io, "const iodineTraceIndices = []; viewerData.iodineTrees.forEach((tree, idx) => { const traceIndex = traces.length; traces.push({ type: 'scatter3d', mode: 'markers', name: tree.name + ' iodine', x: tree.pointX, y: tree.pointY, z: tree.pointZ, marker: { size: tree.pointSize, color: displayColors(tree.conc[0]), cmin: 0, cmax: viewerData.cmax, colorscale: viewerColorscale, opacity: 0.98, showscale: idx === viewerData.iodineTrees.length - 1, colorbar: idx === viewerData.iodineTrees.length - 1 ? {title: 'iodine'} : undefined }, customdata: pointCustomData(tree, tree.conc[0]), hovertemplate: tree.name + '<br>x=%{x:.2f}<br>y=%{y:.2f}<br>z=%{z:.2f}<br>C=%{customdata[0]:.4f} mg/mL<br>L=%{customdata[1]:.3f} mm<br>D=%{customdata[2]:.1f} um<extra></extra>', visible: true }); iodineTraceIndices.push(traceIndex); }); layerTraceMap['iodine'] = iodineTraceIndices;")
        println(io, "const layout = { margin: {l: 0, r: 0, b: 0, t: 0}, paper_bgcolor: '#f7f7f3', plot_bgcolor: '#f7f7f3', scene: { aspectmode: 'data', bgcolor: '#f7f7f3', xaxis: {title: 'X', backgroundcolor: '#f7f7f3', gridcolor: '#dadad2'}, yaxis: {title: 'Y', backgroundcolor: '#f7f7f3', gridcolor: '#dadad2'}, zaxis: {title: 'Z', backgroundcolor: '#f7f7f3', gridcolor: '#dadad2'} }, legend: {orientation: 'h', yanchor: 'bottom', y: 0.02, xanchor: 'left', x: 0.02} };")
        println(io, "Plotly.newPlot(plotDiv, traces, layout, {responsive: true, displaylogo: false});")
        println(io, "slider.max = Math.max(viewerData.times.length - 1, 0);")
        println(io, "function updateTime(frameIdx) { const safeIdx = Math.max(0, Math.min(frameIdx, viewerData.times.length - 1)); slider.value = safeIdx; timeValue.textContent = viewerData.times[safeIdx].toFixed(2) + ' s'; viewerData.iodineTrees.forEach((tree, idx) => { const traceIndex = iodineTraceIndices[idx]; Plotly.restyle(plotDiv, {'marker.color': [displayColors(tree.conc[safeIdx])], 'customdata': [pointCustomData(tree, tree.conc[safeIdx])]}, [traceIndex]); }); }")
        println(io, "const layerConfig = []; viewerData.lineLayers.forEach(layer => layerConfig.push({key: layer.key, label: layer.label, checked: layer.visible})); viewerData.pointLayers.forEach(layer => layerConfig.push({key: layer.key, label: layer.label, checked: layer.visible})); layerConfig.push({key: 'iodine', label: 'Iodine Flow', checked: true});")
        println(io, "layerConfig.forEach(layer => { const label = document.createElement('label'); label.className = 'layer-toggle'; const input = document.createElement('input'); input.type = 'checkbox'; input.checked = layer.checked; input.addEventListener('change', () => { const indices = layerTraceMap[layer.key] || []; if (indices.length === 0) return; Plotly.restyle(plotDiv, {visible: input.checked}, indices); }); label.appendChild(input); label.appendChild(document.createTextNode(layer.label)); layerRow.appendChild(label); });")
        println(io, "slider.addEventListener('input', () => updateTime(Number(slider.value))); let timer = null; playButton.addEventListener('click', () => { if (timer !== null) return; timer = setInterval(() => { const next = (Number(slider.value) + 1) % viewerData.times.length; updateTime(next); }, 80); }); pauseButton.addEventListener('click', () => { if (timer !== null) { clearInterval(timer); timer = null; } }); updateTime(0);")
        println(io, "</script>")
        println(io, "</body>")
        println(io, "</html>")
    end
    return path
end
