import JSON

struct CTViewerData
    title::String
    times::Vector{Float32}
    volumes::Vector{Array{Float32, 3}}
    voxel_size_cm::NTuple{3, Float32}
    value_unit::String
    display_min::Float32
    display_max::Float32
end

function _float_tuple3(x)
    vals = Float32.(collect(x))
    length(vals) == 3 || error("expected length-3 tuple/vector")
    return (vals[1], vals[2], vals[3])
end

function _quantile_value(values::Vector{Float32}, q::Float64)
    isempty(values) && return 0.0f0
    qq = clamp(q, 0.0, 1.0)
    sorted = sort(values)
    idx = clamp(round(Int, 1 + qq * (length(sorted) - 1)), 1, length(sorted))
    return sorted[idx]
end

function _preferred_recon_field(data::AbstractDict; prefer_hu::Bool=true)
    if prefer_hu && haskey(data, "recon_hu") && !(data["recon_hu"] isa Nothing)
        return (Float32.(data["recon_hu"]), "HU")
    elseif haskey(data, "recon_mu")
        return (Float32.(data["recon_mu"]), raw"μ")
    elseif haskey(data, "recon_hu") && !(data["recon_hu"] isa Nothing)
        return (Float32.(data["recon_hu"]), "HU")
    else
        error("scan file does not contain recon_mu or recon_hu")
    end
end

function load_ct_viewer_data(
    run_manifest_path::AbstractString;
    prefer_hu::Bool=true,
    q_low::Float64=0.01,
    q_high::Float64=0.995,
    title::AbstractString="XCAT CT Time-Series Viewer",
)
    run_manifest = JSON.parsefile(run_manifest_path)
    scan_entries = run_manifest["scan_manifests"]
    isempty(scan_entries) && error("scan manifest contains no frames")

    times = Float32[]
    volumes = Array{Float32, 3}[]
    voxel_size_cm = nothing
    value_unit = nothing
    pooled = Float32[]

    sorted_entries = sort(scan_entries; by=entry -> Float64(entry["time_s"]))
    for entry in sorted_entries
        scan_path = String(entry["scan_jld2_path"])
        loaded = JLD2.load(scan_path)
        volume, unit = _preferred_recon_field(loaded; prefer_hu=prefer_hu)
        push!(times, Float32(entry["time_s"]))
        push!(volumes, volume)
        voxel_size_cm = isnothing(voxel_size_cm) ? _float_tuple3(loaded["voxel_size_cm"]) : voxel_size_cm
        value_unit = isnothing(value_unit) ? unit : value_unit
        append!(pooled, vec(volume))
    end

    display_min = _quantile_value(pooled, q_low)
    display_max = _quantile_value(pooled, q_high)
    if !(display_max > display_min)
        display_min = minimum(pooled)
        display_max = maximum(pooled)
        if !(display_max > display_min)
            display_max = display_min + 1.0f0
        end
    end

    return CTViewerData(String(title), times, volumes, voxel_size_cm::NTuple{3, Float32}, String(value_unit), display_min, display_max)
end

function _write_js_cube(io, data::Array{Float32, 3})
    print(io, "[")
    for k in axes(data, 3)
        k != first(axes(data, 3)) && print(io, ",")
        print(io, "[")
        for j in axes(data, 2)
            j != first(axes(data, 2)) && print(io, ",")
            _write_js_vector(io, view(data, :, j, k))
        end
        print(io, "]")
    end
    print(io, "]")
end

function export_ct_viewer_html(
    path::AbstractString,
    viewer::CTViewerData;
    title::AbstractString=viewer.title,
)
    mkpath(dirname(path))
    dims = size(first(viewer.volumes))
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
        println(io, "    #wrap { display: grid; grid-template-rows: auto auto auto auto 1fr; height: 100vh; }")
        println(io, "    #header { padding: 12px 16px 6px 16px; border-bottom: 1px solid #ddd; background: #fbfbf8; }")
        println(io, "    .meta { font-size: 13px; color: #555; }")
        println(io, "    #timeControls, #sliceControls, #windowControls { display: grid; grid-template-columns: auto 1fr auto; gap: 12px; align-items: center; padding: 10px 16px; border-bottom: 1px solid #e4e4df; background: #fdfdf9; }")
        println(io, "    .button-row { display: flex; gap: 8px; align-items: center; }")
        println(io, "    .readout { font-size: 14px; font-variant-numeric: tabular-nums; min-width: 118px; text-align: right; }")
        println(io, "    #plot { width: 100%; height: calc(100vh - 220px); }")
        println(io, "    button, select { border: 1px solid #cfcfc8; background: white; color: #222; border-radius: 6px; padding: 6px 12px; cursor: pointer; }")
        println(io, "    button:hover, select:hover { background: #f2f2eb; }")
        println(io, "    input[type=\"range\"] { width: 100%; }")
        println(io, "    .window-note { font-size: 12px; color: #666; }")
        println(io, "  </style>")
        println(io, "</head>")
        println(io, "<body>")
        println(io, "<div id=\"wrap\">")
        println(io, "  <div id=\"header\">")
        println(io, "    <div><strong>", title, "</strong></div>")
        println(io, "    <div class=\"meta\">Interactive CT time-series viewer for reconstructed GE scan outputs. Use the time slider for t=0..9 s, switch slice plane, and browse the reconstructed volume.</div>")
        println(io, "  </div>")
        println(io, "  <div id=\"timeControls\">")
        println(io, "    <div class=\"button-row\"><button id=\"playButton\" type=\"button\">Play</button><button id=\"pauseButton\" type=\"button\">Pause</button></div>")
        println(io, "    <input id=\"timeSlider\" type=\"range\" min=\"0\" step=\"1\" value=\"0\" />")
        println(io, "    <div class=\"readout\" id=\"timeValue\">0.00 s</div>")
        println(io, "  </div>")
        println(io, "  <div id=\"sliceControls\">")
        println(io, "    <div class=\"button-row\"><label for=\"planeSelect\">Plane</label><select id=\"planeSelect\"><option value=\"axial\">Axial</option><option value=\"coronal\">Coronal</option><option value=\"sagittal\">Sagittal</option></select></div>")
        println(io, "    <input id=\"sliceSlider\" type=\"range\" min=\"0\" step=\"1\" value=\"0\" />")
        println(io, "    <div class=\"readout\" id=\"sliceValue\">Slice 0</div>")
        println(io, "  </div>")
        println(io, "  <div id=\"windowControls\">")
        println(io, "    <div class=\"button-row\"><span>Window</span><span class=\"window-note\">unit: $(viewer.value_unit)</span></div>")
        println(io, "    <input id=\"windowSlider\" type=\"range\" min=\"0\" step=\"1\" value=\"0\" />")
        println(io, "    <div class=\"readout\" id=\"windowValue\"></div>")
        println(io, "  </div>")
        println(io, "  <div id=\"plot\"></div>")
        println(io, "</div>")
        println(io, "<script>")
        print(io, "const viewerData = {\n  title: ")
        _write_js_string(io, viewer.title)
        print(io, ",\n  valueUnit: ")
        _write_js_string(io, viewer.value_unit)
        print(io, ",\n  displayMin: ")
        _write_js_number(io, viewer.display_min)
        print(io, ",\n  displayMax: ")
        _write_js_number(io, viewer.display_max)
        print(io, ",\n  voxelSizeCm: ")
        _write_js_vector(io, collect(viewer.voxel_size_cm))
        print(io, ",\n  dims: ")
        _write_js_vector(io, collect(Float32.(dims)))
        print(io, ",\n  times: ")
        _write_js_vector(io, viewer.times)
        print(io, ",\n  volumes: [\n")
        for (i, volume) in enumerate(viewer.volumes)
            i > 1 && print(io, ",\n")
            print(io, "    ")
            _write_js_cube(io, volume)
        end
        println(io, "\n  ]\n};")
        println(io, "const timeSlider = document.getElementById('timeSlider');")
        println(io, "const sliceSlider = document.getElementById('sliceSlider');")
        println(io, "const windowSlider = document.getElementById('windowSlider');")
        println(io, "const timeValue = document.getElementById('timeValue');")
        println(io, "const sliceValue = document.getElementById('sliceValue');")
        println(io, "const windowValue = document.getElementById('windowValue');")
        println(io, "const planeSelect = document.getElementById('planeSelect');")
        println(io, "const playButton = document.getElementById('playButton');")
        println(io, "const pauseButton = document.getElementById('pauseButton');")
        println(io, "const plotDiv = document.getElementById('plot');")
        println(io, "const windowPresets = [")
        println(io, "  {label: 'Auto', lo: viewerData.displayMin, hi: viewerData.displayMax},")
        println(io, "  {label: 'Tighter', lo: viewerData.displayMin + 0.15 * (viewerData.displayMax - viewerData.displayMin), hi: viewerData.displayMax - 0.10 * (viewerData.displayMax - viewerData.displayMin)},")
        println(io, "  {label: 'Wider', lo: viewerData.displayMin - 0.10 * (viewerData.displayMax - viewerData.displayMin), hi: viewerData.displayMax + 0.10 * (viewerData.displayMax - viewerData.displayMin)}")
        println(io, "]; ")
        println(io, "timeSlider.max = Math.max(viewerData.times.length - 1, 0);")
        println(io, "windowSlider.max = Math.max(windowPresets.length - 1, 0);")
        println(io, "function clamp(v, lo, hi) { return Math.max(lo, Math.min(hi, v)); }")
        println(io, "function currentWindow() { return windowPresets[Number(windowSlider.value)]; }")
        println(io, "function sliceCount(plane) { const dims = viewerData.dims.map(v => Number(v)); if (plane === 'axial') return dims[2]; if (plane === 'coronal') return dims[1]; return dims[0]; }")
        println(io, "function updateSliceBounds() { const plane = planeSelect.value; const count = sliceCount(plane); sliceSlider.max = Math.max(count - 1, 0); if (Number(sliceSlider.value) > count - 1) sliceSlider.value = Math.floor((count - 1) / 2); }")
        println(io, "function getSlice(frameIdx, plane, sliceIdx) { const volume = viewerData.volumes[frameIdx]; const nx = volume[0][0].length; const ny = volume[0].length; const nz = volume.length; if (plane === 'axial') { const k = clamp(sliceIdx, 0, nz - 1); const rows = []; for (let j = ny - 1; j >= 0; --j) rows.push(volume[k][j].slice()); return rows; } if (plane === 'coronal') { const j = clamp(sliceIdx, 0, ny - 1); const rows = []; for (let k = nz - 1; k >= 0; --k) { const row = []; for (let i = 0; i < nx; ++i) row.push(volume[k][j][i]); rows.push(row); } return rows; } const i = clamp(sliceIdx, 0, nx - 1); const rows = []; for (let k = nz - 1; k >= 0; --k) { const row = []; for (let j = 0; j < ny; ++j) row.push(volume[k][j][i]); rows.push(row.reverse()); } return rows; }")
        println(io, "function updateReadouts() { timeValue.textContent = viewerData.times[Number(timeSlider.value)].toFixed(2) + ' s'; const plane = planeSelect.value; const idx = Number(sliceSlider.value); sliceValue.textContent = plane.charAt(0).toUpperCase() + plane.slice(1) + ' slice ' + idx; const win = currentWindow(); windowValue.textContent = win.label + ' [' + win.lo.toFixed(4) + ', ' + win.hi.toFixed(4) + '] ' + viewerData.valueUnit; }")
        println(io, "function drawCurrentSlice() { const frameIdx = Number(timeSlider.value); const plane = planeSelect.value; const sliceIdx = Number(sliceSlider.value); const win = currentWindow(); const slice = getSlice(frameIdx, plane, sliceIdx); const trace = { type: 'heatmap', z: slice, colorscale: 'Gray', reversescale: false, zmin: win.lo, zmax: win.hi, showscale: true, colorbar: {title: viewerData.valueUnit} }; const layout = { margin: {l: 10, r: 20, t: 36, b: 10}, paper_bgcolor: '#f7f7f3', plot_bgcolor: '#f7f7f3', xaxis: {visible: false, constrain: 'domain'}, yaxis: {visible: false, scaleanchor: 'x'}, title: plane.charAt(0).toUpperCase() + plane.slice(1) + ' view' }; Plotly.react(plotDiv, [trace], layout, {responsive: true, displaylogo: false}); updateReadouts(); }")
        println(io, "timeSlider.addEventListener('input', drawCurrentSlice);")
        println(io, "sliceSlider.addEventListener('input', drawCurrentSlice);")
        println(io, "windowSlider.addEventListener('input', drawCurrentSlice);")
        println(io, "planeSelect.addEventListener('change', () => { updateSliceBounds(); drawCurrentSlice(); });")
        println(io, "let timer = null; playButton.addEventListener('click', () => { if (timer !== null) return; timer = setInterval(() => { const next = (Number(timeSlider.value) + 1) % viewerData.times.length; timeSlider.value = next; drawCurrentSlice(); }, 500); });")
        println(io, "pauseButton.addEventListener('click', () => { if (timer !== null) { clearInterval(timer); timer = null; } });")
        println(io, "updateSliceBounds(); sliceSlider.value = Math.floor(sliceCount(planeSelect.value) / 2); drawCurrentSlice();")
        println(io, "</script>")
        println(io, "</body>")
        println(io, "</html>")
    end
    return path
end
