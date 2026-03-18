using Sockets
using Printf

struct ContrastViewerTreeData
    name::String
    line_x::Vector{Float32}
    line_y::Vector{Float32}
    line_z::Vector{Float32}
    point_x::Vector{Float32}
    point_y::Vector{Float32}
    point_z::Vector{Float32}
    point_size::Vector{Float32}
    segment_length_mm::Vector{Float32}
    segment_diameter_um::Vector{Float32}
    concentration_by_time::Matrix{Float32}   # nt_display x n_points
end

struct ContrastViewerData
    times::Vector{Float32}
    trees::Vector{ContrastViewerTreeData}
    cmax::Float32
    blue_reference::Float32
end

function _sample_evenly(ids::Vector{Int}, target::Int)
    target <= 0 && return Int[]
    target >= length(ids) && return copy(ids)
    isempty(ids) && return Int[]

    idx = unique(round.(Int, range(1, length(ids), length=target)))
    chosen = ids[idx]
    if length(chosen) < target
        chosen_set = Set(chosen)
        for id in ids
            if !(id in chosen_set)
                push!(chosen, id)
                push!(chosen_set, id)
                length(chosen) >= target && break
            end
        end
    end
    return chosen[1:target]
end

function _viewer_segment_indices(
    tree::VascularTree;
    min_radius_um::Float64=0.0,
    max_segments::Int=5000,
    large_fraction::Float64=0.25,
    small_fraction::Float64=0.45,
)
    n = tree.segments.n
    radii_um = tree.segments.radius[1:n] .* 1000.0
    keep = findall(r -> r >= min_radius_um, radii_um)
    length(keep) <= max_segments && return keep

    selected = Int[]
    radii_keep = radii_um[keep]

    n_large = clamp(round(Int, large_fraction * max_segments), 1, max_segments)
    n_small = clamp(round(Int, small_fraction * max_segments), 1, max_segments - n_large)

    desc_order = sortperm(radii_keep; rev=true)
    asc_order = sortperm(radii_keep; rev=false)
    append!(selected, keep[desc_order[1:min(n_large, length(desc_order))]])
    append!(selected, keep[asc_order[1:min(n_small, length(asc_order))]])

    selected = unique(selected)
    remaining_needed = max_segments - length(selected)
    if remaining_needed > 0
        selected_set = Set(selected)
        remaining = [id for id in keep if !(id in selected_set)]
        append!(selected, _sample_evenly(remaining, remaining_needed))
    end

    if length(selected) > max_segments
        selected = _sample_evenly(sort(selected), max_segments)
    end

    sort!(selected)
    return selected
end

function _viewer_point_size(
    radius_mm::Float64;
    radius_scale::Float64=10.0,
    min_size::Float64=1.4,
    max_size::Float64=8.0,
)
    size = radius_scale * sqrt(max(radius_mm, 0.0))
    return Float32(clamp(size, min_size, max_size))
end

function prepare_contrast_viewer_data(
    forest::CoronaryForest,
    results::Dict{String, <:ContrastTransportResult};
    min_radius_um::Float64=4.0,
    max_segments_per_tree::Int=2500,
    time_stride::Int=1,
    radius_scale::Float64=10.0,
    blue_reference_mg_mL::Float64=3.0,
)
    tree_names = sort!(collect(keys(forest.trees)))
    trees = ContrastViewerTreeData[]
    cmax = 0.0f0

    for name in tree_names
        tree = forest.trees[name]
        result = results[name]
        seg = tree.segments
        keep = _viewer_segment_indices(tree; min_radius_um=min_radius_um, max_segments=max_segments_per_tree)
        time_idx = collect(1:time_stride:length(result.times))

        line_x = Float32[]
        line_y = Float32[]
        line_z = Float32[]
        point_x = Vector{Float32}(undef, length(keep))
        point_y = Vector{Float32}(undef, length(keep))
        point_z = Vector{Float32}(undef, length(keep))
        point_size = Vector{Float32}(undef, length(keep))
        segment_length_mm = Vector{Float32}(undef, length(keep))
        segment_diameter_um = Vector{Float32}(undef, length(keep))
        conc = Matrix{Float32}(undef, length(time_idx), length(keep))

        for (j, seg_id) in enumerate(keep)
            px = Float32(seg.proximal_x[seg_id])
            py = Float32(seg.proximal_y[seg_id])
            pz = Float32(seg.proximal_z[seg_id])
            dx = Float32(seg.distal_x[seg_id])
            dy = Float32(seg.distal_y[seg_id])
            dz = Float32(seg.distal_z[seg_id])

            push!(line_x, px, dx, NaN32)
            push!(line_y, py, dy, NaN32)
            push!(line_z, pz, dz, NaN32)

            point_x[j] = (px + dx) / 2.0f0
            point_y[j] = (py + dy) / 2.0f0
            point_z[j] = (pz + dz) / 2.0f0
            point_size[j] = _viewer_point_size(seg.radius[seg_id]; radius_scale=radius_scale)
            segment_length_mm[j] = Float32(seg.seg_length[seg_id])
            segment_diameter_um[j] = Float32(seg.radius[seg_id] * 2000.0)

            for (ti, src_ti) in enumerate(time_idx)
                value = Float32(result.concentration[seg_id, src_ti])
                conc[ti, j] = value
                cmax = max(cmax, value)
            end
        end

        push!(trees, ContrastViewerTreeData(name, line_x, line_y, line_z, point_x, point_y, point_z, point_size, segment_length_mm, segment_diameter_um, conc))
    end

    display_times = Float32.(results[first(tree_names)].times[1:time_stride:end])
    return ContrastViewerData(display_times, trees, cmax, Float32(blue_reference_mg_mL))
end

function _write_js_number(io, x::Real)
    if isfinite(x)
        @printf(io, "%.7g", Float64(x))
    else
        print(io, "null")
    end
end

function _write_js_string(io, s::AbstractString)
    esc = replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
    print(io, "\"", esc, "\"")
end

function _write_js_vector(io, data::AbstractVector)
    print(io, "[")
    for i in eachindex(data)
        i != firstindex(data) && print(io, ",")
        if data[i] isa AbstractString
            _write_js_string(io, data[i])
        else
            _write_js_number(io, data[i])
        end
    end
    print(io, "]")
end

function _write_js_matrix_rows(io, data::AbstractMatrix)
    print(io, "[")
    for i in axes(data, 1)
        i != first(axes(data, 1)) && print(io, ",")
        _write_js_vector(io, view(data, i, :))
    end
    print(io, "]")
end

function export_contrast_viewer_html(
    path::AbstractString,
    viewer::ContrastViewerData;
    title::AbstractString="Contrast Transport Viewer",
)
    mkpath(dirname(path))
    open(path, "w") do io
        print(io, """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>$(title)</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body { margin: 0; font-family: Arial, sans-serif; background: #f7f7f3; color: #222; }
    #wrap { display: grid; grid-template-rows: auto auto 1fr; height: 100vh; }
    #header { padding: 10px 16px 6px 16px; border-bottom: 1px solid #ddd; background: #fbfbf8; }
    #controls { display: grid; grid-template-columns: auto 1fr auto; gap: 12px; align-items: center; padding: 10px 16px; border-bottom: 1px solid #e4e4df; background: #fdfdf9; }
    #plot { width: 100%; height: calc(100vh - 120px); }
    .meta { font-size: 13px; color: #555; }
    .time-readout { font-size: 14px; font-variant-numeric: tabular-nums; min-width: 86px; text-align: right; }
    .button-row { display: flex; gap: 8px; }
    button { border: 1px solid #cfcfc8; background: white; color: #222; border-radius: 6px; padding: 6px 12px; cursor: pointer; }
    button:hover { background: #f2f2eb; }
    input[type="range"] { width: 100%; }
  </style>
</head>
<body>
<div id="wrap">
  <div id="header">
    <div><strong>$(title)</strong></div>
    <div class="meta">Dynamic segment-average iodine concentration (mg/mL). Gray lines: vascular skeleton. Colors emphasize the low-concentration range: near 0 = gray, around 3 mg/mL = blue, near peak = red.</div>
  </div>
  <div id="controls">
    <div class="button-row">
      <button id="playButton" type="button">Play</button>
      <button id="pauseButton" type="button">Pause</button>
    </div>
    <input id="timeSlider" type="range" min="0" step="1" value="0" />
    <div class="time-readout" id="timeValue">0.00 s</div>
  </div>
  <div id="plot"></div>
</div>
<script>
const viewerData = {
  times:
""")
        _write_js_vector(io, viewer.times)
        print(io, ",\n  cmax: ")
        _write_js_number(io, viewer.cmax <= 0 ? 1.0 : viewer.cmax)
        print(io, ",\n  blueReference: ")
        _write_js_number(io, viewer.blue_reference)
        print(io, ",\n  trees: [\n")
        for (i, tree) in enumerate(viewer.trees)
            i > 1 && print(io, ",\n")
            print(io, "    { name: ")
            _write_js_string(io, tree.name)
            print(io, ", lineX: "); _write_js_vector(io, tree.line_x)
            print(io, ", lineY: "); _write_js_vector(io, tree.line_y)
            print(io, ", lineZ: "); _write_js_vector(io, tree.line_z)
            print(io, ", pointX: "); _write_js_vector(io, tree.point_x)
            print(io, ", pointY: "); _write_js_vector(io, tree.point_y)
            print(io, ", pointZ: "); _write_js_vector(io, tree.point_z)
            print(io, ", pointSize: "); _write_js_vector(io, tree.point_size)
            print(io, ", lengthMm: "); _write_js_vector(io, tree.segment_length_mm)
            print(io, ", diameterUm: "); _write_js_vector(io, tree.segment_diameter_um)
            print(io, ", conc: "); _write_js_matrix_rows(io, tree.concentration_by_time)
            print(io, " }")
        end
        print(io, """
  ]
};

const plotDiv = document.getElementById('plot');
const slider = document.getElementById('timeSlider');
const timeValue = document.getElementById('timeValue');
const playButton = document.getElementById('playButton');
const pauseButton = document.getElementById('pauseButton');
const viewerColorscale = [
  [0.00, '#8c8c8c'],
  [0.15, '#b2b2b2'],
  [0.55, '#6a8dff'],
  [0.78, '#1048ff'],
  [1.00, '#d7191c']
];

function concentrationToColorValue(concentration) {
  const cmax = Math.max(viewerData.cmax, 1e-6);
  const blueReference = Math.max(Math.min(viewerData.blueReference, cmax), 1e-6);
  const lowFraction = 0.78;
  if (cmax <= blueReference + 1e-6) {
    return concentration;
  }
  if (concentration <= blueReference) {
    return (concentration / blueReference) * (lowFraction * cmax);
  }
  return lowFraction * cmax + ((concentration - blueReference) / (cmax - blueReference)) * ((1.0 - lowFraction) * cmax);
}

function displayColors(concentrations) {
  return concentrations.map(concentrationToColorValue);
}

function pointCustomData(tree, concentrations) {
  return concentrations.map((c, j) => [c, tree.lengthMm[j], tree.diameterUm[j]]);
}

const traces = [];
viewerData.trees.forEach((tree, idx) => {
  traces.push({
    type: 'scatter3d',
    mode: 'lines',
    name: tree.name + ' skeleton',
    x: tree.lineX,
    y: tree.lineY,
    z: tree.lineZ,
    line: {color: 'rgba(120,120,120,0.10)', width: 1.0},
    hoverinfo: 'skip',
    showlegend: false
  });
  traces.push({
    type: 'scatter3d',
    mode: 'markers',
    name: tree.name,
    x: tree.pointX,
    y: tree.pointY,
    z: tree.pointZ,
    marker: {
      size: tree.pointSize,
      color: displayColors(tree.conc[0]),
      cmin: 0,
      cmax: viewerData.cmax,
      colorscale: viewerColorscale,
      opacity: 0.98,
      showscale: idx === viewerData.trees.length - 1,
      colorbar: idx === viewerData.trees.length - 1 ? {title: 'iodine'} : undefined
    },
    customdata: pointCustomData(tree, tree.conc[0]),
    text: tree.pointSize.map((_, j) => tree.name + ' seg ' + (j + 1)),
    hovertemplate: tree.name + '<br>x=%{x:.2f}<br>y=%{y:.2f}<br>z=%{z:.2f}<br>C=%{customdata[0]:.4f} mg/mL<br>L=%{customdata[1]:.3f} mm<br>D=%{customdata[2]:.1f} um<extra></extra>'
  });
});

const layout = {
  margin: {l: 0, r: 0, b: 0, t: 0},
  paper_bgcolor: '#f7f7f3',
  plot_bgcolor: '#f7f7f3',
  scene: {
    aspectmode: 'data',
    xaxis: {title: 'X'},
    yaxis: {title: 'Y'},
    zaxis: {title: 'Z'}
  }
};

slider.max = String(Math.max(viewerData.times.length - 1, 0));

let currentIndex = 0;
let timerId = null;

function updateFrame(index) {
  currentIndex = Math.max(0, Math.min(index, viewerData.times.length - 1));
  slider.value = String(currentIndex);
  timeValue.textContent = viewerData.times[currentIndex].toFixed(2) + ' s';

  viewerData.trees.forEach((tree, idx) => {
    const actual = tree.conc[currentIndex];
    Plotly.restyle(plotDiv, {
      'marker.color': [displayColors(actual)],
      'marker.cmin': [0],
      'marker.cmax': [viewerData.cmax],
      'marker.colorscale': [viewerColorscale],
      'marker.opacity': [0.98],
      'customdata': [pointCustomData(tree, actual)]
    }, [2 * idx + 1]);
  });
}

function startPlayback() {
  if (timerId !== null) {
    return;
  }
  timerId = window.setInterval(() => {
    const nextIndex = currentIndex + 1;
    if (nextIndex >= viewerData.times.length) {
      stopPlayback();
      return;
    }
    updateFrame(nextIndex);
  }, 80);
}

function stopPlayback() {
  if (timerId !== null) {
    window.clearInterval(timerId);
    timerId = null;
  }
}

slider.addEventListener('input', (event) => {
  stopPlayback();
  updateFrame(Number(event.target.value));
});
playButton.addEventListener('click', startPlayback);
pauseButton.addEventListener('click', stopPlayback);

Plotly.newPlot(plotDiv, traces, layout, {responsive: true}).then(() => {
  updateFrame(0);
});
</script>
</body>
</html>
""")
    end
    return path
end

function _http_content_type(path::AbstractString)
    lower = lowercase(path)
    endswith(lower, ".html") && return "text/html; charset=utf-8"
    endswith(lower, ".js") && return "application/javascript; charset=utf-8"
    endswith(lower, ".css") && return "text/css; charset=utf-8"
    endswith(lower, ".json") && return "application/json; charset=utf-8"
    endswith(lower, ".png") && return "image/png"
    endswith(lower, ".jpg") && return "image/jpeg"
    return "application/octet-stream"
end

function _http_response(io::IO, status::Int, body::AbstractVector{UInt8}; content_type::String="text/plain; charset=utf-8")
    status_text = status == 200 ? "OK" :
                  status == 404 ? "Not Found" :
                  status == 403 ? "Forbidden" :
                  status == 405 ? "Method Not Allowed" : "Error"
    write(io, "HTTP/1.1 $status $status_text\r\n")
    write(io, "Content-Type: $content_type\r\n")
    write(io, "Content-Length: $(length(body))\r\n")
    write(io, "Cache-Control: no-store\r\n")
    write(io, "Connection: close\r\n\r\n")
    write(io, body)
end

function serve_static_directory(
    directory::AbstractString;
    host::AbstractString="127.0.0.1",
    port::Int=8008,
)
    root = normpath(abspath(directory))
    server = Sockets.listen(parse(IPAddr, host), port)

    @async begin
        while true
            sock = accept(server)
            @async begin
                try
                    request_line = String(readuntil(sock, "\r\n"))
                    parts = split(chomp(request_line))
                    if length(parts) < 2 || !(parts[1] in ("GET", "HEAD"))
                        while isopen(sock)
                            line = String(readuntil(sock, "\r\n"))
                            isempty(strip(line)) && break
                        end
                        _http_response(sock, 405, Vector{UInt8}(codeunits("Only GET and HEAD are supported.")))
                    else
                        while isopen(sock)
                            line = String(readuntil(sock, "\r\n"))
                            isempty(strip(line)) && break
                        end

                        raw_target = split(parts[2], '?')[1]
                        rel = raw_target == "/" ? "index.html" : lstrip(raw_target, '/')
                        full = normpath(joinpath(root, rel))
                        if !startswith(full, root)
                            _http_response(sock, 403, Vector{UInt8}(codeunits("Forbidden.")))
                        elseif isfile(full)
                            body = read(full)
                            _http_response(sock, 200, body; content_type=_http_content_type(full))
                        else
                            _http_response(sock, 404, Vector{UInt8}(codeunits("Not found.")))
                        end
                    end
                catch err
                    msg = sprint(showerror, err)
                    _http_response(sock, 500, Vector{UInt8}(codeunits(msg)))
                finally
                    close(sock)
                end
            end
        end
    end

    return server
end
