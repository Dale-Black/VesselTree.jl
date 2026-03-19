using Test
using JLD2
using JSON
using VesselTree

@testset "CT viewer" begin
    tmp = mktempdir()
    scan_path = joinpath(tmp, "frame_scan.jld2")
    recon_mu = reshape(Float32.(1:24), 4, 3, 2)
    JLD2.jldsave(scan_path;
        recon_mu=recon_mu,
        recon_hu=nothing,
        voxel_size_cm=(0.1, 0.1, 0.2),
        origin_cm=(0.0, 0.0, 0.0),
        crop_bbox_ijk=(1:4, 1:3, 1:2),
        recon_matrix_size=size(recon_mu),
        frame_manifest_path="dummy",
        mu_water=nothing,
    )
    run_manifest_path = joinpath(tmp, "ge_scan_run.json")
    open(run_manifest_path, "w") do io
        JSON.print(io, Dict(
            "scan_manifests" => [Dict(
                "time_s" => 0.0,
                "scan_jld2_path" => scan_path,
                "scan_manifest_path" => joinpath(tmp, "frame_scan.json"),
            )],
        ))
    end

    viewer = load_ct_viewer_data(run_manifest_path; prefer_hu=true)
    @test length(viewer.times) == 1
    @test viewer.value_unit == raw"μ"
    @test size(first(viewer.volumes)) == (4, 3, 2)
    @test viewer.display_max > viewer.display_min

    html_path = joinpath(tmp, "index.html")
    export_ct_viewer_html(html_path, viewer)
    html = read(html_path, String)
    @test occursin("timeSlider", html)
    @test occursin("planeSelect", html)
    @test occursin("sliceSlider", html)
end
