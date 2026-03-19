using Test
using JLD2
using Unitful
using VesselTree
import XrayAttenuation as XA

@testset "Basis bridge" begin
    tmp = mktempdir()

    raw_data = reshape(UInt8.(1:8), 2, 2, 2)
    raw_path = joinpath(tmp, "toy.raw")
    open(raw_path, "w") do io
        write(io, raw_data)
    end

    loaded = load_xcat_raw_labels(raw_path; dims=(2, 2, 2), reverse_dims=(2,), voxel_size_mm=(0.5, 0.5, 1.0))
    @test size(loaded.mask) == (2, 2, 2)
    @test loaded.mask[1, 1, 1] == raw_data[1, 2, 1]
    @test loaded.voxel_size_mm == (0.5, 0.5, 1.0)

    tree = VascularTree("LAD", 4)
    seg_id = add_segment!(tree, (0.1, 0.5, 0.5), (1.9, 0.5, 0.5), 0.45, Int32(-1))
    tree.segments.flow[seg_id] = 1e-9
    cfg = TreeConfig("LAD", (0.1, 0.5, 0.5), 0.45, (1.0, 0.0, 0.0), 1, 1.0)
    domain = BoxDomain((0.0, 0.0, 0.0), (2.0, 2.0, 2.0))
    territory = initialize_territories(domain, [cfg])
    forest = CoronaryForest(Dict("LAD" => tree), territory, kassab_lad_params())

    result = ContrastTransportResult(Float64[0.0], reshape(Float32[2.5], 1, 1), nothing, Float64[0.1], Float64[0.01])
    base_materials = Dict(1 => XA.Materials.muscle)
    raw_mask = fill(UInt8(1), 3, 3, 3)
    alignment = XCATRawAlignment((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0))

    compiled = compile_basis_material_frame(
        raw_mask,
        base_materials,
        forest,
        Dict("LAD" => result),
        1;
        alignment=alignment,
        blood_fraction_bins=[0.0, 0.5, 1.0],
        iodine_concentration_bins_mg_mL=[0.0, 2.5, 5.0],
        supersample=1,
    )
    @test compiled.overlay_voxels > 0
    @test compiled.dynamic_label_count > 0
    @test maximum(compiled.mask) > 1

    saved = save_basis_frame(
        tmp,
        compiled.mask,
        compiled.materials_dict;
        time_s=0.0,
        voxel_size_cm=(0.1, 0.1, 0.1),
        origin_mm=(0.0, 0.0, 0.0),
        prefix="toy_basis",
    )
    @test isfile(saved.raw_path)
    @test isfile(saved.materials_path)
    @test isfile(saved.csv_path)
    @test isfile(saved.manifest_path)

    manifest = read(saved.manifest_path, String)
    @test occursin("voxel_size_cm", manifest)
    materials = JLD2.load(saved.materials_path, "materials_dict")
    @test length(materials) == length(compiled.materials_dict)
end
