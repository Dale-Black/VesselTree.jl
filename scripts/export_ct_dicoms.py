#!/usr/bin/env python3

import argparse
import json
import shutil
from datetime import datetime, timedelta
from pathlib import Path

import h5py
import numpy as np
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export HU-calibrated GE scan reconstructions to DICOM series."
    )
    parser.add_argument("--run-manifest", required=True, help="Path to ge_scan_run.json")
    parser.add_argument("--output-dir", required=True, help="Output directory for DICOM series")
    parser.add_argument("--overwrite", action="store_true", help="Delete output dir if it exists")
    parser.add_argument("--patient-name", default="XCAT^VMale50")
    parser.add_argument("--patient-id", default="XCAT-VM50")
    parser.add_argument("--study-description", default="Simulated Dynamic Coronary CT")
    parser.add_argument("--series-description-prefix", default="Simulated Coronary CT")
    parser.add_argument("--manufacturer", default="BasisSimulator")
    parser.add_argument("--manufacturer-model-name", default="GE Revolution Apex Elite (simulated)")
    parser.add_argument("--institution-name", default="MolloiLab")
    parser.add_argument("--body-part", default="HEART")
    parser.add_argument("--kvp", type=float, default=120.0)
    parser.add_argument("--window-center", type=float, default=300.0)
    parser.add_argument("--window-width", type=float, default=1200.0)
    return parser.parse_args()


def read_struct_scalar(node):
    value = node[()]
    if isinstance(value, np.void) and value.dtype.names:
        return tuple(value[name].item() for name in value.dtype.names)
    if isinstance(value, np.ndarray) and value.dtype.names and value.shape == ():
        return tuple(value[name].item() for name in value.dtype.names)
    if hasattr(value, "item"):
        return value.item()
    return value


def load_scan_payload(scan_path):
    with h5py.File(scan_path, "r") as f:
        volume = np.asarray(f["recon_hu"][()], dtype=np.float32)
        voxel_size_cm = tuple(float(x) for x in read_struct_scalar(f["voxel_size_cm"]))
        origin_cm = tuple(float(x) for x in read_struct_scalar(f["origin_cm"]))
        matrix_size = tuple(int(x) for x in read_struct_scalar(f["recon_matrix_size"]))
    return {
        "volume_hu": volume,
        "voxel_size_cm": voxel_size_cm,
        "origin_cm": origin_cm,
        "matrix_size": matrix_size,
    }


def build_file_meta():
    meta = FileMetaDataset()
    meta.FileMetaInformationVersion = b"\x00\x01"
    meta.MediaStorageSOPClassUID = CTImageStorage
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    meta.ImplementationClassUID = generate_uid()
    meta.ImplementationVersionName = "V4DICOM1"
    return meta


def make_dataset(*, file_meta, study_uid, series_uid, frame_uid, sop_uid, patient_name, patient_id,
                 study_description, series_description, manufacturer, manufacturer_model_name,
                 institution_name, body_part, kvp, window_center, window_width, time_s,
                 slice_index, nslices, temporal_count, pixel_spacing_mm, slice_thickness_mm,
                 image_position_mm, image_orientation, rows, cols, pixel_array,
                 study_date, study_time):
    ds = FileDataset(None, {}, file_meta=file_meta, preamble=b"\0" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR = False

    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = sop_uid
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.FrameOfReferenceUID = frame_uid

    ds.PatientName = patient_name
    ds.PatientID = patient_id
    ds.PatientSex = "M"

    ds.StudyDate = study_date
    ds.SeriesDate = study_date
    ds.AcquisitionDate = study_date
    ds.ContentDate = study_date

    acquisition_dt = datetime.strptime(study_date + study_time, "%Y%m%d%H%M%S") + timedelta(seconds=float(time_s))
    time_token = acquisition_dt.strftime("%H%M%S")
    ds.StudyTime = study_time
    ds.SeriesTime = time_token
    ds.AcquisitionTime = time_token
    ds.ContentTime = time_token

    ds.Modality = "CT"
    ds.ImageType = ["DERIVED", "SECONDARY", "AXIAL"]
    ds.StudyDescription = study_description
    ds.SeriesDescription = series_description
    ds.ProtocolName = study_description
    ds.Manufacturer = manufacturer
    ds.ManufacturerModelName = manufacturer_model_name
    ds.InstitutionName = institution_name
    ds.BodyPartExamined = body_part

    ds.SeriesNumber = int(round(time_s)) + 1
    ds.AcquisitionNumber = int(round(time_s)) + 1
    ds.InstanceNumber = slice_index + 1
    ds.ImagesInAcquisition = nslices
    ds.TemporalPositionIdentifier = int(round(time_s)) + 1
    ds.NumberOfTemporalPositions = temporal_count
    ds.TriggerTime = int(round(float(time_s) * 1000.0))

    ds.KVP = float(kvp)
    ds.SliceThickness = float(slice_thickness_mm)
    ds.SpacingBetweenSlices = float(slice_thickness_mm)
    ds.PixelSpacing = [float(pixel_spacing_mm[0]), float(pixel_spacing_mm[1])]
    ds.ImagePositionPatient = [float(v) for v in image_position_mm]
    ds.ImageOrientationPatient = [float(v) for v in image_orientation]
    ds.SliceLocation = float(image_position_mm[2])
    ds.PatientPosition = "HFS"

    ds.RescaleIntercept = 0
    ds.RescaleSlope = 1
    ds.RescaleType = "HU"
    ds.WindowCenter = float(window_center)
    ds.WindowWidth = float(window_width)
    ds.ConvolutionKernel = "STANDARD"
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.SamplesPerPixel = 1
    ds.Rows = int(rows)
    ds.Columns = int(cols)
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1
    ds.SmallestImagePixelValue = int(np.min(pixel_array))
    ds.LargestImagePixelValue = int(np.max(pixel_array))
    ds.PixelData = pixel_array.astype(np.int16, copy=False).tobytes()
    return ds


def export_series(scan_entry, output_dir, args, study_uid, frame_uid, study_date, study_time, temporal_count):
    payload = load_scan_payload(scan_entry["scan_jld2_path"])
    volume = payload["volume_hu"]
    voxel_size_mm = tuple(v * 10.0 for v in payload["voxel_size_cm"])
    origin_mm = tuple(v * 10.0 for v in payload["origin_cm"])

    nz, ny, nx = volume.shape
    pixel_spacing_mm = (voxel_size_mm[1], voxel_size_mm[0])
    slice_thickness_mm = voxel_size_mm[2]
    image_orientation = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

    time_s = float(scan_entry["time_s"])
    time_label = f"t{int(round(time_s)):03d}"
    series_dir = output_dir / time_label
    series_dir.mkdir(parents=True, exist_ok=True)

    series_uid = generate_uid()
    series_description = f"{args.series_description_prefix} t={time_s:.1f}s"
    file_meta = build_file_meta()

    written = []
    for iz in range(nz):
        sop_uid = generate_uid()
        file_meta.MediaStorageSOPInstanceUID = sop_uid
        slice_hu = np.clip(np.rint(volume[iz]), -32768, 32767).astype(np.int16)
        image_position_mm = (
            origin_mm[0],
            origin_mm[1],
            origin_mm[2] + iz * slice_thickness_mm,
        )
        ds = make_dataset(
            file_meta=file_meta,
            study_uid=study_uid,
            series_uid=series_uid,
            frame_uid=frame_uid,
            sop_uid=sop_uid,
            patient_name=args.patient_name,
            patient_id=args.patient_id,
            study_description=args.study_description,
            series_description=series_description,
            manufacturer=args.manufacturer,
            manufacturer_model_name=args.manufacturer_model_name,
            institution_name=args.institution_name,
            body_part=args.body_part,
            kvp=args.kvp,
            window_center=args.window_center,
            window_width=args.window_width,
            time_s=time_s,
            slice_index=iz,
            nslices=nz,
            temporal_count=temporal_count,
            pixel_spacing_mm=pixel_spacing_mm,
            slice_thickness_mm=slice_thickness_mm,
            image_position_mm=image_position_mm,
            image_orientation=image_orientation,
            rows=ny,
            cols=nx,
            pixel_array=slice_hu,
            study_date=study_date,
            study_time=study_time,
        )
        out_path = series_dir / f"IM{iz + 1:04d}.dcm"
        ds.save_as(str(out_path), write_like_original=False)
        written.append(str(out_path))

    return {
        "time_s": time_s,
        "series_dir": str(series_dir),
        "series_instance_uid": series_uid,
        "scan_jld2_path": scan_entry["scan_jld2_path"],
        "scan_manifest_path": scan_entry["scan_manifest_path"],
        "slice_count": nz,
        "rows": ny,
        "cols": nx,
        "voxel_size_mm": list(voxel_size_mm),
        "files": written,
    }


def main():
    args = parse_args()
    run_manifest_path = Path(args.run_manifest).resolve()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists() and args.overwrite:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_manifest = json.loads(run_manifest_path.read_text())
    scan_entries = run_manifest["scan_manifests"]
    study_uid = generate_uid()
    frame_uid = generate_uid()
    now = datetime.now()
    study_date = now.strftime("%Y%m%d")
    study_time = now.strftime("%H%M%S")
    temporal_count = len(scan_entries)

    exported = []
    for scan_entry in scan_entries:
        exported.append(export_series(scan_entry, output_dir, args, study_uid, frame_uid, study_date, study_time, temporal_count))

    manifest = {
        "source_run_manifest": str(run_manifest_path),
        "study_instance_uid": study_uid,
        "frame_of_reference_uid": frame_uid,
        "study_date": study_date,
        "study_time": study_time,
        "series": exported,
    }
    manifest_path = output_dir / 'dicom_export_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"DICOM export written to: {output_dir}")
    print(f"Manifest: {manifest_path}")


if __name__ == '__main__':
    main()
