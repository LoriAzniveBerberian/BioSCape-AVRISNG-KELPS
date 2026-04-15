"""
extract_kelp_spectra.py

Extracts per-pixel kelp spectra from AVIRIS-NG NetCDF files using kelp
annotation shapefiles. For each scene with a kelp annotation:
  1. Extract all pixels inside the polygon
  2. Filter to NDVI > threshold to remove water
  3. Write X, Y, NDVI, <wl_1>, <wl_2>, ... to CSV

Outputs:
  outputs/scenes/<scene_id>/exports/<scene_id>_spectra.csv
  outputs/manifests/spectra_extract_log.csv
"""

from pathlib import Path
import sys

REPO     = Path("Y:/personal/lberberian/BioSCape/BioSCape-AVRISNG_Spectra")
RFL_SUB  = REPO / "data" / "rfl_ocean_subset"
SCENES   = REPO / "outputs" / "scenes"
MANIFEST = REPO / "outputs" / "manifests"

ANNOTATION_GLOB  = "*_Kelp.shp"
NDVI_RED_NM      = 670.0
NDVI_NIR_NM      = 750.0
NDVI_THRESHOLD   = 0.0
VALID_REFL_RANGE = (0.0, 1.5)
OVERWRITE_EXISTING = False


def preflight_checks():
    errors, warnings = [], []

    if not REPO.exists():
        errors.append(f"REPO not found: {REPO}")
        return [], errors, warnings

    if not RFL_SUB.exists():
        errors.append(f"rfl_ocean_subset not found: {RFL_SUB}")
    else:
        nc_count = len(list(RFL_SUB.glob("*.nc")))
        if nc_count == 0:
            errors.append(f"No .nc files in {RFL_SUB}")
        else:
            print(f"  .nc files : {nc_count}")

    if not SCENES.exists():
        errors.append(f"scenes folder not found: {SCENES}")
        return [], errors, warnings

    ann_files = sorted(SCENES.rglob(f"annotations/{ANNOTATION_GLOB}"))
    if not ann_files:
        errors.append(f"No annotations matching '{ANNOTATION_GLOB}' under {SCENES}")
        return [], errors, warnings
    print(f"  annotations : {len(ann_files)}")

    n_nc_missing = n_csv_exists = 0
    for ann_path in ann_files:
        scene_id = ann_path.parent.parent.name
        nc_path = RFL_SUB / f"{scene_id}.nc"
        out_csv = SCENES / scene_id / "exports" / f"{scene_id}_spectra.csv"
        exp_dir = SCENES / scene_id / "exports"

        if not nc_path.exists():
            warnings.append(f"{scene_id} -- .nc not found, will skip")
            n_nc_missing += 1
        if ann_path.stat().st_size == 0:
            errors.append(f"Empty annotation: {ann_path}")
        if ann_path.suffix.lower() == ".shp":
            for ext in (".dbf", ".prj", ".shx"):
                if not ann_path.with_suffix(ext).exists():
                    warnings.append(f"{scene_id} -- missing {ext} sidecar")
        if not exp_dir.exists():
            exp_dir.mkdir(parents=True, exist_ok=True)
        if out_csv.exists():
            n_csv_exists += 1

    if n_nc_missing:
        warnings.append(f"{n_nc_missing} scene(s) will be skipped -- no .nc")
    if n_csv_exists and not OVERWRITE_EXISTING:
        warnings.append(f"{n_csv_exists} CSV(s) already exist (OVERWRITE_EXISTING=False)")

    try:
        import bioscape_rfl_tools  # noqa: F401
    except ImportError:
        errors.append(f"Cannot import bioscape_rfl_tools from {REPO}")

    missing = [p for p in ("netCDF4","geopandas","rasterio","pyproj","shapely","numpy","pandas")
               if not _can_import(p)]
    if missing:
        errors.append(f"Missing packages: {', '.join(missing)}")

    return ann_files, errors, warnings


def _can_import(pkg):
    try:
        __import__(pkg)
        return True
    except ImportError:
        return False


# --- main ---

print("=" * 60)
print("extract_kelp_spectra.py -- pre-flight checks")
print("=" * 60)

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

ann_files, errors, warnings = preflight_checks()

for w in warnings:
    print(f"  WARNING: {w}")
if errors:
    print("\nERRORS:")
    for e in errors:
        print(f"  x {e}")
    sys.exit(1)

print(f"\nReady: {len(ann_files)} scenes, NDVI > {NDVI_THRESHOLD}, "
      f"bands {NDVI_RED_NM}/{NDVI_NIR_NM} nm\n")

if input("Proceed? [y/N] ").strip().lower() != "y":
    print("Aborted.")
    sys.exit(0)

# --- extraction ---

import importlib
import bioscape_rfl_tools as brt
importlib.reload(brt)
import pandas as pd

MANIFEST.mkdir(parents=True, exist_ok=True)
results = []

for ann_path in ann_files:
    scene_id = ann_path.parent.parent.name
    nc_path = RFL_SUB / f"{scene_id}.nc"
    out_csv = SCENES / scene_id / "exports" / f"{scene_id}_spectra.csv"

    print(f"\n{scene_id}")

    if not nc_path.exists():
        print("  SKIP: .nc not found")
        results.append({"scene_id": scene_id, "annotation": str(ann_path),
                        "n_pixels_in_mask": None, "n_pixels_kept": None,
                        "ndvi_threshold": NDVI_THRESHOLD, "out_csv": str(out_csv),
                        "error": ".nc not found"})
        continue

    if out_csv.exists() and not OVERWRITE_EXISTING:
        if input(f"  {out_csv.name} exists -- overwrite? [y/N] ").strip().lower() != "y":
            print("  Skipped.")
            results.append({"scene_id": scene_id, "annotation": str(ann_path),
                            "n_pixels_in_mask": None, "n_pixels_kept": None,
                            "ndvi_threshold": NDVI_THRESHOLD, "out_csv": str(out_csv),
                            "error": "skipped (CSV exists)"})
            continue

    try:
        result = brt.extract_kelp_spectra(
            nc_path, ann_path, out_csv,
            ndvi_wls_nm=(NDVI_RED_NM, NDVI_NIR_NM),
            ndvi_threshold=NDVI_THRESHOLD,
            valid_refl_range=VALID_REFL_RANGE,
        )
        result["annotation"] = str(ann_path)

        if result["error"] is None:
            if not out_csv.exists() or out_csv.stat().st_size == 0:
                result["error"] = "CSV missing or empty after extraction"
            else:
                check = pd.read_csv(out_csv, nrows=2)
                missing_cols = [c for c in ("X","Y","NDVI") if c not in check.columns]
                if missing_cols:
                    result["error"] = f"missing columns: {missing_cols}"

        status = (f"  OK: {result['n_pixels_kept']:,} pixels "
                  f"(of {result['n_pixels_in_mask']:,})" if result["error"] is None
                  else f"  FAIL: {result['error']}")
        print(status)

    except Exception as e:
        import traceback
        result = {"scene_id": scene_id, "annotation": str(ann_path),
                  "n_pixels_in_mask": None, "n_pixels_kept": None,
                  "ndvi_threshold": NDVI_THRESHOLD, "out_csv": str(out_csv),
                  "error": str(e)}
        print(f"  ERROR: {e}")
        traceback.print_exc()

    results.append(result)

# --- summary ---

print(f"\n{'=' * 60}")
results_df = pd.DataFrame(results)
n_ok = int(results_df["error"].isna().sum())
n_fail = int(results_df["error"].notna().sum())

print(f"Succeeded: {n_ok} / {len(results_df)}")
print(f"Failed:    {n_fail}")
if n_ok > 0:
    print(f"Total pixels: {int(results_df['n_pixels_kept'].sum()):,}")

log_path = MANIFEST / "spectra_extract_log.csv"
results_df.to_csv(log_path, index=False)
print(f"Log: {log_path}")

if n_fail > 0:
    print("\nFailed:")
    for _, row in results_df[results_df["error"].notna()].iterrows():
        print(f"  {row['scene_id']} -- {row['error']}")
