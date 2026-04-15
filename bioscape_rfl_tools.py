"""
bioscape_rfl_tools.py

Read BioSCape AVIRIS-NG L2A RFL ORT NetCDF files. Build bbox footprints,
intersect with coastline, stage selected files, export quicklooks (RGB + NDVI)
as GeoTIFFs, and extract per-pixel kelp spectra.

NetCDF structure:
  - group: reflectance
      - reflectance (wavelength, northing, easting)
      - wavelength (wavelength)
  - variables: easting, northing
  - grid mapping: transverse_mercator

Hemisphere fix: AVIRIS-NG L2A OE files store false_northing=0 for southern
hemisphere scenes. load_spec() detects this using the Y coordinates and
corrects UTM-North (326XX) to UTM-South (327XX).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import netCDF4 as nc
from pyproj import CRS, Transformer
from shapely.geometry import Polygon

RFL_VAR_PATH = "reflectance/reflectance"
WL_VAR_PATH = "reflectance/wavelength"
X_VAR_NAME = "easting"
Y_VAR_NAME = "northing"
GRIDMAP_VAR = "transverse_mercator"

# UTM-South scenes for SA have Y ~ 6.1M-6.4M due to the 10M false northing.
# Any UTM-North scene at the same latitude would have Y < 4M.
_UTM_SOUTH_Y_THRESHOLD = 4_000_000.0


@dataclass(frozen=True)
class CubeSpec:
    dims: Tuple[str, ...]
    shape: Tuple[int, ...]
    spec_axis: int
    y_axis: int
    x_axis: int
    wl: np.ndarray
    x: np.ndarray
    y: np.ndarray
    crs: Optional[CRS]
    epsg: Optional[int]


# ---------------------------------------------------------------------------
# NetCDF helpers
# ---------------------------------------------------------------------------

def _get_var_by_path(ds: nc.Dataset, var_path: str) -> nc.Variable:
    parts = var_path.split("/")
    grp: Any = ds
    for p in parts[:-1]:
        grp = grp.groups[p]
    return grp.variables[parts[-1]]


def _to_plain_dict(var: nc.Variable) -> Dict[str, Any]:
    d: Dict[str, Any] = {}
    for k in var.ncattrs():
        v = var.getncattr(k)
        if hasattr(v, "item"):
            try:
                v = v.item()
            except Exception:
                pass
        d[k] = v
    return d


# ---------------------------------------------------------------------------
# CRS resolution
# ---------------------------------------------------------------------------

def get_crs_from_grid_mapping(ds: nc.Dataset) -> Tuple[Optional[CRS], Optional[int]]:
    """Best-effort CRS from the transverse_mercator variable.
    Tries: (1) WKT, (2) CF dict, (3) UTM fallback from central meridian.
    Hemisphere correction is deferred to load_spec()."""
    if GRIDMAP_VAR not in ds.variables:
        return None, None

    gm = ds.variables[GRIDMAP_VAR]

    # 1) WKT
    for k in ("spatial_ref", "crs_wkt", "wkt"):
        if k in gm.ncattrs():
            try:
                crs = CRS.from_wkt(gm.getncattr(k))
                return crs, crs.to_epsg()
            except Exception:
                pass

    attrs = _to_plain_dict(gm)

    # 2) CF projection
    try:
        crs = CRS.from_cf(attrs)
        return crs, crs.to_epsg()
    except Exception:
        pass

    # 3) UTM fallback
    lon0 = attrs.get("longitude_of_central_meridian")
    k0 = attrs.get("scale_factor_at_central_meridian")
    fe = attrs.get("false_easting")
    fn = attrs.get("false_northing")

    try:
        lon0 = float(lon0) if lon0 is not None else None
        k0 = float(k0) if k0 is not None else None
        fe = float(fe) if fe is not None else None
        fn = float(fn) if fn is not None else None
    except Exception:
        return None, None

    if lon0 is None:
        return None, None

    if (k0 is not None and abs(k0 - 0.9996) < 1e-3) and \
       (fe is not None and abs(fe - 500000.0) < 2e4):
        zone = int(round((lon0 + 183.0) / 6.0))
        is_south = (fn is not None and abs(fn - 10_000_000) < 1e4)
        epsg = (32700 if is_south else 32600) + zone
        try:
            return CRS.from_epsg(epsg), epsg
        except Exception:
            return None, epsg

    return None, None


def _correct_utm_hemisphere(
    crs: Optional[CRS], epsg: Optional[int], y: np.ndarray
) -> Tuple[Optional[CRS], Optional[int]]:
    """Swap UTM-North (326XX) to UTM-South (327XX) when Y coords indicate
    a southern-hemisphere scene with false_northing incorrectly set to 0."""
    if epsg is None or crs is None:
        return crs, epsg
    if not (32600 <= epsg <= 32660):
        return crs, epsg
    if float(np.nanmedian(y)) <= _UTM_SOUTH_Y_THRESHOLD:
        return crs, epsg
    corrected = epsg + 100
    try:
        return CRS.from_epsg(corrected), corrected
    except Exception:
        return crs, epsg


# ---------------------------------------------------------------------------
# Cube metadata
# ---------------------------------------------------------------------------

def load_spec(ds: nc.Dataset) -> CubeSpec:
    var = _get_var_by_path(ds, RFL_VAR_PATH)
    wl = np.array(_get_var_by_path(ds, WL_VAR_PATH)[:]).astype(float)
    dims = tuple(var.dimensions)
    shape = tuple(var.shape)

    for need in ("wavelength", "northing", "easting"):
        if need not in dims:
            raise ValueError(f"missing dimension '{need}' in {dims}")

    spec_axis = dims.index("wavelength")
    y_axis = dims.index("northing")
    x_axis = dims.index("easting")

    x = np.array(ds.variables[X_VAR_NAME][:]).astype(float)
    y = np.array(ds.variables[Y_VAR_NAME][:]).astype(float)

    crs, epsg = get_crs_from_grid_mapping(ds)
    crs, epsg = _correct_utm_hemisphere(crs, epsg, y)

    return CubeSpec(dims=dims, shape=shape, spec_axis=spec_axis,
                    y_axis=y_axis, x_axis=x_axis, wl=wl, x=x, y=y,
                    crs=crs, epsg=epsg)


# ---------------------------------------------------------------------------
# Band I/O
# ---------------------------------------------------------------------------

def list_nc_files(nc_all_dir: Path, pattern: str = "*_RFL_ORT.nc") -> List[Path]:
    return sorted(Path(nc_all_dir).rglob(pattern))


def closest_band(wl: np.ndarray, target_nm: float) -> int:
    return int(np.nanargmin(np.abs(wl - target_nm)))


def read_band(ds: nc.Dataset, spec: CubeSpec, band_index: int) -> np.ndarray:
    var = _get_var_by_path(ds, RFL_VAR_PATH)
    sl = [slice(None)] * var.ndim
    sl[spec.spec_axis] = band_index
    arr = np.array(var[tuple(sl)]).astype("float32").squeeze()
    expected = (spec.shape[spec.y_axis], spec.shape[spec.x_axis])
    if arr.shape != expected:
        raise ValueError(f"band shape {arr.shape} != expected {expected}")
    return arr


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def bbox_from_xy(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float, float]:
    return (float(np.nanmin(x)), float(np.nanmax(x)),
            float(np.nanmin(y)), float(np.nanmax(y)))


def bbox_to_polygon_in_crs(
    xmin: float, xmax: float, ymin: float, ymax: float,
    src_crs: CRS, dst_crs: CRS,
) -> Polygon:
    t = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
    xs = [xmin, xmax, xmax, xmin, xmin]
    ys = [ymin, ymin, ymax, ymax, ymin]
    X, Y = t.transform(xs, ys)
    return Polygon(zip(X, Y))


def read_coastline_union(coastline_shp: Path):
    import geopandas as gpd
    coast = gpd.read_file(coastline_shp)
    if coast.crs is None:
        raise ValueError("coastline shapefile has no CRS")
    coast_union = coast.geometry.union_all()
    return coast, coast_union, coast.crs


# ---------------------------------------------------------------------------
# File summaries
# ---------------------------------------------------------------------------

def summarize_file_bbox(fn: Path) -> Dict[str, Any]:
    with nc.Dataset(fn, "r") as ds:
        spec = load_spec(ds)
        xmin, xmax, ymin, ymax = bbox_from_xy(spec.x, spec.y)
        return {
            "file": Path(fn).name, "path": str(fn), "epsg_guess": spec.epsg,
            "xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax,
            "shape": str(spec.shape), "dims": ",".join(spec.dims),
            "n_bands": int(spec.shape[spec.spec_axis]),
            "wl_min_nm": float(np.nanmin(spec.wl)),
            "wl_max_nm": float(np.nanmax(spec.wl)),
        }


def summarize_file_bbox_and_intersection(fn: Path, coast_union, coast_crs) -> Dict[str, Any]:
    with nc.Dataset(fn, "r") as ds:
        spec = load_spec(ds)
        if spec.crs is None:
            raise ValueError("could not read CRS from grid mapping")
        xmin, xmax, ymin, ymax = bbox_from_xy(spec.x, spec.y)
        footprint = bbox_to_polygon_in_crs(xmin, xmax, ymin, ymax, spec.crs, coast_crs)
        return {
            "file": Path(fn).name, "path": str(fn), "epsg_guess": spec.epsg,
            "xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax,
            "footprint": footprint,
            "intersects_coast": bool(footprint.intersects(coast_union)),
        }


# ---------------------------------------------------------------------------
# File staging
# ---------------------------------------------------------------------------

def stage_files(src_paths: Sequence[str], dst_dir: Path, mode: str = "hardlink") -> None:
    import os, shutil
    dst_dir = Path(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)
    n_ok = n_skip = n_fail = 0
    for src in src_paths:
        src = Path(src)
        dst = dst_dir / src.name
        if dst.exists():
            n_skip += 1
            continue
        try:
            if mode == "hardlink":
                os.link(src, dst)
            elif mode == "copy":
                shutil.copy2(src, dst)
            else:
                raise ValueError("mode must be 'hardlink' or 'copy'")
            n_ok += 1
        except Exception as e:
            n_fail += 1
            print("failed:", src.name, "|", e)
    print(f"staged -> {dst_dir}  ok:{n_ok} skip:{n_skip} fail:{n_fail}")


# ---------------------------------------------------------------------------
# Array manipulation helpers
# ---------------------------------------------------------------------------

def _fix_xy_orientation(x, y, arr_yx):
    """Make x increasing, y decreasing (row 0 = north)."""
    x2, y2, a2 = x.copy(), y.copy(), arr_yx
    if x2.size >= 2 and x2[1] < x2[0]:
        x2 = x2[::-1]
        a2 = a2[:, ::-1]
    if y2.size >= 2 and y2[1] > y2[0]:
        y2 = y2[::-1]
        a2 = a2[::-1, :]
    return x2, y2, a2


def _transform_from_xy(x, y):
    from rasterio.transform import from_origin
    if x.size < 2 or y.size < 2:
        raise ValueError("x/y vectors too small for transform")
    dx = abs(float(np.nanmedian(np.diff(x))))
    dy = abs(float(np.nanmedian(np.diff(y))))
    west = float(x[0]) - dx / 2.0
    north = float(y[0]) + dy / 2.0
    return from_origin(west, north, dx, dy)


def _as_float32(a):
    if isinstance(a, np.ma.MaskedArray):
        a = a.filled(np.nan)
    return np.array(a, dtype="float32")


def _mask_invalid(a, vmin=None, vmax=None):
    out = a.astype("float32", copy=True)
    m = np.isfinite(out)
    if vmin is not None:
        m &= out >= float(vmin)
    if vmax is not None:
        m &= out <= float(vmax)
    out[~m] = np.nan
    return out


def _stretch_percentile(a, pmin=2.0, pmax=98.0, vmin=None, vmax=None):
    v = a[np.isfinite(a)]
    if vmin is not None:
        v = v[v >= float(vmin)]
    if vmax is not None:
        v = v[v <= float(vmax)]
    if v.size == 0:
        return np.zeros_like(a, dtype="float32")
    lo, hi = np.nanpercentile(v, [pmin, pmax])
    if not np.isfinite(lo) or not np.isfinite(hi) or (hi - lo) <= 1e-12:
        return np.zeros_like(a, dtype="float32")
    out = (a - lo) / (hi - lo)
    out = np.clip(out, 0, 1).astype("float32")
    out[~np.isfinite(a)] = np.nan
    return out


def _scale_to_u8(a01):
    out = a01.copy()
    out[~np.isfinite(out)] = 0.0
    return (np.clip(out, 0, 1) * 255.0).astype("uint8")


# ---------------------------------------------------------------------------
# Quicklook export
# ---------------------------------------------------------------------------

def export_scene_quicklooks(
    nc_path: Path,
    scenes_root: Path,
    rgb_wls_nm=(650.0, 560.0, 470.0),
    ndvi_wls_nm=(670.0, 750.0),
    rgb_stretch_p=(2.0, 98.0),
    rgb_valid_range=(0.0, 1.5),
    ndvi_nodata=-9999.0,
    write_ndvi_u8=True,
    ndvi_u8_mode="fixed",
    ndvi_u8_range=(-0.2, 0.8),
    ndvi_u8_stretch_p=(2.0, 98.0),
) -> Dict[str, Any]:
    """Export RGB u8, NDVI float32, and optional NDVI u8 GeoTIFFs for one scene."""
    import json, rasterio

    nc_path = Path(nc_path)
    scene_id = nc_path.stem
    scene_dir = Path(scenes_root) / scene_id
    qdir = scene_dir / "quicklooks"
    for p in (qdir, scene_dir / "annotations", scene_dir / "exports", scene_dir / "meta"):
        p.mkdir(parents=True, exist_ok=True)

    with nc.Dataset(nc_path, "r") as ds:
        spec = load_spec(ds)
        if spec.crs is None:
            raise ValueError("could not read CRS from NetCDF")

        # read bands
        r_i = closest_band(spec.wl, rgb_wls_nm[0])
        g_i = closest_band(spec.wl, rgb_wls_nm[1])
        b_i = closest_band(spec.wl, rgb_wls_nm[2])
        nd_r_i = closest_band(spec.wl, ndvi_wls_nm[0])
        nd_n_i = closest_band(spec.wl, ndvi_wls_nm[1])

        r = _as_float32(read_band(ds, spec, r_i))
        g = _as_float32(read_band(ds, spec, g_i))
        b = _as_float32(read_band(ds, spec, b_i))
        red = _as_float32(read_band(ds, spec, nd_r_i))
        nir = _as_float32(read_band(ds, spec, nd_n_i))

        # fix orientation
        x_fix, y_fix, red_fix = _fix_xy_orientation(spec.x, spec.y, red)
        _, _, nir_fix = _fix_xy_orientation(spec.x, spec.y, nir)
        _, _, r_fix = _fix_xy_orientation(spec.x, spec.y, r)
        _, _, g_fix = _fix_xy_orientation(spec.x, spec.y, g)
        _, _, b_fix = _fix_xy_orientation(spec.x, spec.y, b)

        transform = _transform_from_xy(x_fix, y_fix)

        # mask invalid reflectance
        r_fix = _mask_invalid(r_fix, *rgb_valid_range)
        g_fix = _mask_invalid(g_fix, *rgb_valid_range)
        b_fix = _mask_invalid(b_fix, *rgb_valid_range)
        red_fix = _mask_invalid(red_fix, *rgb_valid_range)
        nir_fix = _mask_invalid(nir_fix, *rgb_valid_range)

        # NDVI
        ndvi = (nir_fix - red_fix) / (nir_fix + red_fix + 1e-12)
        ndvi = _mask_invalid(ndvi, -1.0, 1.0)

        # RGB stretch
        pmin, pmax = rgb_stretch_p
        r01 = _stretch_percentile(r_fix, pmin, pmax, *rgb_valid_range)
        g01 = _stretch_percentile(g_fix, pmin, pmax, *rgb_valid_range)
        b01 = _stretch_percentile(b_fix, pmin, pmax, *rgb_valid_range)
        rgb_u8 = np.stack([_scale_to_u8(r01), _scale_to_u8(g01), _scale_to_u8(b01)])

        tif_kwargs = dict(driver="GTiff", crs=spec.crs.to_wkt(), transform=transform,
                          compress="deflate", tiled=True, BIGTIFF="IF_SAFER")

        # write NDVI float
        ndvi_path = qdir / f"{scene_id}_ndvi.tif"
        ndvi_out = ndvi.copy()
        ndvi_out[~np.isfinite(ndvi_out)] = float(ndvi_nodata)
        with rasterio.open(ndvi_path, "w", height=ndvi_out.shape[0], width=ndvi_out.shape[1],
                           count=1, dtype="float32", nodata=float(ndvi_nodata), **tif_kwargs) as dst:
            dst.write(ndvi_out, 1)

        # write NDVI u8
        ndvi_u8_path = None
        if write_ndvi_u8:
            if ndvi_u8_mode == "fixed":
                vmin, vmax = ndvi_u8_range
                nd01 = np.clip((ndvi - vmin) / (vmax - vmin + 1e-12), 0, 1).astype("float32")
            elif ndvi_u8_mode == "percentile":
                nd01 = _stretch_percentile(ndvi, *ndvi_u8_stretch_p, vmin=-1.0, vmax=1.0)
            else:
                raise ValueError("ndvi_u8_mode must be 'fixed' or 'percentile'")
            ndvi_u8_path = qdir / f"{scene_id}_ndvi_u8.tif"
            with rasterio.open(ndvi_u8_path, "w", height=nd01.shape[0], width=nd01.shape[1],
                               count=1, dtype="uint8", **tif_kwargs) as dst:
                dst.write(_scale_to_u8(nd01), 1)

        # write RGB
        rgb_path = qdir / f"{scene_id}_rgb_u8.tif"
        with rasterio.open(rgb_path, "w", height=rgb_u8.shape[1], width=rgb_u8.shape[2],
                           count=3, dtype="uint8", photometric="RGB", **tif_kwargs) as dst:
            dst.write(rgb_u8)

        meta = {
            "scene_id": scene_id, "nc_path": str(nc_path), "crs_epsg": spec.epsg,
            "rgb_targets_nm": list(rgb_wls_nm),
            "rgb_wls_nm": [float(spec.wl[r_i]), float(spec.wl[g_i]), float(spec.wl[b_i])],
            "rgb_band_indices": [int(r_i), int(g_i), int(b_i)],
            "ndvi_targets_nm": list(ndvi_wls_nm),
            "ndvi_wls_nm": [float(spec.wl[nd_r_i]), float(spec.wl[nd_n_i])],
            "ndvi_band_indices": [int(nd_r_i), int(nd_n_i)],
            "outputs": {
                "rgb_u8": str(rgb_path), "ndvi": str(ndvi_path),
                "ndvi_u8": str(ndvi_u8_path) if ndvi_u8_path else None,
            },
        }
        (scene_dir / "meta" / f"{scene_id}_bands.json").write_text(json.dumps(meta, indent=2))

    return meta


def batch_export_quicklooks(nc_paths: Iterable[Path], scenes_root: Path, **kwargs) -> List[Dict[str, Any]]:
    metas = []
    for p in nc_paths:
        try:
            metas.append(export_scene_quicklooks(p, scenes_root, **kwargs))
        except Exception as e:
            metas.append({"scene_id": Path(p).stem, "nc_path": str(p), "error": str(e)})
            print("failed:", Path(p).name, "|", e)
    return metas


# ---------------------------------------------------------------------------
# Coastline intersection
# ---------------------------------------------------------------------------

def build_coast_intersection_gdf(nc_files: Sequence[Path], coastline_shp: Path):
    import pandas as pd
    import geopandas as gpd

    coast, coast_union, coast_crs = read_coastline_union(coastline_shp)
    rows = []
    for fn in nc_files:
        try:
            rows.append(summarize_file_bbox_and_intersection(fn, coast_union, coast_crs))
        except Exception as e:
            rows.append({"file": Path(fn).name, "path": str(fn), "error": str(e)})

    df = pd.DataFrame(rows)
    ok = df["error"].isna() if "error" in df.columns else np.ones(len(df), dtype=bool)
    foot_gdf = gpd.GeoDataFrame(df.loc[ok].copy(), geometry="footprint", crs=coast_crs)
    return foot_gdf, coast


# ---------------------------------------------------------------------------
# Kelp spectra extraction
# ---------------------------------------------------------------------------

def extract_kelp_spectra(
    nc_path: Path,
    annotation_shp: Path,
    out_csv: Path,
    ndvi_wls_nm=(670.0, 750.0),
    ndvi_threshold=0.0,
    valid_refl_range=(0.0, 1.5),
) -> Dict[str, Any]:
    """Extract per-pixel spectra inside a kelp annotation polygon.
    Output CSV columns: X, Y, NDVI, <wl_1_nm>, <wl_2_nm>, ..."""
    import geopandas as gpd
    import pandas as pd
    from rasterio.features import geometry_mask

    nc_path = Path(nc_path)
    out_csv = Path(out_csv)
    scene_id = nc_path.stem
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    with nc.Dataset(nc_path, "r") as ds:
        spec = load_spec(ds)
        if spec.crs is None:
            raise ValueError("could not resolve CRS")

        # load and reproject annotation
        ann = gpd.read_file(annotation_shp)
        if ann.crs is None:
            raise ValueError(f"annotation has no CRS: {annotation_shp}")
        ann = ann.to_crs(spec.crs)

        # fix orientation
        dummy = np.ones((spec.shape[spec.y_axis], spec.shape[spec.x_axis]), dtype="float32")
        x_fix, y_fix, _ = _fix_xy_orientation(spec.x, spec.y, dummy)
        transform = _transform_from_xy(x_fix, y_fix)
        n_rows, n_cols = len(y_fix), len(x_fix)

        # rasterize annotation
        inside = ~geometry_mask(
            [geom.__geo_interface__ for geom in ann.geometry if geom is not None],
            out_shape=(n_rows, n_cols), transform=transform, invert=False,
        )
        n_in_mask = int(inside.sum())
        if n_in_mask == 0:
            return {"scene_id": scene_id, "n_pixels_in_mask": 0, "n_pixels_kept": 0,
                    "out_csv": str(out_csv), "error": "annotation does not overlap scene"}

        # NDVI
        red_i = closest_band(spec.wl, ndvi_wls_nm[0])
        nir_i = closest_band(spec.wl, ndvi_wls_nm[1])
        red = _mask_invalid(_as_float32(read_band(ds, spec, red_i)), *valid_refl_range)
        nir = _mask_invalid(_as_float32(read_band(ds, spec, nir_i)), *valid_refl_range)
        _, _, red = _fix_xy_orientation(spec.x, spec.y, red)
        _, _, nir = _fix_xy_orientation(spec.x, spec.y, nir)
        ndvi_arr = _mask_invalid((nir - red) / (nir + red + 1e-12), -1.0, 1.0)

        # pixel selection
        keep = inside & np.isfinite(ndvi_arr) & (ndvi_arr > ndvi_threshold)
        row_idx, col_idx = np.where(keep)
        if len(row_idx) == 0:
            return {"scene_id": scene_id, "n_pixels_in_mask": n_in_mask, "n_pixels_kept": 0,
                    "out_csv": str(out_csv), "error": f"no pixels pass NDVI > {ndvi_threshold}"}

        x_coords = x_fix[col_idx]
        y_coords = y_fix[row_idx]
        ndvi_vals = ndvi_arr[row_idx, col_idx]

        # read all bands for kept pixels
        n_bands = spec.shape[spec.spec_axis]
        band_arrays = np.full((len(row_idx), n_bands), np.nan, dtype="float32")
        for b in range(n_bands):
            band = _mask_invalid(_as_float32(read_band(ds, spec, b)), *valid_refl_range)
            _, _, band = _fix_xy_orientation(spec.x, spec.y, band)
            band_arrays[:, b] = band[row_idx, col_idx]

        # write CSV
        wl_cols = [f"{float(w):.5f}" for w in spec.wl]
        df = pd.DataFrame(band_arrays, columns=wl_cols)
        df.insert(0, "NDVI", ndvi_vals.astype("float64"))
        df.insert(0, "Y", y_coords)
        df.insert(0, "X", x_coords)
        df.to_csv(out_csv, index=False)

    return {"scene_id": scene_id, "n_pixels_in_mask": n_in_mask,
            "n_pixels_kept": len(row_idx), "ndvi_threshold": ndvi_threshold,
            "out_csv": str(out_csv), "error": None}


def batch_extract_kelp_spectra(
    scenes_root: Path, rfl_sub: Path,
    annotation_glob="*_Kelp.shp", ndvi_wls_nm=(670.0, 750.0),
    ndvi_threshold=0.0, valid_refl_range=(0.0, 1.5),
) -> List[Dict[str, Any]]:
    import pandas as pd
    scenes_root, rfl_sub = Path(scenes_root), Path(rfl_sub)
    results = []
    ann_files = sorted(scenes_root.rglob(f"annotations/{annotation_glob}"))
    if not ann_files:
        print(f"no annotations matching '{annotation_glob}' under {scenes_root}")
        return results
    print(f"found {len(ann_files)} annotation(s)")
    for ann_path in ann_files:
        scene_id = ann_path.parent.parent.name
        nc_path = rfl_sub / f"{scene_id}.nc"
        out_csv = scenes_root / scene_id / "exports" / f"{scene_id}_spectra.csv"
        if not nc_path.exists():
            print(f"  SKIP  {scene_id} | .nc not found")
            results.append({"scene_id": scene_id, "out_csv": str(out_csv), "error": ".nc not found"})
            continue
        try:
            result = extract_kelp_spectra(nc_path, ann_path, out_csv,
                                          ndvi_wls_nm=ndvi_wls_nm, ndvi_threshold=ndvi_threshold,
                                          valid_refl_range=valid_refl_range)
            status = f"{result['n_pixels_kept']} pixels" if result["error"] is None else f"ERROR: {result['error']}"
            print(f"  {'OK  ' if result['error'] is None else 'FAIL'} {scene_id} | {status}")
        except Exception as e:
            result = {"scene_id": scene_id, "out_csv": str(out_csv), "error": str(e)}
            print(f"  FAIL {scene_id} | {e}")
        results.append(result)
    return results


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def write_manifest(keep_paths: Sequence[str], out_csv: Path, out_txt: Path) -> None:
    import pandas as pd
    out_csv, out_txt = Path(out_csv), Path(out_txt)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"path": list(keep_paths)}).to_csv(out_csv, index=False)
    out_txt.write_text("\n".join(keep_paths) + "\n")
    print("wrote:", out_csv)
    print("wrote:", out_txt)
