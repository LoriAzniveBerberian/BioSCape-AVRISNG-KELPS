# BioSCape-AVRISNG-KELPS

**Kelp Extraction from L2 Pixel Spectra**

This repository processes AVIRIS-NG L2A hyperspectral imagery from the [BioSCape](https://www.bioscape.io/data) campaign to identify and stage coastal scenes along South Africa's coastline. It exports georeferenced quicklook products (RGB and NDVI GeoTIFFs) and extracts per-pixel reflectance spectra from kelp annotation polygons for downstream spectral analysis.


![description](11.jpg)


## Folder Structure

```
BioSCape-AVRISNG_Spectra/
  data/
    rfl_nc/                          # all staged .nc files
    rfl_ocean_subset/                # coastal/ocean .nc files (working dataset)
  SouthAfricaCoastlineMask/          # coastline shapefile for scene selection
  outputs/
    manifests/                       # processing logs and file lists
    scenes/<scene_id>/
      quicklooks/                    # RGB and NDVI GeoTIFFs
      annotations/                   # kelp annotation shapefiles
      exports/                       # extracted spectra CSVs
      meta/                          # band metadata JSON
    mosaics/
    figures/
    spectra/
```

## Files

| File | Description |
|---|---|
| `BioSCape_AVIRISNG_SpectralAnalysis.ipynb` | Main workflow notebook: scene discovery, coastline filtering, staging, quicklook export, kelp spectra extraction, and visualization |
| `bioscape_rfl_tools.py` | Core library for reading AVIRIS-NG L2A NetCDF files, CRS handling (including UTM hemisphere correction), GeoTIFF export, and spectra extraction |
| `extract_kelp_spectra.py` | Standalone script for batch kelp spectra extraction with pre-flight checks |
| `qc_crs_transforms_metadata.ipynb` | QC notebook: validates CRS, GeoTransforms, pixel values, 8-bit stretch, and processing metadata across all scenes |

## Input Data

This repository expects AVIRIS-NG L2A orthorectified surface reflectance NetCDF files (`*_RFL_ORT.nc`) from the BioSCape campaign. These files contain 425 spectral bands (380-2510 nm at 5 nm intervals) with the reflectance cube stored as `reflectance/reflectance` (wavelength, northing, easting) and coordinates in UTM projection.

The BioSCape (Biodiversity Survey of the Cape) campaign flew AVIRIS-NG over South Africa's coastal and terrestrial ecosystems in October-November 2023. Flight data can be accessed through the [BioSCape Data Portal](https://www.bioscape.io/data).

### South Africa Coastline Shapefile

Scene selection uses a coastline shapefile to identify flight lines that intersect the ocean/coastal zone. The coastline data is from the Global Map dataset produced by South Africa's National Geo-spatial Information (NGI).

> International Steering Committee for Global Mapping & South Africa National Geo-spatial Information. (2016). Coasts, South Africa, 2016 [Shapefile]. Stanford Digital Repository. https://purl.stanford.edu/zj111gb9121

### AVIRIS-NG L2 Data Citation

> Green, R. O., Brodrick, P. G., Chapman, J. W., Eastwood, M., Geier, S., Helmlinger, M., Lundeen, S. R., Olson-Duvall, W., Pavlick, R., Rios, L. M., Thompson, D. R., & Thorpe, A. K. (2023). AVIRIS-NG L2 Surface Reflectance, Facility Instrument Collection, V1 (Version 1). ORNL Distributed Active Archive Center. https://doi.org/10.3334/ORNLDAAC/2110

### AVIRIS-NG Instrument Reference

> Chapman, J. W., Thompson, D. R., Helmlinger, M. C., Bue, B. D., Green, R. O., Eastwood, M. L., Geier, S., Olson-Duvall, W., & Lundeen, S. R. (2019). Spectral and radiometric calibration of the Next Generation Airborne Visible Infrared Spectrometer (AVIRIS-NG). Remote Sensing, 11(18), 2129. https://doi.org/10.3390/rs11182129

## Workflow

1. **Discover** scan a folder of AVIRIS-NG L2A `.nc` files and summarize spatial metadata
2. **Filter** intersect flight line footprints with a South Africa coastline shapefile to select ocean/coastal scenes
3. **Stage** hardlink or copy selected `.nc` files into `data/rfl_ocean_subset/`
4. **Export quicklooks** for each scene, write 8-bit RGB, float32 NDVI, and 8-bit NDVI GeoTIFFs with correct UTM-South CRS
5. **Annotate** manually draw kelp polygons in ArcGIS/QGIS and save to each scene's `annotations/` folder
6. **Extract spectra** for each annotated scene, rasterize the polygon, filter by NDVI threshold, and write per-pixel spectra to CSV (columns: X, Y, NDVI, wavelength_1, wavelength_2, ...)

## Known Issues

Four scenes from flight line `ang20231029t104631` (segments `_000` through `_003`) were processed by JPL with a development software build (`software_build_version: 010200_rdndev`, `product_version: test`) and exhibit a georeferencing offset relative to basemap imagery. The remaining 46 scenes are production data (`software_build_version: 002`, `product_version: 1`). This is an upstream data issue and does not affect the processing code. See the QC notebook for details.

## Requirements

- Python 3.10+
- netCDF4
- numpy
- pandas
- geopandas
- rasterio
- pyproj
- shapely
- folium
- matplotlib
- Pillow
