# Seafloor spreading magnetic test cases

These inputs provide progressively more complex finite-volume magnetic models for an Earth-radius sphere. Run them with `rsphere_km=6371.2`.

## Cases

### `st_paul_reference_body.in`

A single rectangular polygon based on the synthetic validation body in Sichler et al. (2002):

- 4 km across strike
- 10 km along strike
- 1 km thick
- 10 A/m magnetization

The body is centered on the equator and placed from 2.5 to 3.5 km depth.

### `st_paul_symmetric_stripes.in`

A symmetric 16-body block model using the St. Paul Fracture Zone full spreading rate of 31.47 km/Myr and a 0.5 km thick source layer. Each flank therefore grows at 15.735 km/Myr. Stripe boundaries use the following geomagnetic polarity ages:

| boundary | age [Ma] | distance from axis [km] | distance from axis [deg] |
|---|---:|---:|---:|
| Brunhes base | 0.781 | 12.287 | 0.1104 |
| Jaramillo top | 0.988 | 15.545 | 0.1397 |
| Jaramillo base | 1.072 | 16.867 | 0.1515 |
| Cobb Mountain top | 1.173 | 18.456 | 0.1658 |
| Cobb Mountain base | 1.185 | 18.646 | 0.1675 |
| Olduvai top | 1.778 | 27.977 | 0.2513 |
| Olduvai base | 1.945 | 30.603 | 0.2749 |
| Matuyama/Gauss | 2.595 | 40.831 | 0.3668 |

The unusually narrow Cobb Mountain bodies deliberately test thin-volume handling. Normal and reversed blocks use equal 10 A/m amplitudes with opposite global Cartesian inclinations.

### `st_paul_segmented_transform_alteration.in`

An 18-body stress case combining:

- two 12 km long ridge segments
- a 0.08 degree right step across a transform gap
- symmetric normal and reversed stripes through the Jaramillo interval
- an irregular transferred sliver representing ridge propagation or a ridge jump
- an overlapping 8 A/m opposite-polarity correction body that reduces a local normal stripe from 10 to 2 A/m

The 4 km wide correction is inspired by the reduced-magnetization zone at the TAG hydrothermal field. Its 14-point closed outline contains 13 unique vertices and activates the solver's complex-polygon classification at the default 12-vertex threshold.

## References

- Sichler et al. (2002), St. Paul Fracture Zone model: <https://doi.org/10.1029/2001JB000401>
- IODP Expedition 318 geomagnetic polarity table: <https://publications.iodp.org/proceedings/318/102/102_t5.htm>
- Tivey et al. (2003), TAG reduced-magnetization zone: <https://doi.org/10.1029/2002JB001967>
- van Hinsbergen et al. (2013), transferred blocks from ridge jumps and propagation: <https://doi.org/10.1093/GJI/GGT162>

The source papers constrain the dimensions and magnetic contrasts, but these files are synthetic solver tests rather than reconstructions of the cited survey grids. Constant depth replaces the bathymetry-following upper surfaces used in the published models.

## Running

From `examples`:

```bash
./build_gravmag_tools.sh

./gravmag_sphere_point_source 6371.2 \
  seafloor_spreading_cases/st_paul_reference_body.in \
  seafloor_spreading_cases/output/st_paul_reference_body_xyz.txt 1

./gravmag_sphere_point_source 6371.2 \
  seafloor_spreading_cases/st_paul_symmetric_stripes.in \
  seafloor_spreading_cases/output/st_paul_symmetric_stripes_xyz.txt 1

./gravmag_sphere_spectral 6371.2 \
  seafloor_spreading_cases/st_paul_segmented_transform_alteration.in \
  seafloor_spreading_cases/output/st_paul_segmented_transform_alteration_xyz.txt
```

The point_source solver is the clearest reference for the narrow bodies. The segmented case is intended to exercise collective multi-body accumulation, polygon-edge correction, and automatic hybrid evaluation.

## Plotting

`plot_bxyz.py` reads the source geometry from a `.in` file, sums `Bx`, `By`, and `Bz` over repeated body ids in a point_source-solver output, recomputes `Btot` from the summed vector, and overlays the source outlines on four field panels. It also accepts collective spectral output, which already contains one row per observation point.

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python \
  seafloor_spreading_cases/plot_bxyz.py \
  seafloor_spreading_cases/st_paul_symmetric_stripes.in \
  seafloor_spreading_cases/output/st_paul_symmetric_stripes_xyz.txt \
  seafloor_spreading_cases/figs/st_paul_symmetric_stripes_bxyz.png
```

The image path is optional. When omitted, the script writes `<output_stem>_bxyz.png` in the source `.in` directory's `figs/` subfolder. Numeric solver defaults use that same input directory's `output/` folder. Explicit paths override defaults; see [output path conventions](../../docs/output_paths.md).
