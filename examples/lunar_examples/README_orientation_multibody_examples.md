# Orientation + Multi-body Example Inputs

Gravity and magnetic input files for the GravMag Sphere point_source and spectral solvers, using a lunar reference radius of 1737.4 km. Numeric products use this directory's `output/` and figures use `figs/`.

All files follow the naming convention: 
`gravmag_sphere_{numbod}body_{ifield}_{nblim}_inc{incdeg}_dec{decdeg}_{opt:desc}`  
where

- `numbod` = total number of sources defined
- `ifield` = one of (mag, grav), describes field controlled by `ifield` flag
- `nblim` = one of (fixedlim, polygon), describes source geometry controlled by `nblim` flag
- `incdeg` = elevation of magnetization above the global Cartesian xy plane toward +z
- `decdeg` = azimuth of magnetization in the global Cartesian xy plane from +x toward +y
- `opt:desc` = optional further description of example, i.e., weak magnetization, complex large polygon

The frame is body-fixed: +x points to latitude/longitude zero, +y to 90 degrees east on the equator, and +z to the north pole. These are not local north/east/down angles. See [coordinate conventions](../../diagnostics/gravmagsphere_inclination_declination_physical_audit.md).

## Magnetic 1-body pure orientation examples

- `gravmag_sphere_1body_mag_fixedlim_inc90_dec0_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc0_dec0_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc0_dec90_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc-45_dec120_orient.in`

All use:

- magnetic mode (`ifield=2`)
- fixed body limits (`nblim=1`)
- same geometry/depth, only magnetization orientation differs (`incdeg`, `decdeg`)

## Magnetic multi-body polygon examples

- `gravmag_sphere_3body_mag_polygon_incmix_decmix.in`
- `gravmag_sphere_5body_mag_polygon_incmix_decmix.in`

Both use:

- repeated body blocks (e.g., `numbod=3`)
- polygon geometry (`nblim=0`)
- mixed orientations and depths

## Magnetic 1-body complex non-rectangular polygon examples

- `gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in`

This case uses:

- one large irregular polygon with 25 vertices (not rectangular)
- bulk magnetization direction (`incdeg=30`, `decdeg=210`)
- deeper thickness (`depth_top=0.5 km`, `depth_bot=8.0 km`)

## Weak-magnetization examples

- `gravmag_sphere_1body_mag_fixedlim_inc90_dec0_weak.in`
- `gravmag_sphere_3body_mag_polygon_incmix_decmix_weak.in`
- `gravmag_sphere_5body_mag_polygon_incmix_decmix_weak.in`
- `gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_weak.in`

These examples use the same geometry as their corresponding base cases,
but with reduced magnetization amplitudes to test low-amplitude response.

- Per-source amplitudes are specified in Card 5; weak cases use `M=0.05 A/m`

## Gravity examples

- `gravmag_sphere_1body_grav_fixedlim_incna_decna.in`
- `gravmag_sphere_1body_grav_polygon_incna_decna.in`
- `gravmag_sphere_3body_grav_polygon_incna_decna.in`

These examples use:

- gravity mode (`ifield=1`)
- inclination and declination not applicable ( `incdeg=na`, `decdeg=na`)

## Run examples

From `examples/`, run one case through the point_source solver and spherical-component converter:

```bash
./run_gravmag_end_to_end.sh 1737.4 \
  lunar_examples/gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient.in
```

Run every case and generate component plots:

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python run_all_examples_brtp.py \
  --solver point_source --examples-dir lunar_examples
```

The batch runner builds the tools, writes `*_xyz.txt` and `*_brtp.txt` under `lunar_examples/output/`, and saves `*_brtp_2x2.png` under `lunar_examples/figs/`. Choose `--solver spectral` for spectral evaluation. See the [solver guide](../GravMagSphere_README.md) for refinement and spectral controls.
