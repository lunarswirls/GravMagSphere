# GravMag Sphere coordinates and magnetic physics

## Scope

This guide describes source-angle conventions, magnetic kernels, and the interpretation of orientation diagnostics. The numerical comparisons below are reference results for their stated configurations; they do not establish accuracy for every geometry or solver setting.

## Body-fixed Cartesian frame

Both grid solvers use:

```text
x = r cos(latitude) cos(longitude)
y = r cos(latitude) sin(longitude)
z = r sin(latitude)
```

The +x axis intersects latitude/longitude zero, +y intersects the equator at 90 degrees east, and +z points north. This corresponds to the SEL convention for lunar examples, not a sun-pointing frame.

Card 5 maps amplitude and global angles to uniform Cartesian magnetization:

```text
mx = amplitude cos(inclination) cos(declination)
my = amplitude cos(inclination) sin(declination)
mz = amplitude sin(inclination)
```

Inclination is elevation above the global xy plane toward +z. Declination is azimuth from +x toward +y. These are not local geophysical angles. A local north/east/down direction must be rotated into the body-fixed frame before constructing Card 5 or a Python magnetization vector.

The implementations are in [the direct solver](../fortran/gravmag_sphere_point_source.f90) and [the spectral solver](../fortran/gravmag_sphere_spectral.f90). The [LPMAG reference workflow](../docs/architecture.md#lpmag-reference-workflow) describes the matching observation-frame requirements.

## Source kernels

### Direct magnetic solver

The direct solver constructs equivalent surface charges for a uniformly magnetized volume:

```text
q = (magnetization · outward_normal) d_area
d_field = mu0 / (4 pi) * q * displacement / |displacement|^3
```

Top, bottom, and side elements contribute to the field. Side-wall normals are oriented outward using an interior-point test. Kernel implementation is in [gravmag_sphere_physics.f90](../fortran/gravmag_sphere_physics.f90), with source construction in the direct executable.

### Spectral magnetic solver

The spectral workflow samples volume dipoles with moment `magnetization * d_volume`, then fits spherical-harmonic coefficients:

```text
field = mu0 / (4 pi) *
        (3 (moment · displacement) displacement / |displacement|^5
         - moment / |displacement|^3)
```

Finite source sampling, harmonic truncation, joint-component fitting, regularization, local correction, and hybrid settings all affect its agreement with direct surface-charge calculations. Neither solver is an exact reference solely by virtue of its formulation.

### Spherical components

`Br` points radially outward, `Btheta` points toward increasing colatitude (south), and `Bphi` points east. [The converter](../fortran/gravmag_xyz_to_brtp.f90) rotates XYZ tables into this basis. Vector components are summed before calculating total magnitude.

## Orientation diagnostics

The recorded coordinate-roundtrip check has a maximum error of `1.33e-15` nT. The unit-direction check of the global angle formula has a maximum norm error of `2.22e-16`.

At inclination 90 degrees, declination has no effect in exact arithmetic. For the same test geometry with declinations 0 and 137 degrees:

| Diagnostic | Direct | Spectral |
|---|---:|---:|
| Maximum Btot difference, nT | 7.0e-06 | 2.13e-02 |
| RMSE of Btot difference, nT | 3.34e-07 | 1.75e-02 |
| Maximum difference as percent of peak field | 1.41e-05 | 0.116 |

The spectral comparison selects `lmax=18`, `reg_lambda=0.05`, and `reg_power=6` for both directions. These finite-precision differences test orientation sensitivity, not arbitrary-case accuracy.

The direct fixed-limit center-point checks give field directions approximately `[+1, 0, 0]` for `inc=0, dec=0` and `[0, -1, 0]` for `inc=0, dec=90`. Field direction need not equal magnetization direction; the source geometry and observation position determine the response.

## Interpretation and development direction

The direct and spectral implementations share the global angle mapping. The recorded orientation checks do not indicate an axis swap or sign inversion as the dominant cause of the [cross-solver residuals](external_solver_comparison.md). Source discretization and spectral fitting still require convergence and independent-reference tests.

Do not interpret a global Card-5 angle as a local north/east/down angle. Likewise, a small total-magnitude residual can hide component-direction errors. Check individual components, consistent units and frames, baseline refinement, and the declared fitting controls.

The [architecture roadmap](../docs/architecture.md#path-forward-body-independent-multi-solver-equivalent-sources) places explicit vector bases, frame transformations, and cross-backend validation in a body-independent interface. Useful orientation regression cases include vertical-inclination declination invariance, known-axis magnetizations, and Cartesian/spherical roundtrips.
