# RAMSES Dark Matter Heating Patch

This repository contains a minimal RAMSES patch for injecting local heat into gas
from dark matter decay-like energy deposition. It is intended to accompany the
paper draft:

> Using the circumgalactic medium of dwarf galaxies as a calorimetric dark matter detector

The implemented source term follows

```text
dE / (dt dV) = rho_DM * Gamma_heat * c^2
```

where `rho_DM` is built from RAMSES particles with family `FAM_DM` only. Gas,
stars, sinks/cloud particles, debris, and analytic gravity density sources are
not counted as dark matter for the heating term.

## Contents

- `patch/`: RAMSES source overrides. Use this directory as the RAMSES `PATCH`.
- `tests/`: a small particle-based smoke test that verifies the heating path is
  active for actual dark matter particles.

## Compatibility

This patch was developed against RAMSES 3.0 on the `dev` branch, commit
`7050a55b57a71ab85d7cae0710254c80aaaa403e`.

The patch overrides RAMSES internals including `amr_step.f90`, `rho_fine.f90`,
`poisson_commons.f90`, and `init_poisson.f90`. If you use a substantially
different RAMSES version, review these files against your local source tree.

### Custom RAMSES Patches and DICE Initial Conditions

The heating source is intended to be general for RAMSES simulations with
particle dark matter. It uses RAMSES particle families and deposits heat from
particles satisfying `is_DM(typep)`.

This is compatible with the RAMSES DICE initialization patch: in the DICE
Gadget reader, halo particles are assigned `FAM_DM`, while stellar particle
types are assigned `FAM_STAR`. Those stellar particles are therefore not counted
as dark matter heating sources.

One practical caveat is that RAMSES `PATCH=...` selects one patch directory via
the build `VPATH`; it does not compose multiple patch directories automatically.
If your simulation already uses another RAMSES patch, such as `patch/init/dice`,
merge this repository's changes into that patch directory. The files that
usually need attention are:

- `amr_step.f90`: add the `dm_heating_fine(ilevel)` call and include
  `dm_heating_fine.f90`.
- `poisson_commons.f90`: add the `rho_dm` array.
- `init_poisson.f90`: allocate and initialize `rho_dm`.
- `rho_fine.f90`: deposit only `FAM_DM` particles into `rho_dm`.
- `dm_heating_fine.f90`: add this new source file to the patch directory.

## Build

From your RAMSES `bin/` directory:

```sh
make NDIM=3 SOLVER=hydro MPI=0 EXEC=ramses_dmheat PATCH=/path/to/ramses-dm-heating/patch
```

By default the patch uses

```text
Gamma_heat = 1e-27 s^-1
```

To choose another value, pass a preprocessor definition through `USER_FLAGS`:

```sh
make NDIM=3 SOLVER=hydro MPI=0 EXEC=ramses_dmheat \
  PATCH=/path/to/ramses-dm-heating/patch \
  USER_FLAGS="-DDM_HEATING_GAMMA=1d-26"
```

The paper draft considers runs with `Gamma_heat = 1e-27, 1e-26, 1e-25 s^-1`
and a control run without this patch.

## Smoke Test

The included test uses RAMSES' ASCII particle loader. It places eight static
dark matter particles in a low-density gas box and compares a stock run to a
patched run.

From `tests/`:

```sh
mkdir -p run_stock run_dmheat

cd run_stock
/path/to/ramses/bin/ramses_stock3d ../dm_heating_particles.nml > stock.log

cd ../run_dmheat
/path/to/ramses/bin/ramses_dmheat3d ../dm_heating_particles.nml > dmheat.log
```

The patched run should show a positive pressure difference relative to the stock
run while density remains nearly unchanged over the short test.

## Implementation Notes

The patch stores a separate `rho_dm` array in `poisson_commons`. During particle
deposition in `rho_fine.f90`, RAMSES still builds the ordinary Poisson source
`rho` from all relevant gravitating components, but `rho_dm` receives only
particles satisfying `is_DM(typep)`.

`dm_heating_fine.f90` then updates the gas total energy variable:

```fortran
uold(ind_leaf(i),ndim+2) = uold(ind_leaf(i),ndim+2) + rho_dm * c^2 * Gamma_heat * dt
```

The source is operator-split and called after the hydro update and before the
cooling/chemistry block in `amr_step.f90`.

## Limitations

- The patch models local energy deposition only.
- The heating rate is compile-time configurable, not a runtime namelist
  parameter.
- The current implementation is for RAMSES hydro builds with Poisson enabled.
- Analytic gravity density profiles are deliberately not counted as dark matter
  heating sources.

## Citation

If you use this patch, please cite Mintz et al. 2026 (in prep.) and link to
this repository for reproducibility.
