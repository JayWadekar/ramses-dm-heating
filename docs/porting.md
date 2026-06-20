# Porting the Heating Patch into Another RAMSES Patch

RAMSES uses one `PATCH=...` directory through the build `VPATH`. If your run
already uses another patch directory, such as `patch/init/dice`, merge the edits
below into that existing patch directory rather than trying to stack two patch
directories.

## Files to Merge

Copy this new file into the target patch directory:

```text
dm_heating_fine.f90
```

Then merge the changes from these RAMSES files:

```text
amr_step.f90
poisson_commons.f90
init_poisson.f90
rho_fine.f90
```

## Edit Summary

### `amr_step.f90`

Add the heating call after the hydro update/upload block and before the normal
RT/cooling block:

```fortran
  ! Dark matter heating
  if(hydro.and.poisson)call dm_heating_fine(ilevel)
```

At the bottom of the file, include the new routine with a preprocessor include:

```fortran
#include "dm_heating_fine.f90"
```

Use `#include`, not Fortran `include`, so compile-time flags such as
`DM_HEATING_GAMMA` are applied inside `dm_heating_fine.f90`.

### `poisson_commons.f90`

Add a dark-matter-only density array next to the existing Poisson source density:

```fortran
  real(dp),allocatable,dimension(:)  ::rho_dm
```

### `init_poisson.f90`

Allocate and initialize `rho_dm` wherever `rho` is allocated and initialized:

```fortran
  allocate(rho_dm(1:ncell))
  rho_dm=0
```

### `rho_fine.f90`

Reset `rho_dm` in active cells, virtual boundaries, and physical boundaries
where `rho` is reset.

After particle deposition computes the particle mass contribution

```fortran
vol2(j)=mmm(j)*vol(j,ind)/vol_loc
```

also add this contribution to `rho_dm` only for dark matter particles:

```fortran
if(ok(j).and.is_DM(fam(j)))then
   rho_dm(indp(j,ind))=rho_dm(indp(j,ind))+vol2(j)
end if
```

If your RAMSES build may use `TSC`, make the analogous edit in the TSC
deposition routine as well. Keep this separate from the ordinary `rho`
deposition: `rho` remains the full Poisson source, while `rho_dm` is only the
source for dark matter heating.

Finally, update `rho_dm` ghost zones after particle deposition:

```fortran
call make_virtual_reverse_dp(rho_dm(1),ilevel)
call make_virtual_fine_dp   (rho_dm(1),ilevel)
```

## DICE Notes

The RAMSES DICE initialization patch assigns Gadget halo particles to `FAM_DM`
and stellar particle types to `FAM_STAR`. With the edits above, halo particles
source the heating term, while stars and gas do not.

If your initial-condition reader uses custom particle families, check that the
intended dark matter component satisfies:

```fortran
is_DM(typep(ipart))
```

## Quick Check

After merging, build from the RAMSES `bin/` directory:

```sh
make NDIM=3 SOLVER=hydro MPI=0 EXEC=ramses_dmheat \
  PATCH=/path/to/your/merged/patch \
  USER_FLAGS="-DDM_HEATING_GAMMA=1d-27"
```

Then run a short test with actual dark matter particles. A useful sanity check
is that the patched run should show a positive pressure difference relative to
stock, while the density field remains nearly unchanged over a short run.
