# IAP_cas-esm2
This repository contains the model source code and scripts for the IAP CCPP project.
The initial code commit is based on tarfiles received in early 2020 by Lulin from Institute of Atmospheric Physics (IAP).

## Resources
- [Description and Climate Simulation Performance of CAS-ESM Version 2](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002210) paper
- [CCPP Documentation](https://ccpp-techdoc.readthedocs.io)
- CCPP Codebases
  - [CCPP Physics](https://github.com/NCAR/ccpp-physics)
  - [CCPP Framework](https://github.com/NCAR/ccpp-framework)
- [CAM5 Scientific Documentation](https://ncar.github.io/CAM/doc/build/html/cam5_scientific_guide/index.html)

The Community Earth System Model (CESM) Coupler infrastructure is used and the CESM1.0 documentation can be useful.
- [CESM1.0 Documenation](https://www2.cesm.ucar.edu/models/cesm1.0/cesm/)


# Building and Running the Model
## Prerequisites
- C and Fortran Compiler
- NetCDF
- MPI

## Build
After cloning the repository retrieve CCPP Framework and CCPP-Physics repos.
```
$ git submodule update --init --recursive
```

### Setup Machine
Make sure the machine being used is listed under `scripts/ccsm_utils/Machines/` as `Macros.mach` where `mach` is the name of the machine
If the machine in not listed the user will need to create the setup for a new machine by following the CESM1.0 [Porting to a new machine](https://www2.cesm.ucar.edu/models/cesm1.0/cesm/cesm_doc_1_0_6/c2239.html) documentation.

### Setup and Build Case
```
Source and setup environment variables
$ . scripts/ccsm_utils/Machines/env_machopts.machine_name

If Macros.mach file doesn't handle everything export needed values
$ export NETCDF_PATH=$NETCDF

Switch to CCPP to run CCPP Framework Prebuild steps, required everytime subroutines change
$ ./switch_to_ccpp.sh
$ cd scripts

Choose case from list of avaible cases
$ ./create_newcase -list

The following is an example of the FAMIPC5_FD14 case on Derecho
$ ./create_newcase \
    -case FAMIPC5_FD14 \
    -compset FAMIPC5 \
    -res fd14_fd14 \
    -mach derecho \
    -din_loc_root_csmdata path/to/inputdata


Need to give that configre the correct permissions to run
$ chmod +x /glade/work/soren/src/ccpp/iap/iap_cas-esm2/models/atm/cam/bld/configure
$ cd FAMIPC5_FD14
$ ./configure -case

Remove "bnd_topo2 =..." line from Buildconf/cam.input_data_list
NOTE: On Derecho the NetCDF C and Fortran builds are in different locations.
  1. File models/utils/pio/configure was changed to handle environment
     variables netcdf_c_ROOT and netcdf_fortran_ROOT
  2. File scripts/ccsm_utils/Machines/Macros.derecho was changed to handle
     NetCDF environment variables when C and Fortran installs are in different
     locations

$ ./FAMIPC5_FD14.derecho.build
```

### Run Case
```
Change to scratch directory then enter testcase run directory

$ cd FAMIPC5_FD14/run

Create directory for timing and checkpoint information to be written to
$ mkdir -p timing/checkpoints

Run testcase
$ mpiexec -np 128 ./ccsm.exe &> log.txt

```



# Directory Structure
Note, the whole directory structure is not shown, just a few directories important for this project.
```
IAP_cas-esm2/
├──models/
│   └──atm/iap
│       ├──ccpp/
│       │   ├──config/
│       │   ├──framework/
│       │   ├──physics/
│       │   └──suites/
│       └──src/
│           ├──physics/
│           └──dynamics/
└──scripts/
    ├──SourceMods/
    └──ccsm_utils/
```

<!-- # Background Information -->
<!-- <TODO> -->


<!-- ## Adding CCPP Physics Schemes -->
<!-- <TODO> -->
<!-- The following is a guide on adding a CCPP physics scheme. -->
<!-- The [CCPP Documentation](https://ccpp-techdoc.readthedocs.io/en/v7.0.0/) has the most indepth information on this process and might be useful. -->


<!-- ### Instructions -->
<!-- <TODO> -->
<!-- Setting up a physics suite for use with the CCPP framework involves three steps: -->
<!--  - preparing data to be made available to physics through the CCPP -->
<!--  - running the ccpp_prebuild.py script to reconcile SCM-provided variables with physics-required variables -->
<!--  - preparing a suite definition file. -->

# Using CCPP Physics in IAP


## Instructions


## Scheme Tables

| IAP Schemes Added to CCPP Physics |
|-----------------------------------|
| Zhang-McFarlane convection scheme |



| CCPP Physics Scheme Name              |
|---------------------------------------|
| IAP's Zhang-McFarlane deep convection |
|                                       |

### ZM CCPP Complient Notes
The following is a rough guide for how the IAP's ZM convection scheme was
modified to be make its [physics parameterizations CCPP compliant](https://ccpp-techdoc.readthedocs.io/en/latest/CompliantPhysicsParams.html).

The IAP's version of `zm_conv.F90`, found in `models/atm/iap/src/physics/`, and is
a near identical version of CAM's `zm_conv.F90`, found in `models/atm/cam/src/physics/cam/zm_conv.F90`.
A copy of the ZM convection scheme was made CCPP complient and named
`zm_convr.F90` and can be found in a fork of the CCPP Physics repository, in
the branch [feature/iap_dom](https://github.com/NCAR/ccpp-physics/compare/main...climbfuji:ccpp-physics:feature/iap_dom).
This CCPP complient ZM version is then brought back into the IAP physics as an
example of how to integrate CCPP physics scheme into IAP's.
The CCPP complient ZM version is tested in `ccpp_IAP_test_test1_cap.F90`

`#ifdef CCPP` macros have been added to the following files in `models/atm/iap/src/physics`
| files                   | changes                                                    |
|-------------------------|------------------------------------------------------------|
| `./initindx.F90`        | `use physics_register,   only: convect_deep_register`      |
| `./buffer.F90`          | `real(r8), allocatable, target, dimension(:,:) :: &`       |
|                         | `pblht, tpert, tpert2, qpert2`                             |
|                         | `real(r8), allocatable, target, dimension(:,:,:) :: &`     |
|                         | `qpert`                                                    |
| `./cam_diagnostics.F90` | `use phys_control, only: cam_physpkg_is`                   |
|                         | `public: diag_tphysbc`                                     |
|                         | `character(len=16) :: deep_scheme`                         |
|                         | ZM case to add variables using `cam_history`'s `addfld`    |
|                         | subroutine `diag_tphysbc` for output                       |
| `./tphysbc.F90`         | `cam_out, cam_in` arguments to `tphysbc` subroutine        |
|                         | `use ccpp_data,       only: phys_int_ephem`                |
|                         | `use ccpp_static_api, only: ccpp_physics_run`              |
|                         | `use ccpp_types,      only: ccpp_t`                        |
|                         | call `ccpp_physics_run` and copy data back to IAP vars     |
|                         | - do this for `convect_deep_tend, convect_deep_tend2`      |
| `./runtime_opts.F90`    | `call cldwat_readnl_ccpp(nlfilename)`                      |
|                         | `subroutine zmconv_readnl`                                 |
| `./physpkg.F90`         | this file provides the interface for CAM physics packages  |
|                         | subroutines in the `physpkg` module are                    |
|                         | - `phys_inidat`                                            |
|                         | - `phys_init`                                              |
|                         | - `phys_run1`                                              |
|                         | - `phys_run1_adiabatic_or_ideal`                           |
|                         | - `phys_run2`                                              |
|                         | - `phys_final`                                             |
|                         | - `cdata_init`: added by CCPP                              |
|                         | CCPP adds arguments to `phys_{init, run1, final}` routines |
| .`/physics_types.F90`   |                                                            |
|                         | Adds subroutines                                           |
|                         | - `interstitial_ephemeral_create`                          |
|                         | - `interstitial_ephemeral_reset`                           |
|                         | - `interstitial_persistent_associate`                      |
|                         | - `interstitial_persistent_create`                         |
|                         | - `interstitial_persistent_init`                           |
|                         | - `physics_global_init`                                    |


- [PR#1](https://github.com/lulinxue/IAP_cas-esm2/pull/1)
  - [ ] track down email referenced in the PR
<!-- commit messages   -->
<!-- - add CCPP framework and physics as submodules (master branches for now) -->
<!-- - add initial ccpp_prebuild_config.py and SDF (with xsd) -->
<!-- - save work for ZM convection (incomplete) | commiting work-in-progress code as of 20210713 -->
<!-- - update comments for where everything went in zm_conv_intr.F90 for others to follow -->
<!-- - first version of host files that pass ccpp_prebuild.py -->
<!-- - add CCPP API calls; reconfigure ccpp_prebuild_config.py to write compilable caps -->
<!-- - update ccpp/physics submodule pointer, fix error in physpkg.F90, add -DCCPP to Cheyenne macros file -->
<!-- - fix compilation errors in buffer.F90 and physics_types.F90 -->
<!-- - change kind_r8 to shr_kind_r8 in metadata | fix compilation errors -->
<!-- - move data accessed by CCPP from cam_comp.F90 to ccpp_data.F90 -->
<!-- - changes to enable first successful compilation -->
<!-- - Update .gitmodules | Update machine and case config to work for Dom -->
<!-- - Bugfixes for phys_global and cdata_domain -->
<!-- - Remove unused models/atm/cam/src/control/cam_comp.meta -->
<!-- - Make ccpp_data local variables in -->
<!--   models/atm/iap/src/physics/physpkg.F90. Correct metadata for several DDTs in -->
<!--   models/atm/cam/src/dynamics/iap/ccpp_data.F90 -->
<!-- - Add metadata for cam_out_t in models/atm/cam/src/control/camsrfexch_types.* -->
<!-- - Update ccpp-framework and ccpp-physics from NCAR main, update host model metadata -->
<!-- - Major bug fix in cam/src/cpl_mct/atm_comp_mct.F90: import 'cam_in' from -->
<!--   ccpp_data instead of defining it locally  | Update submodule pointer for ccpp-framework -->
<!-- - Remove unused nchnks -->
<!-- - Add metadata for begchunk and endchunk in models/atm/iap/src/physics/ppgrid.meta -->
<!-- - Bug fixes in models/atm/iap/src/physics/buffer.meta: use correct CCPP block numbers and kinds -->
<!-- - Add buffer.F90/meta to models/atm/iap/ccpp/config/ccpp_prebuild_config.py -->
<!-- - Remove unused/wrong metadata table hooks from models/atm/cam/src/cpl_esmf/atm_comp_esmf.F90 -->
<!-- - Temporary or permanent? Use assumed-shape array declarations in models/atm/iap/src/physics/convect_shallow.F90 -->
<!-- - Clean up in models/atm/iap/src/physics/physics_types.* -->
<!-- - Add temporary debugging information in models/atm/iap/src/physics/tphysbc.F90 -->
<!-- - Bug fix in models/atm/iap/src/physics/physpkg.F90: correct thread number in Fortran is 1..N -->
  - Add CCPP prebuild config for IAP AGCM
  - Add Suite Definition File (SDF): `suite_IAP_test.xml`
  - Make [CCPP compliant](https://ccpp-techdoc.readthedocs.io/en/latest/CompliantPhysicsParams.html)
    - Add metadata by adding `.meta` files for `.F90` files that will be
      access by framework
      - `IAP_physics_ephemeral_interstitial_instance_all_blocks`
      - `IAP_physics_persistent_interstitial_instance_all_blocks`
      - `IAP_physics_global_instance`
  - Add `#ifdef CCPP` macros

- [PR#2](https://github.com/lulinxue/IAP_cas-esm2/pull/2) Feature/iap dom 2
I had to CCPP-ize `convect_deep_tend_2` (there was no way to do the first part
in CCPP and the second outside of CCPP), create a separate group in the XML
SDF. I also had to fix several bugs along the way, and update ccpp-framework
and ccpp-physics with the latest code from NCAR main. This made it possible to
get rid of the `dummy_loop` scheme.

After returning from tphysbc - and all state/pbuf variables being identical
for all MPI tasks and all chunks - the code crashes in the newly created CCPP
diag output call. I tested this, if I simply bypass the diag output then it
crashes later on in a nother diag routine from other parts of the model
(somewhere in the dynamics). So what I did is to stop the model after the
first pass through the physics.

The ccpp-framework changes are probably not of interest, but in the
ccpp-physics changes I would look for the changes in the `zm_*.{F90,meta}`
files.

<!-- commit messages   -->
<!-- -- Wrap contents of module ccpp_data in #ifdef CCPP -->
<!-- -- Remove dummy_loop from models/atm/iap/ccpp/config/ccpp_prebuild_config.py, add STATIC_API_SOURCEFILE -->
<!-- -- Add group test2 to models/atm/iap/ccpp/suites/suite_IAP_test.xml -->
<!-- -- Remove temporary testing code from models/atm/iap/src/physics/convect_shallow.F90 -->
<!-- -- Update submodule pointers for ccpp-framework and ccpp-physics -->
<!-- -- Read cldwat namelists for both CCPP and non-CCPP code -->
<!-- -- Add DDT print routines to models/atm/iap/src/physics/physics_types.F90,  -->
<!-- add new DDT members to models/atm/iap/src/physics/physics_types.{F90,meta} -->
<!-- -- Initialize cldwat for both non-CCPP and CCPP physics in models/atm/iap/src/physics/physpkg.F90 -->
<!-- -- Change default debug levels in scripts/ccsm_utils/Case.template/config_definition.xml -->
<!-- -- Change FFLAGS_NOOPT in scripts/ccsm_utils/Machines/Macros.cheyenne -->
<!-- -- Convenience scripts for quickly switching between CCPP and non-CCPP code -->
<!-- -- Temporary change in models/atm/cam/src/physics/cam/phys_buffer.F90 - -->
<!--    initialize phys_buffer fld_ptr to 0.0 instead of NaN so that it can be printed -->
<!--    using tphysbc.F90 -\-> pbuf_print_data -->
<!-- -- Initialize cldwat for CCPP and no CCPP, add temporary stop of code after calling tphysbc for all chunks -->
<!-- -- models/atm/iap/src/physics/tphysbc.F90: CCPP-ize convect_deep_tend_2, -->
<!--    transfer CCPP interstitial data to local arrays, add debugging prints, more -->
<!--    bugfixes -->
<!-- -- One more print statement in models/atm/iap/src/physics/physpkg.F90 -->
<!-- -- Final touches -->
