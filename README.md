# IAP_cas-esm2
This repository contains the model source code and scripts for the IAP CCPP project.
The initial code commit is based on tarfiles received in early 2020 by Lulin from IAP.

## Resources
- [Description and Climate Simulation Performance of CAS-ESM Version 2](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002210) paper
- [CCPP Documentation](https://ccpp-techdoc.readthedocs.io)
- CCPP Codebases
  - [CCPP Physics](https://github.com/NCAR/ccpp-physics)
  - [CCPP Framework](https://github.com/NCAR/ccpp-framework)


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

# Background Information
<TODO>


## Adding CCPP Physics Schemes
<TODO>
The following is a guide on adding a CCPP physics scheme.
The [CCPP Documentation](https://ccpp-techdoc.readthedocs.io/en/v7.0.0/) has the most indepth information on this process and might be useful.


### Instructions
<TODO>
Setting up a physics suite for use with the CCPP framework involves three steps:
 - preparing data to be made available to physics through the CCPP
 - running the ccpp_prebuild.py script to reconcile SCM-provided variables with physics-required variables
 - preparing a suite definition file.
