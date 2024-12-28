# IAP_cas-esm2
This repository contains the model source code and scripts for the IAP CCPP project.
The initial code commit is based on tarfiles received in early 2020 by Lulin from IAP.


## Prerequisites
- C and Fortran Compiler
- NetCDF
- MPI

## Build
After cloning the repository retrieve CCPP Framework and CCPP-Physics repos.
```
$ git submodule update --init --recursive
```

### Setup and Build Case
```
Setup Environment
export NETCDF_PATH=$NETCDF

Switch to CCPP
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

# Adding CCPP Physics Schemes
