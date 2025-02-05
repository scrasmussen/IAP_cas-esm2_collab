#!/bin/csh
#PBS -N test
#PBS -A project_code
#PBS -l walltime=00:20:00
#PBS -q main
#PBS -j oe
#PBS -l select=2:ncpus=128:mpiprocs=128
#


np=256
mpiexec -np ${np} ./ccsm.exe
