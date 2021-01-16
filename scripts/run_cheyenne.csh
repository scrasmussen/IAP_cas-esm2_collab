#!/bin/csh
#PBS -N test 
#PBS -A P48503002
#PBS -l walltime=00:20:00
#PBS -q regular    
#PBS -j oe
#PBS -l select=20:ncpus=36:mpiprocs=36+1:ncpus=32:mpiprocs=32
#


module list

mpiexec_mpt -n 752 ./ccsm.exe


