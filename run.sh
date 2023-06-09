#!/bin/sh
#PBS -l walltime=1:00:00

# Change to working directory
cd ${PBS_O_WORKDIR}

/home/calvarez/ewald_miniapp/openmp/./main
