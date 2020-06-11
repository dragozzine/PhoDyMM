#!/bin/bash

#SBATCH --job-name=Kepler36
#SBATCH --output=sbatch.o
#SBATCH --error=sbatch.o
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2

mpirun demcmc kepler36_longcadence.in kepler36.pldin 


