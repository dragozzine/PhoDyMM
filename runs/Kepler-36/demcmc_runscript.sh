#!/bin/bash 
 
mpirun -np 20 -output-filename outf/sbatch.o --tag-output ./demcmc kepler36_longcadence.in k36.pldin 
