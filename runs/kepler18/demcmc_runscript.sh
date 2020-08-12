#!/bin/bash 
 
mpirun -np 20 -output-filename outf/sbatch.o --tag-output ./demcmc kepler18.in k18.pldin -rv0=kepler18_rvs.txt
