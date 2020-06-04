#!/bin/bash 
 
mpirun -np 20 -output-filename outf/sbatch.o --tag-output ./demcmc k2-19_massprior_flat.in k2-19_initialguess.pldin 
