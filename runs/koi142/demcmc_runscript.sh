#!/bin/bash 
 
#mpirun -np 20 -output-filename outf/sbatch.o --tag-output ./demcmc koi142_2pl_gp.in koi142_2pl_gp.pldin 
#mpirun -np 20 -output-filename outf/sbatch.o --tag-output ./demcmc koi142_2pl.in koi142_2pl.pldin  
#mpirun -np 20 -output-filename outf/sbatch.o --tag-output ./demcmc koi142_2pl_ttvs.in koi142_2pl.pldin
mpirun -np 20 -output-filename outf/sbatch.o --tag-output ./demcmc koi142_3pl_rvs.in koi142_3pl_rvs.pldin -rv0=koi142_rvs.txt
