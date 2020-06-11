#!/bin/bash

# uses Differential Evolution Markov Chain Monte Carol (DEMCMC)
# to attempt to improve the parameters for the Kepler-36 system
# run using 
# ./demcmc_runscript.sh [--deleteold]
# use the --deleteold flag if you want to delete former outputs of the demcmc file
# If you get the error:
# -bash: ./demcmc_runscript.sh: Permission denied
# then make this file executable (chmod u+x lightcurve_runscript.sh)

# Note: whether you save old output files or delete them, this choice
# will be applied to ALL output files in the current directory

# see demcmc.sbatch for an example of how this might be submitted 
# to a SLURM job queue

# If you want to change the input file parameters or the planetary
# system parameters, change the input files on the line below that
# calls mpirun/demcmc or make new input files and call these files instead.


# Take Care of Old Files

# By design, PhoDyMM appends results to the end of various files
# but the analysis tools may not catch this. So, the default is
# to save former output files (which can be recreated using the
# same input files).
# Using the --deleteold flag will delete these files instead.


if [[ "$1" == "--deleteold" ]]
  then
  # deleting files from all previous runs
  echo "Deleting all previous demcmc output files"
  rm mcmc_bestchisq_*.aei
  rm demcmc*.out
  rm gamma*.txt
  rm demcmc.stdout
  rm -rf ./analysis_dir/

fi


if [ -e demcmc.stdout ] 
  then
  # saving files from previous runs

  dt=$(date '+%d%m%Y_%H%M%S')
  echo "saving former lcout files to ./oldrunsde/$dt"

  mkdir -p oldrunsde
  mkdir ./oldrunsde/$dt
  mkdir ./oldrunsde/$dt/analysis_dir/

  mv mcmc_bestchisq_*.aei ./oldrunsde/$dt
  mv demcmc*.out ./oldrunsde/$dt
  mv gamma*.txt ./oldrunsde/$dt
  mv demcmc.stdout ./oldrunsde/$dt

  mv ./analysis_dir/* ./oldrunsde/$dt/analysis_dir   # move all the analysis to the oldrunsde directory too

  # save input files too
  mv kepler36_longcadence.in ./oldrunsde/$dt
  mv kepler36.pldin ./oldrunsde/$dt

fi


echo ""

# Run demcmc code using MPI

# use 2 processors (-np 2); you can change this, but the number of walkers must be a integer multiple of the
# number of processors

mpirun -np 2 -output-filename outf/sbatch.o --tag-output ../../src/demcmc kepler36_longcadence.in kepler36.pldin


echo ""
echo "demcmc complete, starting demcmc_quick_analyze"
echo ""

# Run analysis tools for demcmc output
python ../../src/tools/demcmc_quick_analyze.py kepler36_longcadence.in

echo ""
echo "consider running lightcurve analysis tools on the best fit at the end of mcmc_bestchisq"
