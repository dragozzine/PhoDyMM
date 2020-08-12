#!/bin/bash 
 
# generates a single lightcurve run for the Kepler-36 system
# run using ./lightcurve_runscript.sh [--deleteold]
# use the --deleteold flag if you want to delete former outputs of the lcout file
# If you get the error: 
# -bash: ./lightcurve_runscript.sh: Permission denied
# then make this file executable (chmod u+x lightcurve_runscript.sh)

# Note: whether you save old output files or delete them, this choice
# will be applied to ALL output files in the current directory

# If you want to change the input file parameters or the planetary 
# system parameters, change the input files on the line below that
# calls lcout or make new input files and call these files instead.


# Take Care of Old Files

# By design, PhoDyMM appends results to the end of various files
# but the analysis tools may not catch this. So, the default is 
# to save former output files (which can be recreated using the 
# same input files). 
# Using the --deleteold flag will delete these files instead.

if [[ "$1" == "--deleteold" ]]
  then 
  # deleting files from all previous runs
  echo "Deleting all previous lcout output files"
  rm xyz_*.pldin
  rm tbv*.out
  rm lc_*.lcout
  rm aei_out_*.pldin
  rm *.png

fi


if [ -e lc_*.lcout ]
  then
  # saving files from previous runs

  dt=$(date '+%Y%m%d_%H%M%S')
  echo "saving former lcout files to ./oldrunslc/$dt"

  mkdir -p oldrunslc
  mkdir ./oldrunslc/$dt

  mv xyz_*.pldin ./oldrunslc/$dt
  mv tbv*.out ./oldrunslc/$dt
  mv lc_*.lcout ./oldrunslc/$dt
  mv aei_out_*.pldin ./oldrunslc/$dt
  mv *.png ./oldrunslc/$dt

  # save input files too
  mv kepler36_longcadence.in ./oldrunslc/$dt
  mv kepler36.pldin ./oldrunslc/$dt

fi


echo ""


# Make lightcurve using lcout

../../src/lcout kepler15.in kepler15.pldin -rv0=kepler15_rvs.txt


echo ""
echo "lcout complete"
echo "main output is in lc_kepler36_longcadence.lcout and tbv*.out"


# Run analysis tools for lcout output
python ../../src/tools/lcplot.py
python ../../src/tools/phasefold.py
python ../../src/tools/omcd.py
/
