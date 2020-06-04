#!/bin/bash

## usage:$ ./restart.sh demcmc.sbatch kid009632895.in [rv1.txt] [rv2.txt]
##       $ sbatch demcmc.sbatch.res

if [ $# -lt 2 ]
  then
    echo "At least 2 arguments required"
    echo "You must pass (1) the name of the initial runscript and" 
    echo "              (2) the name of the .in file you originally used"
    echo "Also pass the names of any rv data files used in your fit"
    exit 1
fi

SBATCHF=$1
INFILE=$2
RVFILE0=$3
RVFILE1=$4

if [ -z $SBATCHF ] || [ ! -f $SBATCHF ]; then
    echo "Error: The RUNSCRIPTNAME.sh file must be present in this directory"
    exit 1
fi
if [ -z $INFILE ] || [ ! -f $INFILE ]; then
    echo "Error: The NAME.in file must be present in this directory"
    exit 1
fi

OUTSTR=$(grep "string outstr=" $INFILE | awk '{print $3}')
OUTFILE="demcmc_$OUTSTR.out"
if [ ! -f $OUTFILE ]; then
    echo "Could not find $OUTFILE"
    echo "Error: Both the NAME.in file and the demcmc_NAME.out file must be present in this directory"
    exit 1
fi
NBODIES=$(grep "int nbodies=" $INFILE | awk '{print $3}')
NPL=$(($NBODIES-1))
BSQSTR="mcmc_bestchisq_$OUTSTR.aei"
if [ ! -f $BSQSTR ]; then
    echo "Could not find $BSQSTR"
    echo "Error: Both the NAME.in file and the bestchisq_NAME.aei file must be present in this directory"
    exit 1
fi
REGEN=$(tac $OUTFILE | grep ', 0,' -m 1 | tac | awk '{print $6}')
GAMMAF="gamma_$OUTSTR.txt"
if [ ! -f $GAMMAF ]; then
    echo "Could not find $GAMMAF"
    echo "Error: Both the NAME.in file and the gamma_NAME.txt file must be present in this directory"
    exit 1
fi

starparlinenum=$(( 92 + $NPL ))
starparline=$(sed "${starparlinenum}q;d" $INFILE)
NSTARPAR=$(echo $starparline | sed 's/[^01]//g' | awk '{ print length }')

RESTARTF="$OUTFILE.res"
RESTARTBSQ="$BSQSTR.res"
RESTARTGAMMA="$GAMMAF.res"
RESTARTSBATCH="$SBATCHF.res.sh"

tac $OUTFILE | grep ', 0,' -m 1 -B 9999 -C $(($NPL+$NSTARPAR)) | tac > $RESTARTF
tac mcmc_bestchisq_$OUTSTR.aei | grep 'planet' -m 1 -B 9999 | tac > $RESTARTBSQ
tac $GAMMAF | grep $REGEN -m 1 | tac > $RESTARTGAMMA

if [[ ! -s $RESTARTGAMMA ]] ; then 
    echo "ERROR!!"
    echo "Something went wrong and restarting will fail"
fi ;


####cat $SBATCHF | grep '#SBATCH --ntasks=' -B 9999 > $RESTARTSBATCH
echo "#!/bin/bash " > $RESTARTSBATCH
echo " " >> $RESTARTSBATCH
touch empty.txt

if [ -z "$RVFILE1" ]
then
  if [ -z "$RVFILE0" ] 
  then
    ##echo mpirun ./demcmc $INFILE empty.txt $RESTARTF $RESTARTBSQ $RESTARTGAMMA >> $RESTARTSBATCH
    echo mpirun -np 20 -output-filename outf/sbatch.o --tag-output demcmc $INFILE empty.txt $RESTARTF $RESTARTBSQ $RESTARTGAMMA >> $RESTARTSBATCH 
  else 
    echo mpirun -np 20 -output-filename outf/sbatch.o --tag-output demcmc $INFILE empty.txt -rv0=$RVFILE0 $RESTARTF $RESTARTBSQ $RESTARTGAMMA >> $RESTARTSBATCH 
    ##echo mpirun ./demcmc $INFILE empty.txt -rv0=$RVFILE0 $RESTARTF $RESTARTBSQ $RESTARTGAMMA >> $RESTARTSBATCH
  fi
else
  ##echo mpirun ./demcmc $INFILE empty.txt -rv0=$RVFILE0 -rv1=$RVFILE1 $RESTARTF $RESTARTBSQ $RESTARTGAMMA >> $RESTARTSBATCH
  echo mpirun -np 20 -output-filename outf/sbatch.o --tag-output demcmc $INFILE empty.txt -rv0=$RVFILE0 -rv1=$RVFILE1 $RESTARTF $RESTARTBSQ $RESTARTGAMMA >> $RESTARTSBATCH 
fi

chmod +x $RESTARTSBATCH

