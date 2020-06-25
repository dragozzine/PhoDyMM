#!/bin/bash

# Usage: ./make_ttv_posterior_cloud.sh NAME.in
#
# If submit is 1, we will run a model on each posterior draw
# Otherwise, if submit is 0, we will only make the draw and produce the directory structure
submit=1
# How many posterior draws to make (positive integer)
ndraws=25
# What portion of your chains is the burnin (integer from 2-Ngen)
# 10 means discard first 1/10th of each chain, 20 discards first 1/20th etc 
burnindivisor=10

if [ $# -ne 1 ]
  then
    echo "Exactly 1 argument required"
    echo "You must pass the name of the .in file you wish to create the TTV cloud plot for to this script"
    exit 1
fi

if [ ! -f lcout ]; then
    echo "Error: you must have a compiled forward model version of Phodymm named 'lcout' in this directory"
    exit 1
fi

if [ ! -f make_ttv_cloud_helper.py ]; then
    echo "Error: you must have a copy of the python script 'make_ttv_cloud_helper.py' from example_planets/analysis_tools in this directory"
    exit 1
fi

######
# Process files to get set up

infile=$1
#outstr=$(grep -oP "(?<=string outstr= ).*" $infile)
outstr=$(grep "string outstr=" $infile | awk '{print $3}')

cloudttv_dir=cloudttv_dir

mkdir $cloudttv_dir
cp ttv_*.in $cloudttv_dir 
cp $infile $cloudttv_dir

cd $cloudttv_dir

oldt1=$(grep -oP "(?<=double t1= ).*" $infile)
newt1=$(echo "$oldt1 * 3." | bc -l )

sed -i -e "s/string tfile= /string tfile= \.\.\/\.\.\/\.\.\//g" $infile 
sed -i -e "s/double t1= $oldt1/double t1= $newt1/g" $infile 
sed -i -e "s/int xyzflag= 0/int xyzflag= 4/g" $infile 
perl -pe "s/0/1/g if /^int xyzlist=/" $infile > tmp.txt
mv tmp.txt $infile # sometimes -i doesn't work well with perl



#####
# Generate Posterior Draws
# And optionally run the model on them
 
nchains=$(grep -oP "(?<=long nwalkers= ).*" $infile)
ngen=$(grep -oP "(?<=long nsteps= ).*" $infile)
burnin=$(( $ngen / $burnindivisor ))
nbodies=$(grep -oP "(?<=int nbodies= ).*" $infile)
npl=$(( $nbodies - 1 ))
### do I need another ^ ?
starparlinenum=$(( 92 + $npl ))
starparline=$(sed "${starparlinenum}q;d" $infile)
starpar=$(echo $starparline | sed 's/[^01]//g' | awk '{ print length }')

outputdir=posterior_draws

nskip=$(grep -oP "(?<=int outputinterval= ).*" $infile)

outfile="demcmc_${outstr}.out"

if [ ! -f ../$outfile ]; then
    echo "Error: both the NAME.in file and the demcmc_NAME.out file must be present in this directory"
    cd ../
    exit 1
fi

ngen=$(( $ngen / $nskip ))
burnin=$(( $burnin / $nskip ))

baselines=$(( $starpar + 1 ))

mkdir $outputdir

for i in `seq 1 $ndraws`;
do

  random32=$(( ( ($RANDOM & 3)<<30 | $RANDOM<<15 | $RANDOM ) ))
  random=$(($random32 % ($ngen * $nchains) ))
  while [ $random -lt $(($burnin * $nchains )) ]
  do
    random32=$(( ( ($RANDOM & 3)<<30 | $RANDOM<<15 | $RANDOM ) ))
    random=$(($random32 % ($ngen * $nchains) ))
  done

  beginline=$(( $random * ($baselines + $npl) + 1 ))
  endline=$(($beginline + $baselines + $npl - 1 ))
  qline=$(($endline + 1))

  drawfile=${outfile}.draw_${i}
  ddrawfile=${outputdir}/$drawfile

  sed -n "${beginline},${endline}p;${qline}q" ../$outfile > $ddrawfile

  echo 'planet  period (d)       T0 (d)       sqrtecos(omega)     sqrtesin(omega)      i (deg)      Omega (deg)     mp (mjup)      rpors' | cat - $ddrawfile > temp
  mv temp $ddrawfile

  infilei=${infile}.draw_${i}
  dinfilei=${outputdir}/$infilei

  sed -e "9s/\$/\.draw_${i}/" $infile > $dinfilei

  diri=$outputdir/posteriordir_draw_${i}

  echo $diri
  pwd

  mkdir $diri
  mv ${outputdir}/*.draw_${i} $diri
  for pl in `seq 1 9`; do 
      ttvpl=ttv_0$pl.in
      if [ -f $ttvpl ]; then
          cp $ttvpl $diri 
      fi
  done

  if [ $submit -eq 1 ]; then
    cd $diri
      echo "../../../lcout $infilei $drawfile"
      #../../stability $infilei $drawfile 
      ../../../lcout $infilei $drawfile 
   cd ../..
  fi

done


if [ $submit -eq 0 ]; then
    echo "Succesfully drew $ndraws samples from the posterior"
    echo "To run them, switch submit to 1 instead of 0" 
    exit 1
fi

# In principle some planetary systems may go unstable due to the longer baseline
# We assume that doesn't happen here
stablechains=`seq 1 $ndraws` 


#####
# Collate all of the ttv files into one 

for i in $stablechains; 
do
  for j in `seq 1 $npl`;
  do  
    cat $outputdir/posteriordir_draw_$i/tbv00_0${j}.out >> $outputdir/tbv00_0${j}.out
  done
done


####
# Make some plots
cd $outputdir
python ../../make_ttv_cloud_helper.py 
cd ..


####
# Go home
cd ..

 
