#!/bin/bash
# compilescript.sh
# A script that compiles all the PhoDyMM executables
# Because of the way that PhoDyMM is built, this script is executed in lieu of a Makefile

# Use:
# Edit the System Specific Setup Variables to match your system:
# - CC is the C++ compiler (probably no need to edit)
# - MPIC is the Message Passing Interface C++ compiler (probably no need to edit)
# - CFLAGS are the flags to use (no need to edit)
# - CELERITEINCLUDE is the directory where celerite is installed on your machine (see README for more info)
# - LDFLAGS are the required libaries (no need to edit, though you may need to install GSL)
#
# Then run ./compilescript.sh:
# chmod u+x ./compilescript.sh
# ./compilescript.sh


# System Specific Setup variables
CC="g++"
MPIC="mpic++"
CFLAGS="-w -O3"
MPIFLAGS="-w -Ofast"
CELERITEINCLUDE="-I/home/byu.local/dragozzi/research/celerite/cpp"
LDFLAGS="-lm -lgsl -lgslcblas"

echo "Compiling all the PhoDyMM codes (~1 minute)"


echo "compiling lcout without celerite"
# set #define demcmc_compile 0
sed -i 's/#define demcmc_compile 1/#define demcmc_compile 0/g' phodymm.cpp 
# set #define celerite_compile 0
sed -i 's/#define celerite_compile 1/#define celerite_compile 0/g' phodymm.cpp 

$CC $CFLAGS -o lcout -fpermissive phodymm.cpp $LDFLAGS
chmod u+x lcout


echo "compiling lcout with celerite"
# set #define celerite_compile 0
sed -i 's/#define celerite_compile 0/#define celerite_compile 1/g' phodymm.cpp 
$CC $CFLAGS -o lcout_cele $CELERITEINCLUDE/include $CELERITEINCLUDE/lib/eigen_3.3.3 -fpermissive phodymm.cpp $LDFLAGS
chmod u+x lcout_cele



echo "compiling demcmc without celerite"
# set #define demcmc_compile 1
sed -i 's/#define demcmc_compile 0/#define demcmc_compile 1/g' phodymm.cpp 
# set #define celerite_compile 0
sed -i 's/#define celerite_compile 1/#define celerite_compile 0/g' phodymm.cpp 

$MPIC $MPIFLAGS -o demcmc -fpermissive phodymm.cpp $LDFLAGS
chmod u+x demcmc


echo "compiling demcmc with celerite"
# set #define celerite_compile 0
sed -i 's/#define celerite_compile 0/#define celerite_compile 1/g' phodymm.cpp 
$MPIC $MPIFLAGS -o demcmc_cele $CELERITEINCLUDE/include $CELERITEINCLUDE/lib/eigen_3.3.3 -fpermissive phodymm.cpp $LDFLAGS
chmod u+x demcmc_cele





echo "All done!"
echo "Any compiling errors probably mean you need to edit the compilescript.sh file"
