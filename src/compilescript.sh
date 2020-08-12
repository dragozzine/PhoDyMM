#!/bin/bash

# compilescript.sh
# 

g++ -w -O3 -o lcout -lm -lgsl -lgslcblas -fpermissive phodymm.cpp


# makefile from Mariah, could incorporate pretty easily
# could also, pretty easily, make phodymm_lcout.cpp and phodymm_demcmc.cpp
# that had the different flags already set so that this is instant.
# could add "all"
#lcout: phodymm.cpp
#       g++ -w -O3 -o lcout.c -lm -lgsl -lgslcblas -lrt -fpermissive
#phodymm.cpp demcmc: phodymm.cpp
#       mpiicc -w -O3 -o demcmc -lm -lgsl -lgslcblas -lrt -fpermissive
#/gpfs/group/rxd44/default/Mariah_Sim/K80/phodymm.cpp



