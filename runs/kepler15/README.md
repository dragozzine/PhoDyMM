PhoDyMM/runs/kepler15/README.md

This example of Kepler-15 is one a 1-planet system using a mix of 
long and short cadence data and RVs from 2 different telescopes (Endl 
et al. 2011).

Before running, you must uncompress the photometric data with 
```tar -xvf kplr011359879_1440_2.txt.tar.gz```

See the /runs/kepler36/README file for more detailed information. 
Here we empasize any differences for this particular system.

Note how RV data are sent to lcout and demcmc in the 
lightcurve_runscript.sh and demcmc_runscript.sh files.

One difference is that lcout produces RV outputs in addition to the 
photometric model. There are no analysis tools to plot RV 
data/models.

 
