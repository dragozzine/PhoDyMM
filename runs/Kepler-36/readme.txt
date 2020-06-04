A basic model of a 2-planet system using only long-cadence data.


# Single Forward Model

After phodymm.cpp has been compiled to produce an executable called lcout (see directions in README.md), 
copy the lcout executable to this directory.

Then run the `lightcurve_runscript.sh` script to generate
 (1) A lightcurve model
 (2) Transit times files for each planet
 (3) Ancillary and diagnostic output as described in the Phodymm README.md Appendix


# DEMCMC Model Fit

To run a DEMCMC, compile the demcmc executable (directions in README.md), and copy it to this directory.
Depending on your computing setup, you might then run $ ./demcmc_runscript.sh (if you are already on the machine you wish to run the MCMC on) or submit demcmc.sbatch (in order to submit your job to a slurm queue). If you run demcmc_runscript.sh, you will likely want to use screen or nohup, as an MCMC process will often take a very long time.  

This will generate the file with the state of the MCMC chain every 100 generations which can be used for generating posteriors, as well as several other output files as described in README.md MCMC Fit Output Files section. If you run demcmc_runscript.sh a sub-directory called outf will be created, with diagnostic output from each thread's initialization stored in it, useful for de-bugging. The file demcmc.stdout will track the progress of the MCMC.  

To restart a DEMCMC that has been stopped or crashed, the helper script `restart.sh` is included in the `example_planets/restart_script` directory. After a DEMCMC is run for at least 100 generations (generating all necessary output files) copy it here and run:
$ ./restart.sh demcmc_runscript.sh kepler36_longcadence.in 
will generate several restart files ending in `.res` (see the Optional Input section in README.md) and a script called demcmc_runscript.sh.res.sh, which can be run to restart the MCMC from where it left off 

