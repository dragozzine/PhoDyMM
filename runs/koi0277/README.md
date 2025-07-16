PhoDyMM/runs/kepler36/README.md

This Kepler-36 example is designed to be one of the first for a new user.
It is also an important test that PhoDyMM in installed and running correctly. 
It runs a short test analysis that would require significant extension before being scientifically valuable.
For background on the Kepler-36 system, you can read Carter et al. 2012. 
We will build a basic model of a 2-planet system using only long-cadence data. 

For information on installation and setup of PhoDyMM, see the main README.md. We assume here that you have already completed those instructions. 

# Single Forward Model with lcout

The `lightcurve_runscript.sh` script will run lcout and some analysis tools to generate
 (1) A lightcurve model
 (2) Transit times files for each planet
 (3) Ancillary and diagnostic output as described in the PhoDyMM README.md Appendix
 (4) Multiple output plots (as *.png files) from analysis tools (see that README)
Run this with 
```./lightcurve_runscript.sh```
and see more infomration in the header of that file. 

Becuase PhoDyMM's default is to append to output files, it is better to either move or delete previous files. This script will save previous output 
in a new directory (named after the date and time) files unless the ```--deleteold``` flag is used.

This should take a few seconds, generate some diagnostic output, and produce plots practically identical to those in Ragozzine et al. 2020. If it 
doesn't work or if the plots are different, then you may need to run additional tests. After investigation, consider submitting an "Issue" to the 
PhoDyMM github repository.

You can see the effect of changing planetary parameters by editing kepler36.pldin and rerunning lightcurve_runscript.sh


# DEMCMC Model Fitting with demcmc

PhoDyMM has an automated routine to improve planetary parameters by fitting a photodynamical model to the lightcurve using Differential Evolution 
Markov Chain Monte Carlo (DEMCMC). See Ragozzine et al. 2020 for many more important details. For sufficiently long runs with burnin removed, DEMCMC 
will produce a posterior probability distribution of all the (floating) parameters conditions on the lightcurve data (and the assumed model).

This Kepler-36 example includes a small test run of DEMCMC (20 walkers, 1001 steps). It will use MPI to parallelize over 2 processors/cores. demcmc is 
not designed to be run with a single processor. See the demcmc_runscript.sh header for more information.

Depending on your computing setup, you would run ```./demcmc_runscript.sh``` to run on the machine you are currently using. An example of a sbatch 
file (for a slurm queue) is also given in ```demcmc.sbatch``` but this will likely need to be edited for your specific supercomputing architecture.

As with the lightcurve runscript, you can use --deleteold to remove old files; otherwise they are saved in a subdirectory of oldruns. 

The calculation can take a few minutes. You can use ```tail demcmc.stdout``` (from another terminal) to monitor the progress of the code. (For the 
future, when running long calculations, you may want to use tmux or nohup or a Jupyter Terminal.

This will generate diagnostic output, the file with the state of the MCMC chain every 5 generations which can be used for generating posteriors, as 
well as several other output files as described in PhoDyMM's main README.md (MCMC Fit Output Files section). (Note: thinning by 5 is much smaller 
than you would want to use in practice (100+), but allows our analysis code to have enough numbers to run smoothly.) If you run demcmc_runscript.sh a 
sub-directory called outf will be created, with diagnostic output from each thread's initialization stored in it, potentially useful for de-bugging. 
You can look at and query any ouptut file while you wait (demcmc.stdout showing progress, mcmc_bestchisq_*.aei showing the current best fit, 
demcmc_*.out showing the current state of the thinned chains, etc.)

The demcmc_runscript will also automatically run src/tools/demcmc_quick_analyze which reads in the chains and produces a variety of diagnostic output 
about the DEMCMC runs (e.g., trace plots, corner plots, summary statistics, etc.)


# What next?

If the above codes run correctly and produce the expected output, then the key components of PhoDyMM are working on your machine! Congratulations! 
(If not, try debugging or opening an "Issue" on github.com/dragozzine/PhoDyMM.) Probably the next best thing to do is to go to the other example 
systems and see what they look like. Once you are familiar with running PhoDyMM and have carefully read Ragozzine et al. 2020 for various tips and 
ideas, go forward and do some great science!
