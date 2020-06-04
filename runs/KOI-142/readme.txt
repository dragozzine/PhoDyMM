Models of KOI-142 using various different data inputs and model sophistication for demonstration purposes (may not exactly reproduce published results).
    For a simpler introduction to the code, it is suggested you start with example_planets/Kepler-15, Kepler-18, or Kepler-36.  
- RV data from Weiss et al. 2019 (https://arxiv.org/pdf/1909.02427 and Barros et al. 2014)
- Photometric data must be uncompressed with $ tar -xvf kplr005446285_1440_2_nodetrend.tar.gz
                                             $ tar -xvf kplr005446285_1440_2.tar.gz

# Single Forward Model

After phodymm.cpp has been compiled to produce an executable called lcout (see directions in README.md), 
copy the lcout executable to this directory.

Then run the `lightcurve_runscript_MODELNAME.sh` script to generate
 (1) A lightcurve model
 (2) Transit times files for each planet
 (3) Ancillary and diagnostic output as described in the Phodymm README.md Appendix

The models are:
- 2pl
    This model includes only the photometry short cadence where available, and long cadence otherwise.
- 2pl_ttv
    This demonstrates how to fit measured TTVs directly instead of the lightcurve (or in addition to if desired)
    Fitting in this way is not optimized, so it is suggested to fit TTVs one instead use TTVFast or another code
    However, if you have external timing measurements (e.g., from Spitzer) to go with Kepler/K2 data, this example may be useful for a demonstration 
- 2pl_gp
    Here we don't detrend the stellar noise and instead do a simultaneous N-body planet moedl and Gaussian Process fit using Celerite.
    This is generally slower than a model using the detrended data, and may not improve posterior precision/accuracy appreciably
- 3pl_rvs
    This model simultaneously fits the Kepler photometry along with RV data from Keck and Sophie, including a jitter term and constant offset for each data set.
    The third planet is included in this model only, since it is detected solely via radial velocities. 
    This is the most complete model, likely to be used for publication quality models (by exploring posteriors with DEMCMC)

Note: all models are just first guesses, and are not best-fit or optimal models.  

The scripts in example_planets/analysis_tools can then be copied to this directory to produce some quick diagnostic plots 
However, since there are non-transiting planets in this example, before running these scripts you should $ rm tbv00_0[2-3].out
        otherwise they might crash 


# DEMCMC Model Fit

To run a DEMCMC, compile the demcmc executable (directions in README.md), and copy it to this directory.
Depending on your computing setup, you might then run $ ./demcmc_runscript.sh (if you are already on the machine you wish to run the MCMC on) or submit demcmc.sbatch (in order to submit your job to a slurm queue). If you run demcmc_runscript.sh, you will likely want to use screen or nohup, as an MCMC process will often take a very long time.  
To select which model to run in the demcmc, uncomment the appropriate model in demcmc_runscript.sh 

This will generate the file with the state of the MCMC chain every 100 generations which can be used for generating posteriors, as well as several other output files as described in README.md MCMC Fit Output Files section. If you run demcmc_runscript.sh a sub-directory called outf will be created, with diagnostic output from each thread's initialization stored in it, useful for de-bugging. The file demcmc.stdout will track the progress of the MCMC.  

To restart a DEMCMC that has been stopped or crashed, the helper script `restart.sh` is included in the `example_planets/restart_script` directory. After a DEMCMC is run for at least 100 generations (generating all necessary output files) copy it here and run:
$ ./restart.sh demcmc_runscript.sh koi142_MODELNAME.in 
or
$ ./restart.sh demcmc_runscript.sh koi142_MODELNAME.in koi142_rvs.txt
will generate several restart files ending in `.res` (see the Optional Input section in README.md) and a script called demcmc_runscript.sh.res.sh, which can be run to restart the MCMC from where it left off 


 
