# PhoDyMM: the PhotoDynamical Multiplanet Model

## Attribution

This code is free to use and modify, but you shoud cite
Ragozzine et al. (2020), PASP, submitted 
which is included in the doc directory if you use PhoDyMM. Reading through that paper is strongly recommended before use.
The vast majority of the original code was developed by Sean Mills with support from Daniel Fabrycky. Darin Ragozzine and his group augmented the code and currently maintain it. 

## Requirements

PhoDyMM is written in C and is designed to be executed on UNIX/Linux machines using a Bash shell. Valuable analysis tools require python. 

This code requires: 

  1. GSL (gnu scientific library - https://www.gnu.org/software/gsl/)

If you are using the demcmc rather than just running the forward model, it also requires:

  2. MPI (https://www.open-mpi.org)

[OPTIONAL] If you wish to use Gaussian processes in the fit, it requires:

  3. celerite (https://github.com/dfm/celerite) 

## Running the model

After those are installed, you may compile PhoDyMM from source file src/phodymm.cpp using:
```
$ g++ -w -O3 -o lcout -I/yourpathto/celerite/cpp/include -I/yourpathto/celerite/cpp/lib/eigen_3.3.3 -lm -lgsl -lgslcblas -fpermissive phodymm.cpp
```
where you would replace "yourpathto" with the path to your celerite install.
Or, if you are not using celerite for GPs, you must change the second line from 
```#define celerite_compile 1```
to 
```#define celerite_compile 0```
and then compile with:
```
$ g++ -w -O3 -o lcout -lm -lgsl -lgslcblas -fpermissive phodymm.cpp
```
The `-lm` `-lgsl` and `-lgslcblas` flags should link the compiler to your math libraries (including gsl). You may find you need to include a `-L/yourpathtogsllibraries/lib` before the `-lgsl` link. If you are compiling on a Mac, I recommend you do not use Apple's default g++, but rather one installed from scratch with homebrew/macports to avoid issues with fpermissive. 

This generates an executable called `lcout` (short for light-curve output).
You may use it to run a forward N-body model given a data file, input file (*.in), and initial conditions (*.pldin) file. These files are described below.
The output will be a theoretical lightcurve (in *.lcout) and list of transit times (in tbv*.out for each planet). These files are described below.
For an example, see the readme.txt in `runs/Kepler-36`

## Fitting to Data to Generate Posteriors

To run a demcmc, you must change the first line in phodymm.c from
```#define demcmc_compile 0```
to 
```#define demcmc_compile 1```

Then recompile with
```
$ mpic++ -w -Ofast -o demcmc -I/yourpathto/celerite/celerite/cpp/include -I/yourpathto/celerite/celerite/cpp/lib/eigen_3.3.3 -lm -lgsl -lgslcblas -lmpi -fpermissive phodymm.cpp
```
or, if you are not using celerite, on line 2 set:
```#define celerite_compile 0```
and then compile with:
```
$ mpic++ -w -Ofast -o demcmc -lm -lgsl -lgslcblas -lmpi -fpermissive phodymm.cpp
```

This will create an executable called `demcmc` which can be used to run a differential evolution MCMC (DEMCMC). 

In general, compiled files and bash scripts may need to have executable permissions set using chmod u+x

It is recommended you run the demcmc executable on a computing cluster. Some example scripts to run the DEMCMC are included in 
`runs/Kepler-36`
Input and output files are discussed below. For more details look at the readme.txt in that folder.  

## Step-by-step Instructions

Here are some step-by-step instructions for novice users that also function as a test case. 

1. Clone this repository to your local machine using
``` 
git clone https://github.com/dragozzine/PhoDyMM.git
```

2. Go into the src directory and compile lcout and demcmc

3. Go into the runs/Kepler-36 directory and make all scripts executable (if not already):
```
chmod u+x *.sh
```

4. Examine the lightcurve model using this script:
```
./lightcurve_runscript.sh
```
This should take about 10 seconds and produce multiple output figures. 

5. Do a short test DEMCMC run using this script:
```
./demcmc_runscript.sh
```
You can observe the progress of this DEMCMC run using ```tail demcmc.stdout```.

## Appendix


### Units

The units used internally are AU and days. Generally all quantities should be given in those units. 


### Input Files

1. PHOTOMETRY DATA FILE  

   PhoDyMM is designed to for use on detrended Kepler lightcurves. See kepler_detrend for code to produce these lightcurves automatically. Examples 
in the runs directories have data files, though these may need to be unpacked (using tar -xvzf filename.gz). Any detrended lightcurve will work -- see 
runs/K2-19 for an example from K2 -- but many of the defaults only make sense for Kepler. 

   The file with the list of input times (e.g., runs/Kepler-36/kplr011401755_1440_1.txt) must have the following format:
   ```
   [Line Index] \t [Time (days)] \t [Ignored] \t [Ignored] \t [Flux] \t [Flux Error]
   ```
   In units:
   ```
   [none] \t [days] \t [none] \t [none] \t [Normalized to 1] \t [Relative to Flux Value]
   ```
   With type:
   ```
   [long int] \t [double] \t [numeric] \t [numeric] \t [double] \t [double]
   ```
   The 1st, 3rd, and 4th columns are currently ignored, but must be present. They represent the data point index from Kepler and the 
raw/un-normalized flux and uncertainty values. If short and long cadence data are used simultanously, the data file must also contain an additional 
column indicating the cadence of each point. The column should be of type int and is simply a 1 if the point is long cadence, and 0 if it is short 
cadence.
   
   If only a forward model is being computed, then only the Time column (and cadence if present) are used.
   

2. INPUT FILE  

   This file (*.in) specifies various parameters necessary for completing the fit (e.g., runs/Kepler-36/kepler36_longcadence.in). Examples 
of paramters to edit in this file include the name of the run, the location of the data file, the fitting basis, which parameters to let vary or keep 
fixed, and priors. This file is read in by the C code and must be in the exact format as the example. Commented lines must be retained.

   The file is structured as the description of each variable in 1 or more commented lines beggining with //, followed by the variable name and 
declartion in the format:
   ```
   type name= value
   ```
   Users should change the value, but leave the first two columns unchanged. It is important that spacing is retained.  

   All entries should be sufficiently described in the input file or Ragozzine et al. 2020. 
 

3. INITIAL CONDITIONS  

   The starting state of the planetary system at the epoch specified in the input file must be provided by the user and is called the initial condition file (e.g., runs/Kepler-36/k36.pldin).  
   
   The format of this file is first a line which is ignored describing the column names. The first column is the planet name, the last the planetary radius in units of stellar radius, and the second to last the planetary mass in Jupiter masses. The other columns depend on the choice of input parameter basis selected in the input file (by the xyzflag variable).   
   ```
   [Planet Label] \t  [column name] \t [column name] \t [column name] \t [column name] \t [column name] \t [column name] \t [Mass] \t [Rp/Rstar] 
   [value] \t  [value] \t [value] \t [value] \t [value] \t [value] \t [value] \t [value] \t [value] 
   [[value] \t  [value] \t [value] \t [value] \t [value] \t [value] \t [value] \t [value] \t [value]] 
   ...
   [[value] \t  [value] \t [value] \t [value] \t [value] \t [value] \t [value] \t [value] \t [value]]
   ```
   This is followed by at least 5 rows of single value entries related to the overal system or stellar properties. The following must be included:
   ```
   [Stellar Mass (Solar)]
   [Stellar Radius (solar)]
   [c_1 quadratic limb darkening term]
   [c_2 quadratic limb darkening term]
   [dilution: the fraction of the flux in the system not from the star which the planets are transiting]
   ```
   Several optional rows follow, depending on what fitting is being done. For instance, celerite GP fits require the 4 celerite terms (one per row). RV fits can include an RV jitter term for each set of observations (1 per line), etc. For more details refer to the runs and input file. Extra lines at the end of this file are ignored.  
   
   Typically, the xyzflag in the input file is set to 0, and the basis for the initial condition is specified as:
   ```
   [Planet Label] \t  [period (d)] \t [T0 (d)] \t [e] \t [i (deg)] \t [Omega (deg)] \t [omega(deg)] \t [Mass] \t [Rp/Rstar]
   ```
   Where T0 is the time of conjunction of the planet and the star, i is the incliation, Omega is the nodal angle, and omega is the argument of periastron with the coordinate system such that the sky plane is 0. This initial condition is transformed into the basis selected in the input file before fitting (e.g. {e, omega} -> {sqrt(e)*sin(omega), sqrt(e)*cos(omega)).   
   
   Planets should be entered in increasing order of orbital period. 

   Other options for xyzflag are 1 (Jacobian Cartesian coordinates), 2 (Stellar-centric Cartesian coordinates), or 3 (barycentric Cartesian Coordinates) with format
   ```
   [Planet Label] \t  [x] \t [y] \t [z] \t [v_x] \t [v_y] \t [v_z] \t [Mass] \t [Rp/Rstar]
   ```
   where the units are AU and AU/day. If any of these options are selected, xyzlist must be set to 1 for each planet. An example of usage with a cartesian coordinate style initial conditions instead of the usual orbital elements is found in `runs/Kepler-18/` indicated with filenames containing `cartesian`. 

   You can also choose xyzflag= 5 or xyzflag= 6 to enter elements in a-e-i basis:
   ```
   [Planet Label] \t  [a (AU)] \t [e] \t [i (deg)] \t [Omega (deg)] \t [omega(deg)] \t [angle (deg)] \t [Mass] \t [Rp/Rstar]
   ```
   where angle is the true anomaly for xyzflag= 5, and the mean anomaly for xyzflag= 6. 
   
   Finally, setting xyzflag= 4, means the input is assumed to be in the DEMCMC parameters format, typically:
   ```
   [Planet Label] \t  [period (d)] \t [T0 (d)] \t [sqrt(e)\*cos(omega)] \t [sqrt(e)\*sin(omega)] \t [i (deg)] \t [Omega (deg)] \t [Mass] \t [Rp/Rstar]
   ```
   although this may vary depending on your setup in the .in file. This is useful for using draws from the posterior (demcmc_NAME.out) as initial conditions. Use in this manner is seen in the example posterior TTV cloud script described below (`make_ttv_posterior_cloud.sh`). 


4. [Optional] Radial Velocity Datasets

   Radial velocity data sets can be passed to either the `lcout` or `demcmc` executables for forward model fitting or RV fitting. Files need to have the format:
   ```
   [time (days)] \t [RV (m/s) \t [Error (m/s)] \t [telescope number]
   ```
   The time series must have the same 0 point as the flux time. The telescope number is an integer index starting from 0 indentify each unique telescope used for observations. This information is used to determine different constant RV offsets for different telescopes, and allows for different RV jitters for different telescopes. 

   See `runs/Kepler-15` or `runs/Kepler-18` for example RV data usage. 


5. [Optional, `demcmc` only] Restart Files

   If a demcmc run is stopped it may be resumed from the last recorded generation by passing `demcmc` additional arguments that specify the location of "restart files." These files include the current state of the MCMC chain, current MCMC scaling factor, and the best-fit solution found so far. These files can be generated from the ends of the MCMC output files #1-3 as described below. An example script to automatically generate them correctly is included in `runs/restart_script/restart.sh`. It should be copied to the directory where the MCMC was run, and then invoked with 
   ```
   ./restart.sh demcmc_runscript.sh NAME.in
   ``` 
   and then the MCMC may be restarted by running the generated `demcmc_runscript.sh.res.sh` file. 
   


### Forward Model (`lcout`) Output Files

Note: Most of the output files append to existing files! instead of overwriting to avoid accidental erasure. To make sure your files don't contain duplicate data sets, it is recommended you start from a clean directory, move these files to another directory, or rename the run for each execution of `lcout` or `demcmc`. 

1. Lightcurve file (lc_NAME.lcout), where NAME is the name specified in the input file.  
   
   This file lists the times of output and theoretical output as well as the measured flux and uncertainties. Useful for plotting the best fit.
   4 columns:
   ```
   [time (days)] \t [measured flux (normalized)] \t [model flux (normalized)] \t [measured uncertainty (normalized)]
   ```
2. Transit times files (tbv_XX_YY.out), where XX is the body being transited and YY is the body doing the transiting. Indexes are 00=star, 01=1st planet, 02= 2nd planet, etc.  

   4 columns:
   ```
   [transit index from epoch] \t [transit mid-time (days)] \t [closest approach of centers of both bodies (AU)] \t [transit velocity (AU/day)] 
   ```
 
3. Coordinate Conversion Files
   By default, the lcout command produces transformations of the inputted initial conditions to different coordinate systems. These are divided into two classes (those based on orbital elements, and those based on cartesian coordinates) and outputted into two files: 
   * xyz_NAME.xyzout
   *  aei_NAME.aeiout
   Descriptions within these files specify the coordinate system used for each entry within. 

4. [Optional] RV Output Files, (rv_XX_NAME.rvout), where XX is the index of the body for which RVs are reported (00 = the central star) and NAME is the run name. 
   Format is: 
   ```
   [time (days)] [Measured RV (m/s)] [Modeled RV (m/s)] [Uncertainty (m/s)] [Telescope Index]
   ```  
   The calculated RV offset for each telescope is printed at the bottom of the file.

5. [Optional] Celerite Continuum Fit  

   This file shows the Gaussian Process the celerite fit to the continuum if a celerite fit was chosen. It's format is 
   ```
   [time] [flux]
   ```

6. [Optional] Positions of the Bodies at time = PRINTEPOCH, named xyz_adjusted_NAME.xyzout, where NAME is the run name. Similar to xyz_NAME.xyzout, but at a different epoch time, and only one coordinate system is given. 


### MCMC Fit (`demcmc`) Output Files

1. Current best-fit solution.   

   Every time an MCMC walker encounters a parameter set with highest likelihood of any fit yet found by all of the chains, it outputs it to the file `mcmc_bestchisq_NAME.aei`, where NAME is the run name as specified in the input file. 

2. Chain state

   Every 100 generations, the current state of the chains is printed in the fitting basis to `demcmc_NAME.out`.  

3. Stepsize diagnostics

   Every 10 generations, the fraction of new parameter proposals that are accepted and the current scale factor on the differential evolution vector used to propose new steps is printed. A file is generated with name `gamma_NAME.txt` and format:
   ```
   [generation] \t [proposal acceptance fraction] \t [scale factor]
   ```

4. Diagnostic Output

   Some diagnostic output is saved in `demcmc.stdout` and/or the `outf` subdirectory, depending on how the run is executed. If the run crashes or hangs, these are the best places to look for help.  
 


## Example Systems

Example data files, setup scripts, etc., are included in /runs


## Analysis Tools

Simple plotting and analysis tools written in Python are available in `/src/tools`. These tools require the `numpy`, `pandas`, `matplotlib`, `corner` and possibly other packages. They are designed for use with python 3.5. 

1. Forward Model Plotting

   After a forward model has been run, these scripts provide quick diagnostic plots of the fits and transits. The scripts should be copied to the directory that the forward model was run in as they assume the output files are in the current directory. They may be run with
   ```
   $ python scriptname.py
   ```
   * `lcplot.py` - This script plots a segment of the lightcurve. It overplots the model and the data, and includes the residuals of the fit at the bottom. The modeled location of each planet is indicated at the top of the figure in a unique color. To change the range of the data plotted, the first line of the script may be edited. The script produces a figures titled 'lcplot.png.' A "gp" version can be used that includes the Guassian Process model from celerite.

   * `omc.py` - This script computes the mean period of the planets over the data range by performing a linear fit to the model's transit times. The difference between the modeled transit times and this constant linear ephemeris (period) is known as an O-C (Observed minus calculated) TTV diagram. The O-C for each planet is shown in the outputted figures 'omc_AA_BB.png', where AA indicates the body being transited (the star = 00 in most cases) and BB is the planet index (01 for the innermost planet).  
       
   * `phasefold.py` - This script plots all transits of each planet that do not overlap with any other transits. The innermost (top panel) to outermost (bottom panel) planets are over-plotted by removing their TTVs as computed in the model. The raw data is shown in gray points and black points show the data phase-binned in 15 minute intervals. 

2. DEMCMC Analysis

   * `demcmc_quick_analyze.py` is a python script to perform some standard analysis of the output of a demcmc run. It should be copied to the directory where a demcmc_RUNNAME.out file was created, and should be invoked with:
      ```
      $ python demcmc_quick_analyze.py INFILE.in [burnin]
      ``` 
      where INFILE.in is the name of the .in input file for the run and burnin is an optional parameter which can be set to an integer N to disregard the first N steps when computing posterior information. Running this script produces a new subdirectory called analysis_dir which is populated with:
      
      * Trace plots for each parameter (format=.png)
      * Corner correlation plots for all parameters (format=.png)
      * 1- & 3-sigma confidence intervals for each parameter (format=.txt). These use the median and [.16, .84] and [.0015, .9985] percentiles. 
      * 2-sigma upper limits for each parameter, i.e., the value which 95% of draws lie below (format=.txt). This is useful for understanding the upper bound on paramters like mass, but may be meaningless for other parameters which are more well-defined. 
      * Corner correlation plots for the masses and eccentricities (transformed out of the fitting basis) for each planet (format=.png)
      * 1-d posteriors of all planets m mass, radius, and densities marginalized over all other parameters (format=.png)
      * Gelman-Rubin Rhat statistics for gauging MCMC convergence for each parameter (format=.txt)

   * `make_ttv_posterior_cloud.sh` is a bash script to help visualize a TTV posterior from the output of a demcmc run. It should be copied to the directory where a demcmc_RUNNAME.out file was created, and should be invoked with:
      ```
      $ ./make_ttv_posterior_cloud.sh INFILE.in 
      ```
      It may be necessary to change the first line of the script which specifies the location of bash for the script to run properly. It will create 
a subdirectory structure where individual draws from the posterior are made, modeled, and their TTV output analyzed. The resulting figures can be 
found in `./cloudttv_dir/posterior_draws` and are called `tbv_XX.out.pdf`. These are plots of the TTVs drawn from the posteriors (gray circles), and 
their mean and standard deviation (blue). A few tunable parameters are listed at the head of the `make_ttv_posterior_cloud.sh` script.






