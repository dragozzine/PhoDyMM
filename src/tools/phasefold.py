# phasefold.py
# LAST LAST UPDATED:
# 12 Aug 2022 by Matthew Perez
# LAST UPDATED:
# 22 July 2020 by Sakhee Bhure
# Provides diagnostic plots for output lightcurves from PhoDyMM
# Originally written by Sean Mills. Updated by Darin Ragozzine, Sakhee Bhure, and Matthew Perez
# Designed for python 3. 


import numpy as np
import matplotlib
import matplotlib.axes as axes
matplotlib.use('agg')  # set this up so that it can output plots even without Xwindows
import matplotlib.pyplot as plt
import glob
import sys
import os
import scipy
import bisect
import logging
from scipy import stats
import time as tm




psuedo1 = 0.99999 # we are effectively using this as == 1


# Creating log file
try:
    os.remove('analysis_dir/phasefold.txt')
    file = open('analysis_dir/phasefold.txt', 'w')
except:
    file = open('analysis_dir/phasefold.txt', 'w')


# if no runname is given grab all possibilities in current directory
if len(sys.argv) != 2:
    lcdatafilearr = glob.glob("./lc_*.lcout")
else:
    # if a runname is given, grab it
    lcdatafilearr=glob.glob("./lc_*"+sys.argv[1][0:-3]+"*.lcout")[0]

# automatically takes the most recent lcout file that matches the above
lcdatafile=max(lcdatafilearr, key=os.path.getctime)

print("phasefolding this data file: ", lcdatafile)
lcdata = np.loadtxt(lcdatafile)
time = lcdata[:,0] # first column in the lcout file is the time (PhoDyMM units, BKJD-67 I think)
meas = lcdata[:,1] # second column is the MEASured/observed data from Kepler 
# which data are used isn't preserved directly but was given to the lcout command to make this file
the = lcdata[:,2] # third column is the THEoretical/model lightcurve from PhoDyMM
# the model used to generate this isn't preserved directly but was given to the lcout command to make this file
err = lcdata[:,3] # the ERRors=uncertainites from Kepler (same source as "meas").



# Making sure file was fed exactly one lcout file
uniqtime = np.unique(time) # np.unique removes duplicates
if len(uniqtime) == len(time): 
    pass
else:
    file.write('Warning: Might have multiple outputs in lcout\n')
    time = time[:-len(uniqtime)]
    meas = meas[:-len(uniqtime)]
    the = the[:-len(uniqtime)]
    err = err[:-len(uniqtime)]
    file.write('Only using last lcout output\n')

if len(lcdatafilearr) != 1:
    file.write('URGENT WARNING: no lcout file')


# used to calculate a truncated chi-square, a robust estimate of goodness of fit
truncpercentile=98
resids=((meas-the)/err)**2
#print(np.percentile(resids,95))
notoutliers=np.where(resids < np.percentile(resids,truncpercentile))
truncsum=sum(resids[notoutliers])
arcs=truncsum/len(resids[notoutliers])



# grab the tbv out files to get the true transit times
# these are typicallyoverwritten every time, so it 
# might not be correct if you are using an old lcout with 
# the new tbv.

if len(lcdatafilearr) > 1:
    file.write("Warning: lcout file and tbvout files may not match!\n")
tbvfilelist = glob.glob("./tbv[0-9][0-9]_[0-9][0-9].out")
nfiles = len(tbvfilelist)
npl = nfiles # number of planets

# initialize the figure with one column for each planet, first row is phasefold, second row is corresponding residuals
f, (axes1, axes2) = plt.subplots(2, npl, figsize=(4*npl, 7)) # 2 is number or rows, npl is number columns
axes1 = list(axes1)
axes2 = list(axes2)
axes = [axes1, axes2]


# Failsafe incase file is fed no tbv files        
goodFiles = []        
for i in range(len(tbvfilelist)):  
    if os.stat(tbvfilelist[i]).st_size == 0:
        file.write('Warning: no tbv files recieved')
    else:
        goodFiles.append(i)


# read in the transit time data for each planet
transtimes = [[] for i in range(npl)] # transit times - list of arrays (1 array is 1 planet)
nums = [[] for i in range(npl)] # number of measurement
for i in range(len(goodFiles)):
    data = np.loadtxt(tbvfilelist[goodFiles[i]],ndmin=2)
    transtimes[i] = data[:,1]
    nums[i] = data[:,0]

    
plt.close("all") # closing all plots that have opened in previous iterations of the code
        

phasewidth = [0.4 for i in range(nfiles)]
for i in range(len(goodFiles)):
    # gather values to be plotted here
    phases = [] # "phased" value for this planet in units of days
    fluxes = [] # observed/measured/Kepler flux
    models = [] # theoretical/model fluxes from PhoDyMM
    res = [] # Scaled Residuals
    phasetimes = []
    
    if len(transtimes[i]) > 1:
      # Calculate meanper (mean period)
      meanper, const = np.linalg.lstsq(np.vstack([nums[i], np.ones(len(nums[i]))]).T, transtimes[i], rcond=None)[0] 
      # meanper - mean period (slope of line) [line of lstsq^^^^]    const - y-intercep
      # 3x the duration of an edge-on planet around the sun
      phasewidth[i] = 3.*(13./24.) * ((meanper/365.25)**(1./3.)) # width of plot
      meanper = "{:.3f}".format(meanper) # meanper but with only 3 decimals, just looks better output this way
    
    
    # gather transit times of all other planets
    othertts = transtimes[:i] + transtimes[i+1:]  # transit time of all other planets
    if len(othertts) > 0:
        othertts = np.hstack(np.array(othertts, dtype=object))  # transit time of planet in question
    thistts = np.array(transtimes[i])    
    
    # Calculates collision width
    closest = []
    for tti in thistts:
        closest.append(min(abs(othertts - tti)))
    max_drop_per = 40
    
    
    max_collisionwidth = np.percentile(closest, max_drop_per) 
    collisionwidth = min([phasewidth[i]/2, max_collisionwidth])    

    # for each transit time of this planet
    for tti in thistts:
        # make list of indices to check period within range 
        phasewidthi = phasewidth[i]
        L = bisect.bisect_left(time, (tti - phasewidthi), lo=0, hi=len(time)) # index of beginning of transit
        R = bisect.bisect_right(time, (tti + phasewidthi), lo=L, hi=len(time)) # index of end of transit
        trange = range(L, R) # trange = range of times being plotted
        
        # without other transits, just gather the phases, fluxes, model, and residuals
        if len(othertts) == 0:
            phases.append(time[trange] - tti)
            fluxes.append(meas[trange])
            models.append(the[trange])      
            # with other transits, gather data if there is no collision
            res.append((meas[trange] - the[trange]) / err[trange])
        elif min(abs(othertts - tti)) > collisionwidth:  # checks no other transits are 'close' or 'colliding'
            phases.append(time[trange] - tti)
            fluxes.append(meas[trange])
            models.append(the[trange])
            res.append((meas[trange] - the[trange]) / err[trange])
            phasetimes.append(time[trange])
            

    # reorganize into a plottable numpy array
    phases = np.hstack(phases)
    fluxes = np.hstack(fluxes)
    models = np.hstack(models)
    res = np.hstack(res)
    phasetimes = np.hstack(phasetimes)
    
    # organize all arrays by phase
    phasesort = np.argsort(phases)
    phases = phases[phasesort]
    fluxes = fluxes[phasesort]
    models = models[phasesort]
    res = res[phasesort]
    phasetimes = phasetimes[phasesort]

    # use 10 minute bins, set up bins
    binwidth = 1./1440. * 10.
    nbins = int(2*phasewidth[i] / binwidth)
    binedges = np.arange(nbins+1, dtype=float)*binwidth - phasewidth[i]
    bincenters = np.arange(nbins, dtype=float)*binwidth - phasewidth[i] + binwidth/2.

    # calculate binned values of observed fluxes and theoretical models
    j=0 # iterates over data
    k=0 # iterates over bins
    mbinned = np.ones(nbins)
    mmbinned = np.ones(nbins)
    while j < len(phases):
        mbinvals = []
        mmbinvals = []
        while phases[j] < binedges[k+1]:
            mbinvals.append(fluxes[j])
            mmbinvals.append(models[j])
            j += 1
            if j >= len(phases):
                break
        if len(mbinvals) > 0:
            mbinned[k] = np.mean(mbinvals)
            mmbinned[k] = np.mean(mmbinvals)
        else:
            mbinned[k] = np.nan
            mmbinned[k] = np.nan
        k += 1
        if k >= nbins:
            break
            
    
    # Calculating the model phasewidth  
    try:
        transit_min = np.where(abs(phases) == min(abs(phases)))
        if abs(models[transit_min[0]]) < psuedo1:
            lst_of_ind = []
            right_index = 0
            for r in range(len(models)):
                if r < transit_min[0]:
                    if models[transit_min[0]-r] > psuedo1: # and models[x[0]-r-1] < 0.99999:
                        lst_of_ind.append(transit_min[0] - r)
                elif r > transit_min[0]:
                    if models[r] > psuedo1:
                        right_index = r 
                        break
            left_index = lst_of_ind[0]
            model_phasewidth = phases[right_index] - phases[left_index] # width of transit
            # Check to make sure model_phasewidth (transit width) is < collisionwidth
            if model_phasewidth >= 2 * collisionwidth:
                file.write('URGENT WARNING: Model Phasewidth is wider than Collisionwidth!\n')   
            if model_phasewidth < 0.01:
                file.write('Warning: phasewidth is too small, something went wrong')
                model_phasewidth = [1]
        else:
            file.write('\nWarning: Planet ' + str(i+1) + ' is lost in Kansas (Period = ' + str(meanper) + ')')             
            model_phasewidth = [1]  
    except:
        file.write('\nWarning: Planet ' + str(i+1) + ' has no minimum model value (Period = ' + str(meanper) + ')')              
        model_phasewidth = [1]       
                            
    
    # Output Metrics to phasefold.txt
    file.write('\nPeriod = ' + str(meanper) + ' Days' + '\nphasewidth: ' + str(phasewidth[i]) + '\nmodel_phasewidth: '
              + str(model_phasewidth[0]) + '\ncollisionwidth: ' + str(collisionwidth) + '\n')
    
    
    # plot observed and model data (Phasefold Plot)
    axes1[i].scatter(phases, fluxes, s=0.01, color='grey', alpha=0.5) # Kepler data
    axes1[i].scatter(phases, models, s=0.02, color='m', alpha=0.5)  # PhoDyMM Theoretical data
    

    # add the binned values to the plot
    axes1[i].scatter(bincenters, mbinned, s=1.0, color='k') 
    axes1[i].scatter(bincenters, mmbinned, s=0.8, color='r')   
    
    # Set vertical lines at collision width
    axes1[i].axvline(x = collisionwidth)
    axes1[i].axvline(x = -collisionwidth)
#    These plot the model_phasewidth as vertical red bars - used during creation of file
#    axes1[i].axvline(x = model_phasewidth[0]/2, color = 'r') 
#    axes1[i].axvline(x = -model_phasewidth[0]/2, color = 'r')
    
    # Set plot scale and labels 
    axes1[i].set_xlim((-phasewidth[i], phasewidth[i]))
    axes1[i].set_xlabel('Phase (days)')
    axes1[i].set_title('Period = ' + str(meanper) + ' Days')
    # Set Y limits as model +1.25*depth_of_transit and model -1.5*depth_of_transit
    ylimit = 1 - min(models[abs(phases) < model_phasewidth[0]/2]) 
    
    if max(abs(mbinned)) > 1 + (1.5 * ylimit):
            file.write('Plot of planet ' + str(i+1) + ' axes are scaled differently (Period = ' + str(meanper) + ')\n')
            ylimit = max(abs(1 - mbinned))
            axes1[i].set_ylim(1 - ylimit, 1 + ylimit)
    else:
        axes1[i].set_ylim(1 - (1.5 * ylimit), 1 + (1.25 * ylimit))

    axes1[i].set_ylabel('Normalized Flux')
    
    # add scaled residual plots
    axes2[i].scatter(phases, res, s=0.01, c='b', alpha=1) # Scaled Phase Residuals
    axes2[i].set_xlim((-phasewidth[i], phasewidth[i]))
    axes2[i].set_xlabel('Phase (days)')
    axes2[i].set_ylabel('Residuals')


# for PhoDyMM analysis, the ldcatafile is typically the KIC number for the star

#kic=lcdatafile[5:-6].split('_')[0]
kic=lcdatafile.split('_')[-1][:-6]


plt.title(kic+' and rcsi is '+str(round(arcs,3)))

f.tight_layout()
f.savefig('analysis_dir/PhaseFolded_'+kic+'.png')



# BINNED SUMMARIES

nfmbins = int(30)
epsilon = 0.001

modelwidth = max(the) - min(the)
fmbinwidth = modelwidth / nfmbins
fmbinedges = min(the) + np.arange(nfmbins+1, dtype=float)*fmbinwidth
fmbincenters = min(the) + np.arange(nfmbins, dtype=float)*fmbinwidth + fmbinwidth/2.

fmbinedges[0] = fmbinedges[0]-epsilon
fmbinedges[30] = fmbinedges[30]+epsilon

def mad(arr):
    med=np.median(arr)
    return np.median(np.abs(arr-med))

j1=0
k1=0

np.sort(the)
modelsort = np.argsort(the)
fluxes1 = meas[modelsort]
models1 = the[modelsort]
mfmbinned = np.zeros(nfmbins)
fmerr = np.zeros(nfmbins)

ksbins = []

while j1 < len(models1):
    mfmbinvals = []
    while models1[j1] <= fmbinedges[k1+1]:
        mfmbinvals.append(fluxes1[j1])
        j1 += 1
        if j1 >= len(models1):
            break
    if len(mfmbinvals) > 0:
        mfmbinned[k1] = np.median(mfmbinvals)
        fmerr[k1] = mad(mfmbinvals)
        ksbins.append(stats.kstest(mfmbinvals, 'norm'))
    else:
        mfmbinned[k1] = np.nan
    k1 += 1
    if k1 >= nfmbins:
        break

# print(ksbins)
# Now make the flux-model comparison plot

aline = np.linspace(min(meas), max(meas), len(meas))
plt.figure()

plt.plot(the,meas,'.',markersize=1.,alpha=0.2, color='steelblue')
plt.plot(aline,aline,'r-', alpha=0.5)
plt.errorbar(fmbincenters, mfmbinned, yerr=fmerr, xerr=None, fmt='.',markersize=0.5, linewidth=0.35, capsize=0., linestyle="None", color='k')

plt.xlabel('Modeled Flux')
plt.ylabel('Observed Flux')
#plt.title(lcdatafile[5:-6].split('_')[0])
plt.title('Observed vs Modeled Flux: '+kic)

#plt.savefig('FluxModelComp_'+lcdatafile[5:-6]+'.png')
plt.savefig('analysis_dir/FluxModelComp_'+kic+'.png')

#plt.savefig('FluxModelComp_'+lcdatafile[5:-6]+'.pdf')
#plt.savefig('FluxModelComp_'+kic+'.pdf')


# quality of fits
residuals=(meas-the)/err

# K-S test
stats.kstest(residuals, 'norm')

bline=np.zeros(len(meas)) # for the horizontal line
plt.figure()
plt.xlabel('Observed Flux')
plt.ylabel('Residuals')
plt.title('Residuals vs Observed Flux: '+kic)

plt.plot(meas,residuals, '.',markersize=1.0,alpha=0.2,c='g')
plt.plot(meas,bline,'r-')# horizontal line at y = 0
plt.savefig('analysis_dir/ResidualMeas_'+kic+'.png')


# plotting residuals vs theoretical
# plotting it separately for now, needs to be combined with Residuals vs Measured if it looks good

cline=np.zeros(len(the))

plt.figure()
plt.xlim(min(the),1.001) # earlier, the plot showed something beyond 1.000, now it's blank...?
plt.xlabel('Theoretical Flux')
plt.ylabel('Residuals')
plt.title('Residuals vs Theoretical Flux: '+kic)

plt.plot(the,residuals,'.',markersize=1.0,alpha=0.2,c='teal')
plt.plot(the,cline,'r-') # horizontal line
plt.savefig('analysis_dir/ResidualThe_'+kic+'.png')


#plt.show()
