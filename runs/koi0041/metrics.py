# metrics.py
# LAST LAST UPDATED:
# 12 Aug 2022 by Matthew Perez
# LAST UPDATED:
# 22 July 2020 by Sakhee Bhure
# Provides 

# Originally written by Sean Mills. Updated by Darin Ragozzine, Sakhee Bhure, and Matthew Perez
# Designed for python 3. 


import numpy as np
import glob
import sys
import os
import scipy
import bisect
from scipy import stats
import time as tm




psuedo1 = 0.99999 # we are effectively using this as == 1

# Creating .csv file for metrics
try:
    os.remove('analysis_dir/metrics.csv')
    file_csv = open('analysis_dir/metrics.csv', 'w')
except:
    file_csv = open('analysis_dir/metrics.csv', 'w')
    
    
# Creating .txt file for warnings
try:
    os.remove('analysis_dir/metrics.txt')
    file = open('analysis_dir/metrics.txt', 'w')
except:
    file = open('analysis_dir/metrics.txt', 'w')
  

# Set up .csv file
file_csv.write('Period,Mean-(in_transit_res),Mean-(out_transit_res),Standard_Deviation-(in_transit_res),Standard_Deviation-(out_transit_res),Correlation_Coefficient_Residuals_and_Flux,P-value_Residuals_and_Flux,Correlation_Coefficient_Residuals_and_Model,P-value_Residuals_and_Model\n')

    
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


# Making sure file was only fed one lcout file
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


phasewidth = [0.4 for i in range(nfiles)]
for i in range(len(goodFiles)):
    # gather values to be plotted here
    phases = [] # "phased" value for this planet in units of days
    fluxes = [] # observed/measured/Kepler flux
    models = [] # theoretical/model fluxes from PhoDyMM
    res = [] # Scaled Residuals
    
    if len(transtimes[i]) > 1:
      # Calculate meanper (mean period)
      meanper, const = np.linalg.lstsq(np.vstack([nums[i], np.ones(len(nums[i]))]).T, transtimes[i], rcond=None)[0] 
      # meanper - mean period (slope of line) [line of lstsq^^^^]    const - y-intercep
      # 3x the duration of an edge-on planet around the sun
      phasewidth[i] = 3.*(13./24.) * ((meanper/365.25)**(1./3.)) # width of plot
      meanper = "{:.3f}".format(meanper) # meanper but with only 3 decimals
    
    
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
    
    max_collisionwidth = np.percentile(closest, max_drop_per) # try 100 - max_drop_per
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
            

    # reorganize into a plottable numpy array
    phases = np.hstack(phases)
    fluxes = np.hstack(fluxes)
    models = np.hstack(models)
    res = np.hstack(res)

    # organize all arrays by phase
    phasesort = np.argsort(phases)
    phases = phases[phasesort]
    fluxes = fluxes[phasesort]
    models = models[phasesort]
    res = res[phasesort]
                
    
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
                file.write('\nURGENT WARNING: Model Phasewidth is wider than Collisionwidth!\n')  

            # Differentiating residuals: whether they are in-transit or out-of-transit residuals
            in_transit_res = [] # list of in transit residuals
            out_transit_res = [] # list of out of transit residuals
            in_transit_flux = [] # list of in transit fluxes
            in_transit_model = [] # list of in transit models
            for m in range(int(left_index), int(right_index)):
                in_transit_res.append(res[m])
                in_transit_flux.append(fluxes[m])
                in_transit_model.append(models[m])

            for n in range(len(models)):
                    if models[n] >= 1:
                        out_transit_res.append(res[n])


            # Calculates the mean and St.D of in and out of transit residuals           
            if len(in_transit_res) >= 1:
                mean_in_transit = np.mean(in_transit_res)
                std_in_transit = np.std(in_transit_res)

            if len(out_transit_res) >= 1:
                mean_out_transit = np.mean(out_transit_res)
                std_out_transit = np.std(out_transit_res)


            # Calculate Correlation coeff between in_transit_res and in_transit_flux/in_transit_model and their pvalue's
            corr_flux = scipy.stats.spearmanr(in_transit_res, in_transit_flux)
            corr_model = scipy.stats.spearmanr(in_transit_res, in_transit_model)    


            # Outputing to metrics.csv
            file_csv.write(str(meanper) + ',' + str(mean_in_transit) + ',' + str(mean_out_transit) + ',' + str(std_in_transit) + 
                       ',' + str(std_out_transit) + ',' + str(corr_flux[0]) + ',' + str(corr_flux[1]) + ',' + str(corr_model[0]) + 
                       ',' + str(corr_model[1]) + '\n')

        else:
            file.write('\nWarning: Planet ' + str(i+1) + ' is lost in Kansas and metrics are tainted (Period = ' + str(meanper) + ')\n') 
            # Outputing to metrics.csv
            file_csv.write('0,0,0,0,0,0,0,0,0\n')
    
    except:
        file.write('\nWarning: Planet ' + str(i+1) + ' has no minimum model value (Period = ' + str(meanper) + ')')              
        model_phasewidth = [1]   

    
# quality of fits
residuals=(meas-the)/err

# K-S test
stats.kstest(residuals, 'norm')

