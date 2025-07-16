# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use('agg')  # set this up to be able to save figures without X tunneling
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd
import fileinput
import sys
import csv

TTVDIR="/home/danielkj/nobackup/autodelete/PhoDyMM_v3/ttv_data/"
MAPFILE = "/home/danielkj/nobackup/autodelete/PhoDyMM_v3/koiprops_20230830_v2.csv"


# get TTV data for a planet with an approximate period "per" and specific "kic" value
# returns an array of all transit times with 3 columns: time, TTV, and errors (in days)
# all np.zeros(100,3) are returned if we can't figure this out
def get_ttv_data(per,koiName):
  koiData = pd.read_csv(MAPFILE)
  koiData = koiData.to_numpy()
  
  koiNum = int(koiName[3:])
  indices = []
  for i in range(len(koiData)):
      if int(koiData[i][1]) == koiNum:
          indices.append(i)
  
  koiIndex = -1
  for index in indices:
      koiPer = koiData[index][3]
      if abs(koiPer-per) < 0.05 * per:
          koiIndex = index
  
  koiName_full = str(koiData[koiIndex][1])
  koiName_full = 'koi' + koiName_full.rjust(7,'0')
  
  koifilename = TTVDIR + koiName_full + '.tt'
  
  print(koifilename)
  
  try: 
    ttvdata=np.loadtxt(koifilename,ndmin=2)
  except:
    ttvdata=np.zeros((1,1))
  if len(ttvdata) < 3:
    # didn't work for some reason
    ttvdata=np.zeros((100,3))

  # number of ttvs vs. [time, ttv, ttverr] in days
  return ttvdata  

"""
Try to find if the data is actually being grabbed
"""

# read in a tbv out file produced by PhoDyMM's lcout
def read_in_tbv(fname):
  data = np.loadtxt(fname,ndmin=2)
  n = data[:,0]
  tt = data[:,1]
  sorti = np.argsort(n)
  n = n[sorti]
  tt = tt[sorti]
  return n, tt

# compute the Observed minus Calculated for the model (in minutes)
def compute_omc(n, tt):
  # O-C requires removable of the best fit linear model
  fit = np.polyfit(n, tt, 1)
  b = fit[0] # best fit period
  a = fit[1] # best fit epoch
  c = b*n + a # best fit model
  
  # calculate difference between transit time and best-fit line 
  # and convert to minutes
  omc = 24.*60.*(tt-c)

  return tt, omc, b  # transit times, O-C, and period


# make the main O-C plot
# take in transit times, O-C, period, and optional output filename
def plot_omc(tt, omc, b, koiName,fname=None):
  if fname is None:
    fname = 'omcd.png' # default file name

  # calculate the stadard deviation for putting on the plot 
  sigma = np.std(omc/60./24.)

  # set up figure
  f = plt.figure()
  ax = f.add_subplot(111)
  # Plot transit times and O-C=TTV
  ax.scatter(tt, omc)

  # estimate period and gather TTV data
  per = np.median(np.diff(tt))
  ttvdata=get_ttv_data(per,koiName)
  
  # make the plot with error bars of the TTV data
  # when there's no TTV data, this is just a single point at (0,0)
  plt.errorbar(ttvdata[:,0],ttvdata[:,1]*1440.0,ttvdata[:,2]*1440.0,capsize=0,alpha=0.5,fmt='r.')
  # scale the plot to cover Model, not data
  plt.ylim([min(omc)*1.5,max(omc)*1.7])

  # Annotate the plot
  plt.text(0.02, 0.95, "Average Period = %.5lf days" % b, transform=ax.transAxes)
  plt.text(0.02, 0.90, "Sigma = %.5lf days" % sigma, transform=ax.transAxes)
  plt.xlabel('Time (days)')
  plt.ylabel('TTV (minutes)')
  f.savefig(fname)
  plt.close('all')
  
###########################################################################
# This section is for outputting a tdv file, previous sections were for ttv
###########################################################################



# This function reads in all of the tbv, as they are all necessary for tdvs
# The columns are then sorted because the negatives were at the back
def read_in_tbv2(fname):
  data = np.loadtxt(fname,ndmin=2)
  n = data[:,0]
  tt = data[:,1]
  b = data[:,2]
  v = data [:,3]
  sorti = np.argsort(n)
  n = n[sorti]
  tt = tt[sorti]
  b = b[sorti]
  v = v[sorti]
  return n, tt, b, v

# This function pulls in radius of the star in solar radii from the first pldin
def read_in_Rstar():
    Rstar = 0
    #pldinFileList = glob.glob("./*.pldin")
    pldinFileList = pldinFile
    #if len(pldinFileList) > 1:
    #    print ("Warning: There is more than one pldin! Assuming the first:")
    #    print(pldinFileList[0])
    with open(pldinFileList) as infile:
        for line in infile:
            if "Rstar" in line:
                Rstar = line.split()[0]
            if "rstar" in line:
                Rstar = line.split()[0]
        infile.close()
    return Rstar

def read_in_rpors(fname):
    planet = 0
    rpors = 0.0
    if fname[-4:] != ".out":
        print ("Warning, trying to read in tbv but it is not a .out file")
    planet = fname[-5:-4]
    #pldinFileList = glob.glob("./*.pldin")
    pldinFileList = pldinFile
    with open(pldinFileList) as infile:
        flag = "0." + str(planet)
        for line in infile:
            if flag in line:
                parameters = line.split()
                # this is necessary so it doesnt try to grab something like limb darkening
                # or the rpors of a different line because some other number started
                # with 0.1
                if len(parameters) == 9:
                    if parameters[0] == flag:
                        rpors = parameters[8]
        infile.close()
    return rpors

def compute_durations(Rstar, rpors, b, v):
    durations = list() # maybe try and array and push back?
    # debugging
    for i in range(len(b)):
        #sqrt((Rs(1+rpors)^2 +b^2))/ v^2
        #d = ((((Rstar * (1 + rpors))^2 +b[i]^2) ^ (1 / 2)) * 2) / v[i]^2
        #24 multiplied to convert days into hours
        #print(Rstar/0.00465047)
        #print(rpors)
        #print(i, b[i], v[i])
        #print(i, b[i]/0.00465, v[i]/0.00465)
        #print(np.sqrt(Rstar*Rstar - b[i]*b[i])/v[i]*48)
        #print()
        d = 24.0 * 2 * ((((Rstar * (1.0 + rpors)) ** 2.0) - (b[i] ** 2.0)) ** (1.0 / 2.0)) / (v[i])
        #print(d)
        #print()
        durations.append(d)
    return durations

def plot_tdv(tt, d, fname=None): # pass what you need
  if fname is None:
    fname = 'tdv.png' # default file name

  # calculate the stadard deviation for putting on the plot 
  #sigma = np.std(omc/60./24.)

  # set up figure
  f = plt.figure()
  ax = f.add_subplot(111)
  # Plot transit times and O-C=TTV
  ax.scatter(tt, d)

  # estimate period and gather TTV data
  #per = np.median(np.diff(tt))
  #ttvdata=get_ttv_data(per)
  
  # make the plot with error bars of the TTV data
  # when there's no TTV data, this is just a single point at (0,0)
  #plt.errorbar(ttvdata[:,0],ttvdata[:,1]*1440.0,ttvdata[:,2]*1440.0,capsize=0,alpha=0.5,fmt='r.')
  # scale the plot to cover Model, not data
  #plt.ylim([min(omc)*1.5,max(omc)*1.7])

  # Annotate the plot
  #plt.text(0.02, 0.95, "Average Period = %.5lf days" % b, transform=ax.transAxes)
  #plt.text(0.02, 0.90, "Sigma = %.5lf days" % sigma, transform=ax.transAxes)
  ax.set_xlabel('Time (days)')
  ax.set_ylabel('Duration (hours)')
  plt.tight_layout()
  f.savefig(fname)
  plt.close('all')


def save_data(fileName,time,variation,parameter):
    output = [['Time',parameter]]
    
    for i in range(len(time)):
        row = [time[i],variation[i]]
        output.append(row)
        
    
    with open(fileName, 'w+',newline='') as f:
        file_writer = csv.writer(f,delimiter =',',quotechar = "'")
        for i in range(len(output)):
            file_writer.writerow(output[i])

# expecting input files like lcout
# print(sys.argv)

inFileName = glob.glob('*_cadence.in')[0]
koiName = inFileName.split('_')[0]

bestPldinPath = koiName + '_best_parameters.pldin'
pldinPath = koiName + '.pldin'

if os.path.exists(bestPldinPath):
    pldinFile = bestPldinPath
else:
    pldinFile = pldinPath 

for f in glob.glob("./tbv[0-9][0-9]_[0-9][0-9].out"):
    n, tt, b, v = read_in_tbv2(f)  # read in the model results from lcout
    r = float(read_in_Rstar()) * 0.00465047 # converts Rstar into AU
    k = float(read_in_rpors(f)) # reads rpors as a number instead of a string
    durations = compute_durations(r, k, b, v) # computes a list of the durations
    fname = 'tdv_'+f[5:10]+'.png' # make a filename for the output to keep each planet separate
    plot_tdv(tt, durations, fname=fname) # make the plot, including current data
    save_data('tdv_'+f[5:10]+'.csv',tt,durations,'Duration')
 
  ###########################################################################



# main code
# grab all tbv out files (using glob) and make a plot for each
# using the above functions
for f in glob.glob("./tbv[0-9][0-9]_[0-9][0-9].out"):
  n, tt = read_in_tbv(f)  # read in the model results from lcout
  tt, omc, b = compute_omc(n, tt) # calculate omc values
  fname = 'omcd_'+f[5:10]+'.png' # make a filename for the output to keep each planet separate
  plot_omc(tt, omc, b,koiName, fname=fname) # make the plot, including current data
  save_data('omcd_'+f[5:10]+'.csv',tt,omc,'OMC')

