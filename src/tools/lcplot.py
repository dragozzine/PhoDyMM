### Edit this range to plot a different segment of data
trange = [100,150]
#

import numpy as np
import matplotlib.pyplot as plt
import glob

colorlist = ['b', 'r', 'g', 'y', 'c', 'm', 'midnightblue', 'yellow'] 

def plotlc(time, meas, the, err, trange, outputname=None, autoadjust=True):

  f = plt.figure()
  f, (a0, a1) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[3, 1]})
  a0.scatter(time, meas, s=0.1, c='k') 
  a0.plot(time, the, marker="None", c=colorlist[0])
  
  a1.scatter(time, meas-the, s=0.1, c='k')
  #a1.errorbar(time, meas, yerr=err, marker="None", linestyle="None", c='k')
  
  inrange = np.where((time > trange[0]) & (time < trange[1]))
  if autoadjust and (np.any(inrange) == False):
    while (np.any(inrange) == False):
      trange = [tr+50. for tr in trange]
      inrange = np.where((time > trange[0]) & (time < trange[1]))
    print("Warning: trange changed to [%7.2f, %7.2f]\n" % (trange[0], trange[1]))

  maxf0 = np.max(meas[inrange])
  minf0 = np.min(meas[inrange])
  yrange0 = [minf0, maxf0]
  maxf1 = np.max(meas[inrange]-the[inrange])
  minf1 = np.min(meas[inrange]-the[inrange])
  yrange1 = [1.1*minf1, 1.1*maxf1]
  
  a0.set_xlim((trange[0], trange[1]))
  a0.set_ylim((yrange0[0], yrange0[1]))
  a1.set_xlim((trange[0], trange[1]))
  a1.set_ylim((yrange1[0], yrange1[1]))
  
  ylim = a0.get_ylim()
  yspan = ylim[1]-ylim[0]
  
  i=0 
  for fname in glob.glob("./tbv[0-9][0-9]_[0-9][0-9].out"):
    i+=1
    data = np.loadtxt(fname)
    tt = data[:,1]
    for tti in tt:
      a0.plot([tti, tti], [ylim-0.05*yspan, ylim], c=colorlist[i], marker="None") 
  
  plt.xlabel('Time (days)')
  savestr = 'lcplot'
  if outputname is not None:
    savestr += '_'
    savestr += outputname
    savestr += '.png'
  f.savefig(savestr)

####

lcdatafile = glob.glob("./lc_*.lcout") 
if len(lcdatafile) == 0:
  print("This script must be run in the directory containing the lc_RUNNAME.lcout")
  print("    file produced with the 'lcout' command. No such file was not found here.")
  print("    Aborting")
  exit()
if len(lcdatafile) > 1:
  print("Warning: Multiple lc_RUNNAME.lcout files found in this directory")
  print("    The default behavior is to plot the first one alphabetically")

lcdata = np.loadtxt(lcdatafile[0])
time = lcdata[:,0]
meas = lcdata[:,1]
the = lcdata[:,2]
err = lcdata[:,3]

####

plotlc(time, meas, the, err, [50., 1570.], outputname='AllData')
plotlc(time, meas, the, err, trange, outputname='zoom', autoadjust=True)



