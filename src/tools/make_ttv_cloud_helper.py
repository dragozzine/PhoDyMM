import glob
import numpy as np
import matplotlib.pyplot as plt

# Generate average TTV data
ttvfiles = glob.glob("tbv00_*.out")
if len(ttvfiles) == 0:
    print("Error: This script is not meant to be called directly")
    print("       It should be invoked by make_ttv_posterior_cloud.sh")
    exit()
for ttvfile in ttvfiles:
    ttvdata = np.loadtxt(ttvfile)
    n = ttvdata[:,0]
    t = ttvdata[:,1]
    alln = []
    avgt = []
    sigt = []
    for i in range(int(np.min(n)), int(np.max(n))):
        mask = np.where(n == i)
        tmask = t[mask]
        avg = np.mean(tmask)
        sig = np.std(tmask)
        alln.append(i)
        avgt.append(avg)
        sigt.append(sig)

    np.savetxt(ttvfile+".avg", np.transpose([alln, avgt, sigt]))
     

# Plot TTV data
avgttvfiles = glob.glob("tbv00_*.out.avg")
for avgttvfile in avgttvfiles:
    avgttvdata = np.loadtxt(avgttvfile)
    avgn = avgttvdata[:,0]
    avgt = avgttvdata[:,1]
    avge = avgttvdata[:,2]
    ttvfile = avgttvfile[:-4]
    ttvdata = np.loadtxt(ttvfile)
    n = ttvdata[:,0]
    t = ttvdata[:,1]
  
    mint = min(t)
    maxt = max(t/2.)

    linreg = np.polyfit(avgn, avgt, 1, w = 1./avge)
    slope = linreg[0]
    const = linreg[1]


    omc = t - (slope * n + const)
    avgomc = avgt - (slope * avgn + const)
    plt.scatter(t, omc*1440, alpha=0.2, edgecolors='gray', marker='o', facecolors='none', s=3)
    plt.errorbar(avgt, avgomc*1440, yerr=avge*1440, color='blue', ms=0, ls='none', capsize=1, elinewidth=0.3, markeredgewidth=0.3)
    plt.axvline(x=maxt/3., color='black')
   
    plt.xlim((mint, maxt))
    plt.xlabel('Time (d)')
    plt.ylabel('TTV O-C (m)')
    plt.savefig(ttvfile+'.pdf')
    plt.clf()
 


