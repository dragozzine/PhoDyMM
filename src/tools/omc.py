import numpy as np
import matplotlib.pyplot as plt
import glob

def read_in_tbv(fname):
  data = np.loadtxt(fname)
  n = data[:,0]
  tt = data[:,1]
  sorti = np.argsort(n)
  n = n[sorti]
  tt = tt[sorti]
  return n, tt

def compute_omc(n, tt):
  fit = np.polyfit(n, tt, 1)
  b = fit[0]
  a = fit[1]
  c = b*n + a

  omc = 24.*60.*(tt-c)

  return tt, omc, b

def plot_omc(tt, omc, b, fname=None):
  if fname is None:
    fname = 'omc.png'

  sigma = np.std(omc/60./24.)

  f = plt.figure()
  ax = f.add_subplot(111)
  ax.scatter(tt, omc)
  plt.text(0.02, 0.95, "Average Period = %.5lf days" % b, transform=ax.transAxes)
  plt.text(0.02, 0.90, "Sigma = %.5lf days" % sigma, transform=ax.transAxes)
  plt.xlabel('Time (days)')
  plt.ylabel('TTV (minutes)')
  f.savefig(fname)
        

tbvflist = glob.glob("./tbv[0-9][0-9]_[0-9][0-9].out")
if len(tbvflist) == 0:
  print("Warning: no tbvXX_YY.out files found in this directory")

for f in tbvflist:
  n, tt = read_in_tbv(f)
  tt, omc, b = compute_omc(n, tt)
  fname = 'omc_'+f[5:10]+'.png'
  plot_omc(tt, omc, b, fname=fname)




