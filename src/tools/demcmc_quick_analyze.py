# DEMCMC analysis
# Usage:
# $ python demcmc_quick_analyze.py kepler36_longcadence.in [burnin]
# 
# burnin is an optional long integer of the number of generations to discard at the beginning of the demcmc output file
#

##
# additional uncertainty in stellar mass to include 
msigma = 0.0
##

import sys
import os
import pandas as pd
import numpy as np
import corner
import matplotlib.pyplot as plt
plt.switch_backend('agg')

if len(sys.argv) == 1:
  print("ERROR: You must pass this script a .in file")
  print("e.g.")
  print("$ python demcmc_quick_analyze.py kepler36_longcadence.in")
  exit()

fname = sys.argv[1]

burnin = 0
if len(sys.argv) == 3:
  burnin = int(sys.argv[2])
burnin //= 100


f=open(fname)
lines=f.readlines()

def get_in_value(linenumber):
  line = lines[linenumber]
  line = line.split()
  return line[2]

runname = get_in_value(8)
#runname = runname[:-1]
demcmcfile='demcmc_'+runname+'.out'
nchain = int(get_in_value(14))
npl = int(get_in_value(10)) - 1
parlist = lines[91+npl]
npar = len(parlist.split())
ndead=1



# lines per output in demcmc_ file
nper=npl+npar+ndead
# string to save files with
#savestr = demcmcfile.partition('demcmc_')[2]
#savestr = savestr.partition('.out')[0]
savestr = runname

print(nper, npl, npar, ndead)

# Read in file
try:
  df1 = pd.read_csv(demcmcfile, sep='\t', header=None, skiprows=lambda x: x % nper >= npl, comment=';')
except IOError:
  print("Error: This script must be run in the same directory as the 'demcmc' output file:")
  print("    "+demcmcfile)
  exit() 
except Exception:
  print("The .in file was apparently parsed incorrectly, or your demcmc_NAME.out file has been corrupted.")
  print("   The .in file parsing has obtained the following values:")
  print("   N planets = %i" % npl)
  print("   N non-planet parameters (5 stellar + jitter + GPs) = %i" % npar)
  print("   N total rows per parameter set in the demcmc_NAME.out file = %i" % nper)
  print("")
  exit() 
df2 = pd.read_csv(demcmcfile, sep='\t', header=None, skiprows=lambda x: (x % nper < npl or x % nper >= npl+npar), comment=';')


# Number of parameters for each planet
pperplan=df1.iloc[0].size

# Restructure data into dataframe that has row indexes corresponding to generation
#   and columns corresponding to different parameters
allchs=[]
for k in range(nchain):
  gens = [df1.iloc[k*npl+i::npl*nchain] for i in range(npl)]
  [i.reset_index(drop=True, inplace=True) for i in gens]
  chk = pd.concat(gens, axis=1, ignore_index=True) 
  
  gens_extra = [df2.iloc[k*npar+i::npar*nchain] for i in range(npar)]
  [i.reset_index(drop=True, inplace=True) for i in gens_extra]
  chk=pd.concat([chk]+gens_extra, axis=1, ignore_index=True)

  chk=pd.concat([chk]+[pd.DataFrame(k, index=np.arange(chk.shape[0]), columns=['Chain#'])], axis=1, ignore_index=True)

  allchs = allchs + [chk]

allparam = pd.concat(allchs)


## Construct list of column names
##pparamlist=["Planet", "Period", "T0", "sqrt(e)*cosw", "sqrt(e)*sinw", "i", "Omega", "Mjup", "Rp/Rs"]
##if pperplan > 9:
##  pparamlist += ["bright", "c_1", "c_2"]
##sparamlist=["Ms", "Rs", "c_1", "c_2", "dilute"]
##extralist=["sigma"]
## Construct list of column names with fancy text
pparamlist=["Planet", "Period", "T$_0,$", r"$\sqrt{e}\cos \omega$", r"$\sqrt{e}\sin \omega$", "i", r"$\Omega$", "M$_{jup,}$", "R$_p$/R$_s$"]
if pperplan > 9:
  pparamlist += ["bright", "c$_1$", "c$_2$"]
sparamlist=["M$_s$", "R$_s$", "c$_1$", "c$_2$", "dilute"]
extralist=["$\sigma$"]

fullplist=[]
alphabet = 'bcdefghijklmnopqrstuvwxyz'
for i in range(npl):
  fullplist += [pi+'$_'+alphabet[i]+'$' for pi in pparamlist]
fullplist += sparamlist
for i in range(npar - len(sparamlist)):
  fullplist += [a+"$_"+alphabet[i]+'$' for a in extralist] 
fullplist += ["Chain#"]


#allparam.to_csv('analysis_dir/out_4_all_preb.txt', sep='\t')
# Remove burnin
allparam.drop(range(burnin), inplace=True)

allparam.columns = fullplist

if not os.path.exists('analysis_dir'):
  os.makedirs('analysis_dir')   

# Find and print median and error bars of each parameter
sigma1o2=0.682689492137/2*100
sigma3o2=0.997300203936740/2*100
sigma2ul=0.954499736103642*100
percentiles=[50.-sigma3o2,50.-sigma1o2,50.,50.+sigma1o2,50.+sigma3o2,sigma2ul]

for p in fullplist:
  if (p[:6] != 'Planet' and p!= 'Chain#'):
    plistp = np.percentile(allparam[p], percentiles)
    with open('analysis_dir/fits_1sigma_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s\t%20.10f + %20.10f - %20.10f \n' %('{0: <16}'.format(p), plistp[2], plistp[3]-plistp[2], plistp[2]-plistp[1]))
    with open('analysis_dir/fits_3sigma_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s\t%20.10f + %20.10f - %20.10f \n' %('{0: <16}'.format(p), plistp[2], plistp[4]-plistp[2], plistp[2]-plistp[0]))
    with open('analysis_dir/fits_2sigmaUpperLimits_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s < %20.10f \n' %('{0: <16}'.format(p), plistp[5])) 


for p in fullplist:
  if (p[:6] != 'Planet' and p!= 'Chain#'):
    # compute Rubin-Gelman Statistic
    num_samples = allparam.groupby('Chain#')[p].size()[0]
    B = num_samples * allparam.groupby('Chain#')[p].mean().var(ddof=1)
    W = allparam.groupby('Chain#')[p].var(ddof=1).mean()
    Vhat = W * (num_samples - 1) / num_samples + B / num_samples
    Rhat = np.sqrt(Vhat / W)

    with open('analysis_dir/Rhat_'+savestr+'.txt', 'a') as pctfile:
      #pctfile.write('%s\t%20.10f \t %20.10f \n' %('{0: <16}'.format(p), Rhat, autocorrlen))
      pctfile.write('%s\t%20.10f\n' %('{0: <16}'.format(p), Rhat))
    

# make corner plot of parameters
colstrue = [c for c in allparam.columns if (c[:6] != 'Planet' and c != 'Chain#')]
alltrueparam = allparam[colstrue]
ranges=[(min(alltrueparam.iloc[:,i]), max(alltrueparam.iloc[:,i])) for i in range(alltrueparam.iloc[0].size)]
ranges=[(i[0],i[0]+1) if i[0] == i[1] else i for i in ranges]
figure = corner.corner(alltrueparam, range=ranges, top_ticks=True)
figure.savefig('analysis_dir/cornerplots_'+savestr+'.png')


# Produce trace plots for all parameters
fig, axes = plt.subplots(nrows=alltrueparam.shape[1], ncols=1, sharex=True, figsize=(5,1.0*alltrueparam.iloc[0].size))
pcount = 0
for p in range(allparam.shape[1]-1):
  if (fullplist[p][:6] != 'Planet'):
    for key, grp in allparam.groupby(['Chain#']):
      grp.plot(y=fullplist[p], use_index=True, ax=axes[pcount], legend=False)
      axes[pcount].text(0.5, 0.9, fullplist[p], transform=axes[pcount].transAxes, ha="center", va="center")
    pcount +=1
fig.savefig('analysis_dir/traces_'+savestr+'.png')


MEOMJ=0.00314636
for i in range(npl):
  allparam['M$_%s$' % alphabet[i]] = allparam['M$_{jup,}$$_%s$' % alphabet[i]]/MEOMJ * (1.+np.random.randn(len(allparam['M$_{jup,}$$_%s$' % alphabet[i]]))*msigma)


if npl > 4:
  fig, axes = plt.subplots(nrows=sum(divmod(npl, 2)), ncols=2, figsize=(6, 2*sum(divmod(npl,2))))
  for i in range(npl):
    allparam.hist('M$_%s$' % alphabet[i], ax=axes[i/2, i % 2], normed=1)
  axes[sum(divmod(npl,2))-1,0].set_xlabel('M$_{Earth}$', fontsize=14)
  axes[npl/2-1,1].set_xlabel('M$_{Earth}$', fontsize=14)
  fig.savefig('analysis_dir/m_'+savestr+'.png')
else:
  fig, axes = plt.subplots(nrows=npl, ncols=1, figsize=(6, 3*npl))
  for i in range(npl):
    allparam.hist('M$_%s$' % alphabet[i], ax=axes[i], normed=1)
  axes[npl-1].set_xlabel('M$_{Earth}$', fontsize=14)
  fig.savefig('analysis_dir/m_'+savestr+'.png')


REORS=0.009171
for i in range(npl):
  allparam['R$_%s$' % alphabet[i]] = allparam['R$_p$/R$_s$$_%s$' % alphabet[i]]*allparam['R$_s$']/REORS

if npl > 4:
  fig2, axes2 = plt.subplots(nrows=sum(divmod(npl, 2)), ncols=2, figsize=(6, 2*sum(divmod(npl,2))))
  for i in range(npl):
    allparam.hist('R$_%s$' % alphabet[i], ax=axes2[i/2, i % 2], normed=1)
  axes2[sum(divmod(npl,2))-1,0].set_xlabel('R$_{Earth}$', fontsize=14)
  axes2[npl/2-1,1].set_xlabel('R$_{Earth}$', fontsize=14)
  plt.tight_layout()
  fig2.savefig('analysis_dir/r_'+savestr+'.png')
else:
  fig2, axes2 = plt.subplots(nrows=npl, ncols=1, figsize=(6, 3*npl))
  for i in range(npl):
    allparam.hist('R$_%s$' % alphabet[i], ax=axes2[i], normed=1)
  axes2[npl-1].set_xlabel('R$_{Earth}$', fontsize=14)
  plt.tight_layout()
  fig2.savefig('analysis_dir/r_'+savestr+'.png')


MEG=5.9721986*10**27 # earth in grams
RECM=6371.008*10**5  # earth in cm
for i in range(npl):
  allparam[r'$\rho$$_%s$' % alphabet[i]] = allparam['M$_%s$' % alphabet[i]]*MEG / (4./3.*np.pi * (RECM*allparam['R$_%s$' % alphabet[i]])**3. ) 

if npl>4:
  fig3, axes3 = plt.subplots(nrows=sum(divmod(npl, 2)), ncols=2, figsize=(6, 2*sum(divmod(npl,2))))
  for i in range(npl):
    allparam.hist(r'$\rho$$_%s$' % alphabet[i], ax=axes3[i/2, i % 2], normed=1)
  axes3[sum(divmod(npl,2))-1,0].set_xlabel('g cm$^{-3}$', fontsize=14)
  axes3[npl/2-1,1].set_xlabel('g cm$^{-3}$', fontsize=14)
  plt.tight_layout()
  fig3.savefig('analysis_dir/rho_'+savestr+'.png')
else:
  fig3, axes3 = plt.subplots(nrows=npl, ncols=1, figsize=(6, 3*npl))
  for i in range(npl):
    allparam.hist(r'$\rho$$_%s$' % alphabet[i], ax=axes3[i], normed=1)
  axes3[npl-1].set_xlabel('g cm$^{-3}$', fontsize=14)
  plt.tight_layout()
  fig3.savefig('analysis_dir/rho_'+savestr+'.png')

for i in range(npl):
  allparam['$e_%s$' % alphabet[i]] = allparam[r'$\sqrt{e}\cos \omega$$_%s$' % alphabet[i]]**2. + allparam[r'$\sqrt{e}\sin \omega$$_%s$' % alphabet[i]]**2.



mecols=[]
for i in range(npl):
  mecols.append('M$_%s$' % alphabet[i])
for i in range(npl):
  mecols.append('$e_%s$' % alphabet[i])
print(mecols)
# corner plot es and ms
mepars = allparam[mecols]

ranges=[(min(mepars.iloc[:,i]), max(mepars.iloc[:,i])) for i in range(mepars.iloc[0].size)]
ranges=[(i[0],i[0]+1) if i[0] == i[1] else i for i in ranges]

figure3 = corner.corner(mepars, range=ranges)
figure3.savefig('analysis_dir/me_'+savestr+'.png')




# Find and print median and error bars of each parameter
sigma1o2=0.682689492137/2*100
sigma3o2=0.997300203936740/2*100
sigma2ul=0.954499736103642*100
percentiles=[50.-sigma3o2,50.-sigma1o2,50.,50.+sigma1o2,50.+sigma3o2,sigma2ul]

fullplist = allparam.columns 
print(fullplist)

for p in fullplist:
  if (p[:6] != 'Planet' and p!= 'Chain#'):
    plistp = np.percentile(allparam[p], percentiles)
    with open('analysis_dir/fits_1sigma_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s\t%20.10f + %20.10f - %20.10f \n' %('{0: <16}'.format(p), plistp[2], plistp[3]-plistp[2], plistp[2]-plistp[1]))
    with open('analysis_dir/fits_3sigma_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s\t%20.10f + %20.10f - %20.10f \n' %('{0: <16}'.format(p), plistp[2], plistp[4]-plistp[2], plistp[2]-plistp[0]))
    with open('analysis_dir/fits_2sigmaUpperLimits_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s < %20.10f \n' %('{0: <16}'.format(p), plistp[5])) 

