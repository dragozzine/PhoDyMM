# demcmc_quick_analyze.py
# Mostly written by Sean Mills
# Updates by Darin Ragozzine and students
# Designed for python 3(.5) 
#
#
# Analyzes output from PhoDyMM DEMCMC calculations to provide information and 
# diagnostics useful for understanding the outputs. (Note: many of these are
# generally useful for DEMCMC.)
#
# Usage:
# $ python demcmc_quick_analyze.py kepler36_longcadence.in [burnin]
# where the *.in file is the one that was used to produce the demcmc*.out file. 
# burnin is an optional integer (default=0) of the number of generations (BEFORE) thinning
# to remove from the beginning of the DEMCMC chain. The choice of burnin is 
# complicated (Hogg & Foreman-Mackey 2017, Ragozzine et al. 2020), but some 
# information on this can be gleaned from the outputs of this file. It is common
# to run this code with 0 burnin at first and use the diagnostic plots to 
# choose (and iterate towards) better values.
#
# Warning: if there are multiple different calculations in the same directory
# it is possible that this code will get these confused. If you are doing multiple
# runs, be sure to keep things organized.
# 


##
# additional uncertainty in stellar mass beyond PhoDyMM to include 
# valuable if you want to include additional systematic mass uncertinty
# in your results. Better to include this using PhoDyMM mass constraints.
# Will also NOT preserve mass correlation between planets appropriately. 
msigma = 0.0
##


# packages needed; can be downloaded with pip or conda
import sys
import os
import pandas as pd
import numpy as np
import corner
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# check for .in file argument
if len(sys.argv) == 1:
  print("ERROR: You must pass this script a .in file")
  print("e.g.")
  print("$ python demcmc_quick_analyze.py kepler36_longcadence.in")
  exit()

# open and read in *.in file
fname = sys.argv[1]
f=open(fname)
lines=f.readlines()


# function that gets values from *.in file for python
def get_in_value(linenumber):
  line = lines[linenumber]
  line = line.split()
  return line[2]





# use runname in input file to get appropriate demcmc output file
# this can go wrong if there are many overlapping runs, so be careful
runname = get_in_value(8)
demcmcfile='demcmc_'+runname+'.out'

nchain = int(get_in_value(14))  # number of chains = walkers used
npl = int(get_in_value(10)) - 1 # number of planets (nbodies - 1)



# identify how much burnin to remove  
if len(sys.argv) == 3:
  burnin = int(sys.argv[2])
else:
  print("No burnin given, using 0.")
  burnin=0
outputinterval = int(get_in_value(110+2*npl))   # how much the outputs were thinned
burnin = int(burnin / outputinterval) # thin the burnin appropriately



# use parameter list in the input file to figure out how many parameters there are
parlist = lines[91+npl]         
npar = len(parlist.split())
ndead=1

# lines per output in demcmc_ file
nper=npl+npar+ndead
# string to save files with; could change for more complex cases
savestr = runname

print("")
print("number of lines per output in "+demcmcfile+", number of planets, number of extra parameters, number of lines to skip : ")
print(nper, npl, npar, ndead)
print("")


# Read in demcmc output file, skipping lines relevant to planets
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

# Read in demcmc output file, skipping lines not relevant to planets
df2 = pd.read_csv(demcmcfile, sep='\t', header=None, skiprows=lambda x: (x % nper < npl or x % nper >= npl+npar), comment=';')


# Creating a third data frame to hold the chisq/likelihood values. Ben Proudfoot (BP) 05/19/20
df3 = pd.read_csv(demcmcfile, sep=',', header=None, skiprows=lambda x: (x % nper) != (nper - 1))
df3vals = df3.values
chisq = np.empty((df3.values.shape[0],df3.values.shape[1]))
for i in range(df3vals.shape[0]):
  chisq[i,0] = float(df3vals[i,0].split()[-1])
for i in range(1,df3vals.shape[1]):
  chisq[:,i] = df3vals[:,i]
df3 = pd.DataFrame(chisq[:,0])



# Number of parameters for each planet
pperplan=df1.iloc[0].size

# Restructure data into dataframe that has row indexes corresponding to generation
#   and columns corresponding to different parameters, keeping everything organized
allchs=[]
for k in range(nchain):
  gens = [df1.iloc[k*npl+i::npl*nchain] for i in range(npl)]
  [i.reset_index(drop=True, inplace=True) for i in gens]
  chk = pd.concat(gens, axis=1, ignore_index=True) 
  
  gens_extra = [df2.iloc[k*npar+i::npar*nchain] for i in range(npar)]
  [i.reset_index(drop=True, inplace=True) for i in gens_extra]
  chk=pd.concat([chk]+gens_extra, axis=1, ignore_index=True)

  # Adding another column to the df containing the likelihood/chisq vals. BP 05/19/20
  gens_ll = [df3.iloc[k::nchain]]
  [i.reset_index(drop=True, inplace=True) for i in gens_ll]
  chk=pd.concat([chk]+gens_ll, axis=1, ignore_index=True)

  chk=pd.concat([chk]+[pd.DataFrame(k, index=np.arange(chk.shape[0]), columns=['Chain#'])], axis=1, ignore_index=True)

  allchs = allchs + [chk]

allparam = pd.concat(allchs)


if int(get_in_value(69)) != 1:
  print("CALCULATIONS FOR ECCENTRICITY WILL BE WRONG because demcmc_quick_analyze assumes SQRTE= 1")
  print("see documentation and SQRTE value in *.in file")

## Construct list of column names with fancy text
# planets
pparamlist=["Planet", "Period", "T$_0,$", r"$\sqrt{e}\cos \omega$", r"$\sqrt{e}\sin \omega$", "i", r"$\Omega$", "M$_{jup,}$", "R$_p$/R$_s$"]
if pperplan > 9:
  pparamlist += ["bright", "c$_1$", "c$_2$"]

# stellar and other parameters
# doesn't handle celerite parameters?
sparamlist=["M$_s$", "R$_s$", "c$_1$", "c$_2$", "dilute"]
extralist=["$\sigma$"]

fullplist=[]
alphabet = 'bcdefghijklmnopqrstuvwxyz'
for i in range(npl):
  fullplist += [pi+'$_'+alphabet[i]+'$' for pi in pparamlist]
fullplist += sparamlist
for i in range(npar - len(sparamlist)):
  fullplist += [a+"$_"+alphabet[i]+'$' for a in extralist] 
fullplist += ["chisq"]
fullplist += ["Chain#"]


# Remove burnin
allparam.drop(range(burnin), inplace=True)

# Give the allparam dataframe the list of column names
allparam.columns = fullplist


# create the output directory
if not os.path.exists('analysis_dir'):
  os.makedirs('analysis_dir')   


# save all this work into a pandas dataframe
if os.path.isfile("analysis_dir/dqa_allparam.csv"):
  print("Warning! Overwriting dqa_allparam.csv")

allparam.to_csv("analysis_dir/dqa_allparam.csv")





# Now provide multiple files and plots to support DEMCMC analysis

# Begin by studying exactly the parameters as they are used in DEMCMC
# for DEMCMC diagnostics


"""
# Start with the Gelman-Rubin "R-hat" statistic. This statistic from 
# Gelman & Rubin 1992 has been used to test convergence by seeing whether
# multiple Markov Chains are well-mixed (using intra vs. inter chain variance).
# Conventional wisdom holds that the chains have converged when the largest
# R-hat for any parameter is less than some threshold (roughly 1.1 for pretty converged,
# 1.05 for quite converged, 1.01 for very converged).
# HOWEVER! Gelman-Rubin statistics are for independent Markov Chains which is
# NOT the case in DEMCMC since the chains "talk" to each other (ter Braak 2006). 
# It is likely that these values of R-hat are not great metrics of convergence
# and should ALWAYS be used with other convergence assessments. Because Gelman-
# Rubin values are used by the community (including past uses of PhoDyMM), we 
# leave these calculations intact for convenience.
# print("NOTE! Gelman-Rubin statistics are probably not accurate for DEMCMC convergence tests.")

grheader = 0
for p in fullplist:
  if (p[:6] != 'Planet' and p!= 'Chain#'):
    # compute Rubin-Gelman Statistic
    num_samples = allparam.groupby('Chain#')[p].size()[0]
    B = num_samples * allparam.groupby('Chain#')[p].mean().var(ddof=1)
    W = allparam.groupby('Chain#')[p].var(ddof=1).mean()
    Vhat = W * (num_samples - 1) / num_samples + B / num_samples
    Rhat = np.sqrt(Vhat / W)

    with open('analysis_dir/Rhat_'+savestr+'.txt', 'a') as grfile:
      if grheader == 0:
        grfile.write("NOTE! Gelman-Rubin statistics are probably not accurate for DEMCMC convergence tests.\n")
        grheader =1

      grfile.write('%s\t%20.10f\n' %('{0: <16}'.format(p), Rhat))

grfile.close()
""" 


# now make plots of posterior distributions of parameters actually used by DEMCMC ("true")
# as opposed to derived parmaeters added below

colstrue = [c for c in allparam.columns if (c[:6] != 'Planet' and c != 'Chain#')]
alltrueparam = allparam[colstrue]

# downselect only to columns that have variation (e.g., floating parameters)
removefixedcolumns=True #flag to remove fixed columns

if removefixedcolumns:
  goodcols=[]
  # remove columns where there's no variation
  for i in range(alltrueparam.iloc[0].size):
    if (max(alltrueparam.iloc[:,i])-min(alltrueparam.iloc[:,i]))!=0.0:
      # floating column
      goodcols.append(i)
  floatcolstrue=[colstrue[j] for j in goodcols]

floattrueparam = allparam[floatcolstrue]

# make cornerplot
if removefixedcolumns:
  ranges=[(min(floattrueparam.iloc[:,i]), max(floattrueparam.iloc[:,i])) for i in range(floattrueparam.iloc[0].size)]
  figure = corner.corner(floattrueparam, range=ranges, top_ticks=True)
  figure.savefig('analysis_dir/cornerplots_'+savestr+'.png')
else:
  ranges=[(min(alltrueparam.iloc[:,i]), max(alltrueparam.iloc[:,i])) for i in range(alltrueparam.iloc[0].size)]
  ranges=[(i[0],i[0]+1) if i[0] == i[1] else i for i in ranges]
  figure = corner.corner(alltrueparam, range=ranges, top_ticks=True)
  figure.savefig('analysis_dir/cornerplots_'+savestr+'.png')


# Produce trace plots for all parameters (including fixed)
fig, axes = plt.subplots(nrows=alltrueparam.shape[1], ncols=1, sharex=True, figsize=(5,1.0*alltrueparam.shape[1]))
pcount = 0
for p in range(allparam.shape[1]-1):
  if (fullplist[p][:6] != 'Planet'):
    for key, grp in allparam.groupby(['Chain#']):
      grp.plot(y=fullplist[p], use_index=True, ax=axes[pcount], legend=False)
      axes[pcount].text(0.5, 0.9, fullplist[p], transform=axes[pcount].transAxes, ha="center", va="center")
    pcount +=1
fig.savefig('analysis_dir/traces_'+savestr+'.png')









# Now add columns to allparam for most interesting derived parameters (mass, radius, density, eccentricity)
# and make corresponding plots


MEOMJ=0.00314636 # conversion factor for Jupiter to Earth masses
# convert masses to Earth masses and add extra stellar uncertainty using msigma defined above
# and put into allparam
for i in range(npl):
  allparam['M$_%s$' % alphabet[i]] = allparam['M$_{jup,}$$_%s$' % alphabet[i]]/MEOMJ * (1.+np.random.randn(len(allparam['M$_{jup,}$$_%s$' % alphabet[i]]))*msigma)

# make figure showing histogram of mass values in Earth masses
# (if more than 4 planets, put into two columns)
if npl > 4:
  fig, axes = plt.subplots(nrows=sum(divmod(npl, 2)), ncols=2, figsize=(6, 2*sum(divmod(npl,2))))
  for i in range(npl):
    allparam.hist('M$_%s$' % alphabet[i], ax=axes[i/2, i % 2], density=1)
  axes[sum(divmod(npl,2))-1,0].set_xlabel('M$_{Earth}$', fontsize=14)
  axes[npl/2-1,1].set_xlabel('M$_{Earth}$', fontsize=14)
  fig.savefig('analysis_dir/m_'+savestr+'.png')
else:
  fig, axes = plt.subplots(nrows=npl, ncols=1, figsize=(6, 3*npl))
  for i in range(npl):
    allparam.hist('M$_%s$' % alphabet[i], ax=axes[i], density=1)
  axes[npl-1].set_xlabel('M$_{Earth}$', fontsize=14)
  fig.savefig('analysis_dir/m_'+savestr+'.png')


# make figure showing histogram of radius values in Earth radii

REORS=0.009171 # convert solar radii to Earth radii
# add radius to allparam
for i in range(npl):
  allparam['R$_%s$' % alphabet[i]] = allparam['R$_p$/R$_s$$_%s$' % alphabet[i]]*allparam['R$_s$']/REORS

if npl > 4:
  fig2, axes2 = plt.subplots(nrows=sum(divmod(npl, 2)), ncols=2, figsize=(6, 2*sum(divmod(npl,2))))
  for i in range(npl):
    allparam.hist('R$_%s$' % alphabet[i], ax=axes2[i/2, i % 2], density=1)
  axes2[sum(divmod(npl,2))-1,0].set_xlabel('R$_{Earth}$', fontsize=14)
  axes2[npl/2-1,1].set_xlabel('R$_{Earth}$', fontsize=14)
  plt.tight_layout()
  fig2.savefig('analysis_dir/r_'+savestr+'.png')
else:
  fig2, axes2 = plt.subplots(nrows=npl, ncols=1, figsize=(6, 3*npl))
  for i in range(npl):
    allparam.hist('R$_%s$' % alphabet[i], ax=axes2[i], density=1)
  axes2[npl-1].set_xlabel('R$_{Earth}$', fontsize=14)
  plt.tight_layout()
  fig2.savefig('analysis_dir/r_'+savestr+'.png')



# make figure showing histogram of density values in g/cm^3
MEG=5.9721986*10**27 # earth in grams
RECM=6371.008*10**5  # earth in cm
# add density to allparam
for i in range(npl):
  allparam[r'$\rho$$_%s$' % alphabet[i]] = allparam['M$_%s$' % alphabet[i]]*MEG / (4./3.*np.pi * (RECM*allparam['R$_%s$' % alphabet[i]])**3. ) 

if npl>4:
  fig3, axes3 = plt.subplots(nrows=sum(divmod(npl, 2)), ncols=2, figsize=(6, 2*sum(divmod(npl,2))))
  for i in range(npl):
    allparam.hist(r'$\rho$$_%s$' % alphabet[i], ax=axes3[i/2, i % 2], density=1)
  axes3[sum(divmod(npl,2))-1,0].set_xlabel('g cm$^{-3}$', fontsize=14)
  axes3[npl/2-1,1].set_xlabel('g cm$^{-3}$', fontsize=14)
  plt.tight_layout()
  fig3.savefig('analysis_dir/rho_'+savestr+'.png')
else:
  fig3, axes3 = plt.subplots(nrows=npl, ncols=1, figsize=(6, 3*npl))
  for i in range(npl):
    allparam.hist(r'$\rho$$_%s$' % alphabet[i], ax=axes3[i], density=1)
  axes3[npl-1].set_xlabel('g cm$^{-3}$', fontsize=14)
  plt.tight_layout()
  fig3.savefig('analysis_dir/rho_'+savestr+'.png')

# add eccentricity to allparam
for i in range(npl):
  allparam['$e_%s$' % alphabet[i]] = allparam[r'$\sqrt{e}\cos \omega$$_%s$' % alphabet[i]]**2. + allparam[r'$\sqrt{e}\sin \omega$$_%s$' % alphabet[i]]**2.


# TODO: add mutual inclinations to allparam




# Lithwick et al. 2012 and others have shown that masses and eccentricities are often
# degenerate from Kepler TTVs. Other studies (e.g., Holman et al. 2010) have shown that
# masses of planets are often correlated in systems. As a result, it is useful to 
# investigate mass and eccentricity relationships directly in a smaller corner plot.

# gather the masses and ecentricities 
mecols=[]
for i in range(npl):
  mecols.append('M$_%s$' % alphabet[i])
for i in range(npl):
  mecols.append('$e_%s$' % alphabet[i])
mepars = allparam[mecols]

ranges=[(min(mepars.iloc[:,i]), max(mepars.iloc[:,i])) for i in range(mepars.iloc[0].size)]
ranges=[(i[0],i[0]+1) if i[0] == i[1] else i for i in ranges]

# make corner plot
figure3 = corner.corner(mepars, range=ranges)
figure3.savefig('analysis_dir/me_'+savestr+'.png')



# Make files that provide summaries of each parameter, including derived parameters.
# These use CREDIBLE INTERVALS based on percentiles that make no assumption about shape
# Note that these are only correct approximations of the posterior distribution
# if the Effective Sample Size is big enough to sample the tails of the 
# distribution. In particular, the 3-sigma limits are unlikely to be accurate
# if ESS < ~1000, though it may still be useful to see the range of outputs.

# "fits_1sigma" gives median and 16-84 percentile credible interval (NOT Guassian)
# "fits_3sigma" gives median and 0.13-99.86 percentile credible interval (NOT Guassian, possibly inaccurate)
# "fits_2sigmaUpperLimits" gives 95.4% upper limit on parameters 


# Find and print median and error bars of each parameter
sigma1o2=0.682689492137/2*100
sigma3o2=0.997300203936740/2*100
sigma2ul=0.954499736103642*100
percentiles=[50.-sigma3o2,50.-sigma1o2,50.,50.+sigma1o2,50.+sigma3o2,sigma2ul]

fullplist = allparam.columns 
#print(fullplist)

for p in fullplist:
  if (p[:6] != 'Planet' and p!= 'Chain#'):
    plistp = np.percentile(allparam[p], percentiles)
    with open('analysis_dir/fits_1sigma_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s\t%20.10f + %20.10f - %20.10f \n' %('{0: <16}'.format(p), plistp[2], plistp[3]-plistp[2], plistp[2]-plistp[1]))
    with open('analysis_dir/fits_3sigma_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s\t%20.10f + %20.10f - %20.10f \n' %('{0: <16}'.format(p), plistp[2], plistp[4]-plistp[2], plistp[2]-plistp[0]))
    with open('analysis_dir/fits_2sigmaUpperLimits_'+savestr+'.txt', 'a') as pctfile:
      pctfile.write('%s < %20.10f \n' %('{0: <16}'.format(p), plistp[5])) 




