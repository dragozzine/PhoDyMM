#
#	phodymmplots.py
#
#	Makes likelihood vs walker position plots for PhoDyMM
#
#	Benjamin Proudfoot
#	05/18/20
#
#
#	Quick note:
#
# So I looked at demcmc_quick_analyze.py and I couldn't really understand the structure of 
# the allparam dataframe. I think it has each walker basically concatenated after the last
# walker, but I'm really unsure if my reading of the code is right. With a touch more direction
# I could easily make the follow code be easily inputted into demcmc_quick_analyze.py. 
#

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner

def unpack_df(allparam):
	# First unpack the allparam df into a 3d chain with dimensions
	# (# of generations) x (# of chains) x (# of parameters)
	data = allparam.values

	nparams = data.shape[1] - 2
	nwalkers = int(data[:,-1].max() + 1)
	ngens = int(data.shape[0] / nwalkers)

	chain = np.empty((ngens, nwalkers, nparams))
	for i in range(nwalkers):
		chain[:,i,:] = data[i*ngens:((i + 1)*ngens), 0:nparams]
	
	likelihood = np.empty((ngens, nwalkers))
	for i in range(nwalkers):
		likelihood[:,i] = data[i*ngens:((i+1)*ngens), -2]

	return chain, likelihood

def likelihood_plots(allparam):
	# Getting a useable chain from unpack_df
	chain, likelihood = unpack_df(allparam)

	# Also converting into a flattened chain with dimensions
	# (# of generations times # of chains) x (# of parameters)

	data = allparam.values
	nparams = chain.shape[2]
	nwalkers = chain.shape[1]
	flatchain = data[:,0:nparams]
	flatlikelihood = data[:,-2]
	chainnumbers = data[:,-1]
	names = list(allparam)
	
	# Likelihood plots??
	for i in range(nparams):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		for j in range(nwalkers):
			plt.hist(chain[:,j,i].flatten(), bins = 20, histtype = "step", 
				color = "black",
				alpha = 0.4, density = True)
		plt.hist(chain[:,:,i].flatten(), bins = 20, histtype = "step", color = "black", density = True)
		plt.subplot(223)
		plt.scatter(flatchain[:,i].flatten(), flatlikelihood.flatten(), 
			    c = chainnumbers.flatten()/chainnumbers.max(), cmap = "nipy_spectral", 
			    edgecolors = "none")
		plt.xlabel(names[i])
		plt.ylabel("Log(L)")
		plt.subplot(224)
		plt.hist(flatlikelihood[~np.isnan(flatlikelihood)], bins = 20,
			 orientation = "horizontal",histtype = "step", color = "black")
		plt.savefig("./testing/phodymm/" + names[i] + "_likelihood.png")
		plt.close("all")

# Autocorrelation functions

def auto_window(taus, c):
        m = np.arange(len(taus)) < c * taus
        if np.any(m):
                return np.argmin(m)
        return len(taus) - 1

def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i

def autocorr_func_1d(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf

def autocorr_new(y, c = 5.0):
        f = np.zeros(y.shape[1])
        for yy in y:
                f += autocorr_func_1d(yy)
        f /= len(y)
        taus = 2.0 * np.cumsum(f) - 1.0
        window = auto_window(taus, c)
        return taus[window]

def autocorrelation(allparam, thin):
        # Getting chain
        chain, likelihood = unpack_df(allparam)

	# Pulling out needed values from the chain
        nwalkers = chain.shape[1]
        ndims = chain.shape[2]
        nsteps = chain.shape[0]

        # Converting parameters df to array of names
        names = list(allparam)

        # Calculating values to calculate tau for
        # This chould be changed eventually
        N = np.exp(np.linspace(np.log(1), np.log(nsteps), 10)).astype(int)

        # Setting up array for tau estimates
        tau = np.empty( (len(N), ndims) )

        # Calculating tau for each value in N for each parameter
        for i in range(ndims):
                thischain = chain[:,:,i].T
                for j, n in enumerate(N):
                        tau[j,i] = autocorr_new(thischain[:, :n])

        # Setting up to plot curves for tau in a grid
        x = 3
        y = ndims
        nrows = 0
        ncols = 3
        while x <= y:
                y = y - x
                nrows += 1

        # Plotting
        fig, ax = plt.subplots(nrows = nrows, ncols = ncols, sharey=True,
                               gridspec_kw={'wspace': 0},
                               figsize = (6.4*(ncols-1),4.8*(nrows)),
                               squeeze = False)
        fig.suptitle("Autocorrelation estimates")
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.xlabel("number of samples, $N$")
        plt.ylabel(r"$\tau$ estimates")
        for i in range(nrows):
                for j in range(ncols):
                        dim = i*ncols + j
                        taus = ax[i,j].loglog(N, tau[:,dim], "o-", label="new")
                        line = ax[i,j].plot(N, N / 50.0, "--k", label=r"$\tau = N/50$")
                        ax[i,j].text(N[0], 100, names[dim])
        fig.savefig("autocorr.png")
        # Outputting the emcee calculated autocorrelation time as an additional check
#        print(
#                "Mean autocorrelation time: {0:.3f} steps".format(
#                        np.mean(sampler.get_autocorr_time())
#                )
#        )
#        print("Estimated autocorrelation time for each parameter:")
#        print(sampler.get_autocorr_time())

