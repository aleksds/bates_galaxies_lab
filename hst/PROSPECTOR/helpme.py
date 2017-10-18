from __future__ import print_function, absolute_import
import time, sys, os
import h5py
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
#matplotlib inline

# re-defining plotting defaults
from matplotlib.font_manager import FontProperties
from matplotlib import gridspec
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'xtick.color': 'k'})
rcParams.update({'ytick.color': 'k'})
rcParams.update({'font.size': 30})
import sys
sys.path.append('/Users/sgottlie/github/prospector')
sys.path.append('/Users/sgottlie/github/sedpy')
sys.path.append('/Users/sgottlie/github/python-fsps')
sys.path.append('/Users/sgottlie/github/corner')
import fsps
import sedpy
import prospect
from prospect.models import model_setup

# maybe if i import the stuff from the beginning of corner, the dumb part about quantile will work

hist2d_kwargs = {}
smooth = None

import logging
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator, NullLocator
from matplotlib.colors import LinearSegmentedColormap, colorConverter
from matplotlib.ticker import ScalarFormatter
#okay that is in, let's try this bullshit

## HERE IS A LIST OF THINGS THAT I REFUSE TO HARD CODE INTO THE PROGRAM!!!
#paramfile = '/Users/sgottlie/github/prospector/demo/demo_mock_params1.py'
paramfile = '/demo_mock_params1.py'
dmp1 = 'demo_mock_params1.py'
photfile ='/demo_photometry1.dat'
### OKAY ALL DONE

clargs = {'param_file':paramfile}
run_params = model_setup.get_run_params(argv=dmp1, **clargs)
print(run_params)


# load sps model (default)
sps = model_setup.load_sps(**run_params)

# load noise model (none)
spec_noise, phot_noise = model_setup.load_gp(**run_params)

# demo model
model = model_setup.load_model(**run_params)

# demo data (generated from the script)
obs = model_setup.load_obs(phottable = photfile,**run_params)
print('Mock S/N={}'.format(obs['mock_snr']))
if run_params['add_noise']:
    print('Noise realization added to mock photometry')
else:
    print('No noise added to mock photometry')


    from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

def quantile(x, q, weights=None):
    """
    Compute sample quantiles with support for weighted samples.

    Note
    ----
    When ``weights`` is ``None``, this method simply calls numpy's percentile
    function with the values of ``q`` multiplied by 100.

    Parameters
    ----------
    x : array_like[nsamples,]
       The samples.

    q : array_like[nquantiles,]
       The list of quantiles to compute. These should all be in the range
       ``[0, 1]``.

    weights : Optional[array_like[nsamples,]]
        An optional weight corresponding to each sample. These

    Returns
    -------
    quantiles : array_like[nquantiles,]
        The sample quantiles computed at ``q``.

    Raises
    ------
    ValueError
        For invalid quantiles; ``q`` not in ``[0, 1]`` or dimension mismatch
        between ``x`` and ``weights``.

    """
    x = np.atleast_1d(x)
    q = np.atleast_1d(q)

    if np.any(q < 0.0) or np.any(q > 1.0):
        raise ValueError("Quantiles must be between 0 and 1")

    if weights is None:
        return np.percentile(x, list(100.0 * q))
    else:
        weights = np.atleast_1d(weights)
        if len(x) != len(weights):
            raise ValueError("Dimension mismatch: len(weights) != len(x)")
        idx = np.argsort(x)
        sw = weights[idx]
        cdf = np.cumsum(sw)[:-1]
        cdf /= cdf[-1]
        cdf = np.append(0, cdf)
        return np.interp(q, cdf, x[idx]).tolist()
    
def lnprobfn(theta):
    """Given a parameter vector, a dictionary of observational data 
    and a model object, return the ln of the posterior. 
    This requires that an sps object (and if using spectra 
    and gaussian processes, a GP object) be instantiated.
    """

    # Calculate prior probability and exit if not within prior
    lnp_prior = model.prior_product(theta)
    if not np.isfinite(lnp_prior):
        return -np.infty
        
    # Generate mean model
    t1 = time.time()
    try:
        spec, phot, x = model.mean_model(theta, obs, sps=sps)
    except(ValueError):
        return -np.infty
    d1 = time.time() - t1
    vectors = {}  # This would be used for noise model weight functions

    # Calculate likelihoods
    t2 = time.time()
    lnp_spec = lnlike_spec(spec, obs=obs, spec_noise=spec_noise)
    lnp_phot = lnlike_phot(phot, obs=obs, phot_noise=phot_noise)
    d2 = time.time() - t2
    if verbose:
        write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

    return lnp_prior + lnp_phot + lnp_spec
#RUNNING PROSPECTOR
# Outputs

from prospect.io import write_results

outroot = "{0}_{1}".format(run_params['outfile'], int(time.time()))
try:
    hfilename = outroot + '_mcmc.h5'
    hfile = h5py.File(hfilename, "a")
    print("Writing to file {}".format(hfilename))
    write_results.write_h5_header(hfile, run_params, model)
    write_results.write_obs_to_h5(hfile, obs)
except:
    hfile = None

print ('Free params:', model.free_params)
print ('Fixed params:', model.fixed_params)

# SED Preview

wspec = sps.csp.wavelengths # *restframe* spectral wavelengths
a = 1.0 + model.params.get('zred', 0.6) # cosmological redshifting
# photometric effective wavelengths
wphot = np.array([f.wave_effective for f in obs['filters']])
# initial parameters based on `init` values of model_params
initial_theta = model.rectify_theta(model.initial_theta)
# generate model
out_init = model.mean_model(initial_theta, obs, sps=sps) 
mspec_init, mphot_init, mextra_init = out_init

# establish bounds
xmin, xmax = np.min(wphot)*0.8, np.max(wphot)/0.8
temp = np.interp(np.linspace(xmin,xmax,10000), wspec * a, mspec_init)
ymin, ymax = temp.min()*0.8, temp.max()/0.8
figure(figsize=(16,8))
ymax = 1e-7

# plot model + data
fig = plt.figure()

plt.loglog(wspec * a, mspec_init, label='Model spectrum', 
       lw=0.7, color='navy', alpha=0.7)
plt.errorbar(wphot, mphot_init, label='Model photometry', 
         marker='s',markersize=10, alpha=0.8, ls='', lw=3,
         markerfacecolor='none', markeredgecolor='blue', 
         markeredgewidth=3)
plt.errorbar(wphot, obs['maggies'], yerr=obs['maggies_unc'], 
         label='Observed photometry',
         marker='o', markersize=10, alpha=0.8, ls='', lw=3,
         ecolor='red', markerfacecolor='none', markeredgecolor='red', 
         markeredgewidth=3)

# plot Filters
for f in obs['filters']:
    print(str(f))
    w, t = f.wavelength.copy(), f.transmission.copy()
    while t.max() > 1:
        t /= 10.
    t = 0.1*(ymax-ymin)*t + ymin
    plt.loglog(w, t, lw=3, color='gray', alpha=0.7)

# prettify
plt.xlabel('Wavelength [A]')
plt.ylabel('Flux Density [maggies]')
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.legend(loc='best', fontsize=20)
plt.tight_layout()

plt.savefig('hi.png')

# Minimization Step
# We can attempt to initialize our model reasonably close to the data by using some numerical minimization routines.

from prospect import fitting
from scipy.optimize import least_squares

from prospect.likelihood import chi_spec, chi_phot
def chivecfn(theta):
    """A version of lnprobfn that returns the simple uncertainty 
    normalized residual instead of the log-posterior, for use with 
    least-squares optimization methods like Levenburg-Marquardt.
    """
    lnp_prior = model.prior_product(theta)
    if not np.isfinite(lnp_prior):
        return -np.infty

    # Generate mean model
    t1 = time.time()
    try:
        spec, phot, x = model.mean_model(theta, obs, sps=sps)
    except(ValueError):
        return -np.infty
    d1 = time.time() - t1

    chispec = chi_spec(spec, obs)
    chiphot = chi_phot(phot, obs)
    return np.concatenate([chispec, chiphot])

verbose = False # don't output function calls

# start minimization
min_method = 'levenburg-marquardt'
nmin = 5 # We'll start from 5 places, 4 of which are drawn from the prior
ts = time.time()
pinitial = fitting.minimizer_ball(model.initial_theta.copy(), nmin, model)
guesses = []
for i, pinit in enumerate(pinitial): #loop over initial guesses
    print(str(i))
    res = least_squares(chivecfn, pinit, method='dogbox', x_scale='jac',
                        xtol=1e-16, ftol=1e-16)
    guesses.append(res)

# Calculate chi-square of the results, and choose the best one
chisq = [np.sum(r.fun**2) for r in guesses]
best = np.argmin(chisq)
initial_center = fitting.reinitialize(guesses[best].x, model,
                        edge_trunc=run_params.get('edge_trunc', 0.1))
initial_prob = None
pdur = time.time() - ts

# output results
print('done {0} in {1}s'.format(min_method, pdur))
print('best {0} guess: {1}'.format(min_method, initial_center))
print('best {0} chi-sq: {1}'.format(min_method, chisq[best]))

'''Note that creating a new model with FSPS is somewhat time-intensive,
but once the relevant model(s) have been loaded they are subsequently
stored in cache so similar models can be generated much more quickly.
Now let's see how our model looks.'''

# initial parameters
theta = model.rectify_theta(initial_center)
# generate model
mspec, mphot, mextra = model.mean_model(theta, obs, sps=sps)

# establish bounds
xmin, xmax = wphot.min()*0.8, wphot.max()/0.8
temp = np.interp(np.linspace(xmin,xmax,10000), wspec*a, mspec)
ymin, ymax = temp.min()*0.8, temp.max()/0.8
figure(figsize=(16,8))
ymax = 1e-7

# plot Data, models, and old models
loglog(wspec * a, mspec_init, label='Old model spectrum',
       lw=0.7, color='gray', alpha=0.5)
errorbar(wphot, mphot_init, label='Old model Photometry', 
         marker='s', markersize=10, alpha=0.6, ls='', lw=3, 
         markerfacecolor='none', markeredgecolor='gray', 
         markeredgewidth=3)
loglog(wspec * a, mspec, label='Model spectrum', 
       lw=0.7, color='navy', alpha=0.7)
errorbar(wphot, mphot, label='Model photometry', 
         marker='s', markersize=10, alpha=0.8, ls='', lw=3,
         markerfacecolor='none', markeredgecolor='blue', 
         markeredgewidth=3)
errorbar(wphot, obs['maggies'], yerr=obs['maggies_unc'],
         label='Observed photometry', 
         marker='o', markersize=10, alpha=0.8, ls='', lw=3, 
         ecolor='red', markerfacecolor='none', markeredgecolor='red', 
         markeredgewidth=3)

# plot filter transmission curves
for f in obs['filters']:
    w, t = f.wavelength.copy(), f.transmission.copy()
    while t.max() > 1:
        t /= 10.
    t = 0.1*(ymax-ymin)*t + ymin
    loglog(w, t, lw=3, color='gray', alpha=0.7)

# Prettify
xlabel('Wavelength [A]')
ylabel('Flux Density [maggies]')
xlim([xmin, xmax])
ylim([ymin, ymax])
legend(loc='best', fontsize=20)
tight_layout()

# Sampling the Posterior
'''Now that we're somewhat burned in, we can begin sampling from the
posterior using Markov Chain Monte Carlo (MCMC). Prospector by default
uses emcee, and will try to parallelize the process over multiple
cores when available through MPI and mpi4py. In this interactive
notebook though we will assume single-threaded operation. Let's go
ahead and start sampling!'''

postkwargs = {} #any keyord arguments to the lnpostfn would go here.

fout = sys.stdout
fnull = open(os.devnull, 'w')
sys.stdout = fnull

tstart = time.time()  # time it
out = fitting.run_emcee_sampler(lnprobfn, initial_center, model,
                                postkwargs=postkwargs, 
                                initial_prob=initial_prob,
                                pool=None, hdf5=hfile, **run_params)
esampler, burn_p0, burn_prob0 = out
edur = time.time() - tstart

sys.stdout = fout

print('done emcee in {0}s'.format(edur))

'''Now that everything's all set, let's save our results to disk. These
will be written to 2 or 3 files beginning with the value of outroot.'''

write_results.write_pickles(run_params, model, obs, esampler, guesses,
                            outroot=outroot, toptimize=pdur, tsample=edur,
                            sampling_initial_center=initial_center,
                            post_burnin_center=burn_p0,
                            post_burnin_prob=burn_prob0)

if hfile is None:
    hfile = hfilename
write_results.write_hdf5(hfile, run_params, model, obs, esampler, 
                         guesses,
                         toptimize=pdur, tsample=edur,
                         sampling_initial_center=initial_center,
                         post_burnin_center=burn_p0,
                         post_burnin_prob=burn_prob0)

print('Finished')

# Visualizing the Results
'''There are a few basic plotting tools available to do a quick check on
the results available in prospect.io.read_results and prospect.utils.plotting.
We'll hack a few of these together in plot_utils here in the demo folder to
make them a bit more amenable to plotting in this notebook.'''

import plot_utils as pread
from prospect.io.read_results import results_from

# grab results, powell results, and our corresponding models
res, pr, mod = results_from("{}_mcmc.h5".format(outroot))

# To see how our MCMC samples look, we can examine a few traces
choice = np.random.choice
tracefig = pread.param_evol(res, figsize=(20,10), 
                            chains=choice(128, size=10, replace=False))

# Show samples in a triangle or corner plot
theta_truth = np.array([run_params[i] 
                        for i in ['mass','logzsol','tau','tage','dust2']])
theta_truth[0] = np.log10(theta_truth[0])
parnames = np.array(res['theta_labels'])
flatchain = res['chain'][:, 0::5, :]
flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1], flatchain.shape[2])
# logify mass
if 'mass' in parnames:
    midx = [l=='mass' for l in parnames]
    flatchain[:,midx] = np.log10(flatchain[:,midx])
    parnames[midx] = 'logmass'
flatchain = np.ma.masked_invalid(flatchain)
range =[[x.min(), x.max()] for x in flatchain]
weights = None

#import triangle
theta_truth = np.array([run_params[i] 
                        for i in ['mass','logzsol','tau','tage','dust2']])
theta_truth[0] = np.log10(theta_truth[0])
# I AM GOING TO PUT SOME STUFF IN HERE THAT IS FROM PLOT_UTILS.PY
# cornerfig = pread.subtriangle(res, start=0, thin=5, truths=theta_truth, fig=subplots(5,5,figsize=(27,27))[0])

import corner as triangle

sample_results = res
start = 0
thin = 5
truths = theta_truth
fig = subplots(5,5,figsize=(27,27))[0]
parnames = np.array(sample_results['theta_labels'])
### have a bunch of stuff that is needed in plot utils
showpars = None
trim_outliers = None

flatchain = sample_results['chain'][:, start::thin, :]
flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],flatchain.shape[2])

###
if 'mass' in parnames:
    midx = [l=='mass' for l in parnames]
    flatchain[:,midx] = np.log10(flatchain[:,midx])
    parnames[midx] = 'logmass'
flatchain = np.ma.masked_invalid(flatchain)
# restrict to parameters you want to show
if showpars is not None:
    ind_show = np.array([p in showpars for p in parnames], dtype=bool)
    flatchain = flatchain[:, ind_show]
    #truths = truths[ind_show]
    parnames = parnames[ind_show]
if trim_outliers is not None:
    trim_outliers = len(parnames) * [trim_outliers]
#fig = triangle.corner(flatchain, labels=parnames, truths=truths,  verbose=False, quantiles=[0.16, 0.5, 0.84], range=trim_outliers, **kwargs)
# Lets recreate what is happening above
xs = flatchain
bins = 20
weights = None
range = None
color = 'k'
smooth1d = None
labels = parnames
label_kwargs = None
show_titles = None
title_fmt = ".2f"
truth_color = "#4682b4"
scale_hist = False
quantiles = [0.16, 0.5, 0.84]
verbose=False
range=trim_outliers
max_n_ticks = 5
top_ticks=False
use_math_text=False
reverse = False
hist_kwargs=None
title_kwargs = dict()
label_kwargs = dict()
#### here is triangle from corner

if quantiles is None:
    quantiles = []
if title_kwargs is None:
    title_kwargs = dict()
if label_kwargs is None:
    label_kwargs = dict()

# Try filling in labels from pandas.DataFrame columns.
if labels is None:
    try:
        labels = xs.columns
    except AttributeError:
        pass

# Deal with 1D sample lists.
xs = np.atleast_1d(xs)
if len(xs.shape) == 1:
    xs = np.atleast_2d(xs)
else:
    assert len(xs.shape) == 2, "The input sample array must be 1- or 2-D."
    xs = xs.T
assert xs.shape[0] <= xs.shape[1], "I don't believe that you want more " \
                                   "dimensions than samples!"

# Parse the weight array.
if weights is not None:
    weights = np.asarray(weights)
    if weights.ndim != 1:
        raise ValueError("Weights must be 1-D")
    if xs.shape[1] != weights.shape[0]:
        raise ValueError("Lengths of weights must match number of samples")

# Parse the parameter ranges.
if range is None:
#    if "extents" in hist2d_kwargs:
#        logging.warn("Deprecated keyword argument 'extents'. "
#                     "Use 'range' instead.")
#        range = hist2d_kwargs.pop("extents")
#    else:
    range = [[x.min(), x.max()] for x in xs]
    # Check for parameters that never change.
    m = np.array([e[0] == e[1] for e in range], dtype=bool)
    if np.any(m):
        raise ValueError(("It looks like the parameter(s) in "
                          "column(s) {0} have no dynamic range. "
                          "Please provide a `range` argument.")
                         .format(", ".join(map(
                             "{0}".format, np.arange(len(m))[m]))))

else:
    # If any of the extents are percentiles, convert them to ranges.
    # Also make sure it's a normal list.
    range = list(range)
    for i, _ in enumerate(range):
        try:
            emin, emax = range[i]
        except TypeError:
            q = [0.5 - 0.5*range[i], 0.5 + 0.5*range[i]]
            range[i] = quantile(xs[i], q, weights=weights)

if len(range) != xs.shape[0]:
    raise ValueError("Dimension mismatch between samples and range")

# Parse the bin specifications.
try:
    bins = [int(bins) for _ in range]
except TypeError:
    if len(bins) != len(range):
        raise ValueError("Dimension mismatch between bins and range")

# Some magic numbers for pretty axis layout.
K = len(xs)
factor = 2.0           # size of one side of one panel
if reverse:
    lbdim = 0.2 * factor   # size of left/bottom margin
    trdim = 0.5 * factor   # size of top/right margin
else:
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.2 * factor   # size of top/right margin
whspace = 0.05         # w/hspace size
plotdim = factor * K + factor * (K - 1.) * whspace
dim = lbdim + plotdim + trdim

# Create a new figure if one wasn't provided.
if fig is None:
    fig, axes = pl.subplots(K, K, figsize=(dim, dim))
else:
    try:
        axes = np.array(fig.axes).reshape((K, K))
    except:
        raise ValueError("Provided figure has {0} axes, but data has "
                         "dimensions K={1}".format(len(fig.axes), K))

# Format the figure.
lb = lbdim / dim
tr = (lbdim + plotdim) / dim
fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                    wspace=whspace, hspace=whspace)

# Set up the default histogram keywords.
if hist_kwargs is None:
    hist_kwargs = dict()
hist_kwargs["color"] = hist_kwargs.get("color", color)
if smooth1d is None:
    hist_kwargs["histtype"] = hist_kwargs.get("histtype", "step")

for i, x in enumerate(xs):
    # Deal with masked arrays.
    if hasattr(x, "compressed"):
        x = x.compressed()

    if np.shape(xs)[0] == 1:
        ax = axes
    else:
        if reverse:
            ax = axes[K-i-1, K-i-1]
        else:
            ax = axes[i, i]
    # Plot the histograms.
    if smooth1d is None:
        n, _, _ = ax.hist(x, bins=bins[i], weights=weights,
                          range=np.sort(range[i]), **hist_kwargs)
        #                  range = range[i], **hist_kwargs)
    else:
        if gaussian_filter is None:
            raise ImportError("Please install scipy for smoothing")
        n, b = np.histogram(x, bins=bins[i], weights=weights,
                            range=np.sort(range[i]))
        n = gaussian_filter(n, smooth1d)
        x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
        y0 = np.array(list(zip(n, n))).flatten()
        ax.plot(x0, y0, **hist_kwargs)

    if truths is not None and truths[i] is not None:
        ax.axvline(truths[i], color=truth_color)

    # Plot quantiles if wanted.
    if len(quantiles) > 0:
        qvalues = quantile(x, quantiles, weights=weights)
        for q in qvalues:
            ax.axvline(q, ls="dashed", color=color)

        if verbose:
            print("Quantiles:")
            print([item for item in zip(quantiles, qvalues)])

    if show_titles:
        title = None
        if title_fmt is not None:
            # Compute the quantiles for the title. This might redo
            # unneeded computation but who cares.
            q_16, q_50, q_84 = quantile(x, [0.16, 0.5, 0.84],
                                        weights=weights)
            q_m, q_p = q_50-q_16, q_84-q_50

            # Format the quantile display.
            fmt = "{{0:{0}}}".format(title_fmt).format
            title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
            title = title.format(fmt(q_50), fmt(q_m), fmt(q_p))

            # Add in the column name if it's given.
            if labels is not None:
                title = "{0} = {1}".format(labels[i], title)

        elif labels is not None:
            title = "{0}".format(labels[i])

        if title is not None:
            if reverse:
                ax.set_xlabel(title, **title_kwargs)
            else:
                ax.set_title(title, **title_kwargs)

    # Set up the axes.
    ax.set_xlim(range[i])
    if scale_hist:
        maxn = np.max(n)
        ax.set_ylim(-0.1 * maxn, 1.1 * maxn)
    else:
        ax.set_ylim(0, 1.1 * np.max(n))
    ax.set_yticklabels([])
    if max_n_ticks == 0:
        ax.xaxis.set_major_locator(NullLocator())
        ax.yaxis.set_major_locator(NullLocator())
    else:
        ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))
        ax.yaxis.set_major_locator(NullLocator())

    if i < K - 1:
        if top_ticks:
            ax.xaxis.set_ticks_position("top")
            [l.set_rotation(45) for l in ax.get_xticklabels()]
        else:
            ax.set_xticklabels([])
    else:
        if reverse:
            ax.xaxis.tick_top()
        [l.set_rotation(45) for l in ax.get_xticklabels()]
        if labels is not None:
            if reverse:
                ax.set_title(labels[i], y=1.25, **label_kwargs)
            else:
                ax.set_xlabel(labels[i], **label_kwargs)

        # use MathText for axes ticks
        ax.xaxis.set_major_formatter(
            ScalarFormatter(useMathText=use_math_text))

    for j, y in enumerate(xs):
        if np.shape(xs)[0] == 1:
            ax = axes
        else:
            if reverse:
                ax = axes[K-i-1, K-j-1]
            else:
                ax = axes[i, j]
        if j > i:
            ax.set_frame_on(False)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        elif j == i:
            continue

        # Deal with masked arrays.
        if hasattr(y, "compressed"):
            y = y.compressed()
        '''
        h0tp00ps = ax[i]
        hist2d(y, x, ax=h0tp00ps, range=[range[j], range[i]], weights=weights, color=color, smooth=smooth, bins=[bins[j], bins[i]], **hist2d_kwargs)

        if truths is not None:
            if truths[i] is not None and truths[j] is not None:
                ax.plot(truths[j], truths[i], "s", color=truth_color)
            if truths[j] is not None:
                ax.axvline(truths[j], color=truth_color)
            if truths[i] is not None:
                ax.axhline(truths[i], color=truth_color)

        if max_n_ticks == 0:
            ax.xaxis.set_major_locator(NullLocator())
            ax.yaxis.set_major_locator(NullLocator())
        else:
            ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks,
                                                   prune="lower"))
            ax.yaxis.set_major_locator(MaxNLocator(max_n_ticks,
                                                   prune="lower"))

        if i < K - 1:
            ax.set_xticklabels([])
        else:
            if reverse:
                ax.xaxis.tick_top()
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            if labels is not None:
                ax.set_xlabel(labels[j], **label_kwargs)
                if reverse:
                    ax.xaxis.set_label_coords(0.5, 1.4)
                else:
                    ax.xaxis.set_label_coords(0.5, -0.3)

            # use MathText for axes ticks
            ax.xaxis.set_major_formatter(
                ScalarFormatter(useMathText=use_math_text))

        if j > 0:
            ax.set_yticklabels([])
        else:
            if reverse:
                ax.yaxis.tick_right()
            [l.set_rotation(45) for l in ax.get_yticklabels()]
            if labels is not None:
                if reverse:
                    ax.set_ylabel(labels[i], rotation=-90, **label_kwargs)
                    ax.yaxis.set_label_coords(1.3, 0.5)
                else:
                    ax.set_ylabel(labels[i], **label_kwargs)
                    ax.yaxis.set_label_coords(-0.3, 0.5)

            # use MathText for axes ticks
            ax.yaxis.set_major_formatter(
                ScalarFormatter(useMathText=use_math_text))
'''
#return fig

### all done with triangle from corner
# This is the end of plot_utils / subtriangle:
outname = 'help'
if outname is not None:
    fig.savefig('{0}.triangle.png'.format(outname))
    #pl.close(fig)
#else:
#    return fig

## hist2d(y, x, range=[range[j], range[i]], weights=weights, color=color, smooth=smooth, bins=20, **hist2d_kwargs)

# HERE IS THE END OF THE INTERACTIVE DEMO OMG

# plot transmission curves
fig = plt.figure()
for f in obs['filters']:
    w, t = f.wavelength.copy(), f.transmission.copy()
    while t.max() > 1:
        t /= 10.
    t = 0.1*(ymax-ymin)*t + ymin
    plt.loglog(w, t, lw=3, color='gray', alpha=0.7)


for f in obs['filters']:
    w,t = f.wavelength.copy(),f.transmission.copy()
    while t.max() > 1:
        t /= 10.
    t = 0.1*(ymax-ymin)*t +ymin
    plt.loglog(w,t,lw=3,color='gray',alpha=0.7)
plt.savefig('somewaves.png')
