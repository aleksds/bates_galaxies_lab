# In [1]
import time, sys, os
import h5py
import numpy as np
import matplotlib.pyplot as plt

# In [2]
import fsps
import sedpy
import prospect
import emcee
import sg_params as params

# In [3]
from prospect.models import model_setup

#paramfile = 'demo_mock_params1.py'
paramfile = 'sg_params.py'
#photfile ='demo_photometry1.dat'

# In [4]
clargs = {'param_file':paramfile}
run_params = model_setup.get_run_params(argv=paramfile, **clargs)
print(run_params)

# In [5]
# load sps model (default)
sps = model_setup.load_sps(**run_params)
print(sps)

# In [6]
# load noise model (none)
spec_noise, phot_noise = model_setup.load_gp(**run_params)

# In [7]
# demo model
model = model_setup.load_model(**run_params)

# In [8]
# demo data (generated from the script)
# note: the phottable=photfile line is an addiiton
#obs = model_setup.load_obs(phottable = photfile,**run_params)
obs = params.load_obs(**run_params)
#obs = model_setup.load_obs(**run_params)
#print('Mock S/N={}'.format(obs['mock_snr']))
#if run_params['add_noise']:
#    print('Noise realization added to mock photometry')
#else:
#    print('No noise added to mock photometry')

# In [9]
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

# In [10]
def lnprobfn(theta):
    """Given a parameter vector, a dictionary of observational data 
    and a model object, return the ln of the posterior. 
    This requires that an sps object (and if using spectra 
    and gaussian processes, a GP object) be instantiated.
    """

    print('lnprobfn loves pizza')
    # Calculate prior probability and exit if not within prior
    lnp_prior = model.prior_product(theta)
    if not np.isfinite(lnp_prior):
        print('oh shit prior')
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

# In [11]
from prospect.io import write_results

# In [12]
outroot = "{0}_{1}".format(run_params['outfile'], int(time.time()))
try:
    hfilename = outroot + '_mcmc.h5'
    hfile = h5py.File(hfilename, "a")
    print("Writing to file {}".format(hfilename))
    write_results.write_h5_header(hfile, run_params, model)
    write_results.write_obs_to_h5(hfile, obs)
except:
    hfile = None

# In [13]
print ('Free params:', model.free_params)
print ('Fixed params:', model.fixed_params)

# In [14]
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

# In [15]
from prospect import fitting
from scipy.optimize import least_squares

# In [16]
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
print('your mom', pinitial)
for i, pinit in enumerate(pinitial): #loop over initial guesses
    print(str(i), pinit)
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

print('i love pizza')

# In [17]
# initial parameters
theta = model.rectify_theta(initial_center)
# generate model
mspec, mphot, mextra = model.mean_model(theta, obs, sps=sps)

# establish bounds
xmin, xmax = wphot.min()*0.8, wphot.max()/0.8
temp = np.interp(np.linspace(xmin,xmax,10000), wspec*a, mspec)
ymin, ymax = temp.min()*0.8, temp.max()/0.8
ymax = 1e-7

fig = plt.figure()
# plot Data, models, and old models
plt.loglog(wspec * a, mspec_init, label='Old model spectrum',
       lw=0.7, color='gray', alpha=0.5)
plt.errorbar(wphot, mphot_init, label='Old model Photometry', 
         marker='s', markersize=10, alpha=0.6, ls='', lw=3, 
         markerfacecolor='none', markeredgecolor='gray', 
         markeredgewidth=3)
plt.loglog(wspec * a, mspec, label='Model spectrum', 
       lw=0.7, color='navy', alpha=0.7)
plt.errorbar(wphot, mphot, label='Model photometry', 
         marker='s', markersize=10, alpha=0.8, ls='', lw=3,
         markerfacecolor='none', markeredgecolor='blue', 
         markeredgewidth=3)
plt.errorbar(wphot, obs['maggies'], yerr=obs['maggies_unc'],
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
    plt.loglog(w, t, lw=3, color='gray', alpha=0.7)

# Prettify
plt.xlabel('Wavelength [A]')
plt.ylabel('Flux Density [maggies]')
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.legend(loc='best', fontsize=20)
plt.tight_layout()

plt.savefig('ho.png')

print('you love pizza')

# In [18] --> this is the step where the code breaks
postkwargs = {} #any keyord arguments to the lnpostfn would go here.

print('pizza sucks')

#fout = sys.stdout
#fnull = open(os.devnull, 'w')
#sys.stdout = fnull

tstart = time.time()  # time it
print('your mom loves pizza')
out = fitting.run_emcee_sampler(lnprobfn, initial_center, model, postkwargs=postkwargs, initial_prob=initial_prob, pool=None, hdf5=hfile, **run_params)

print('they love pizza')

esampler, burn_p0, burn_prob0 = out
edur = time.time() - tstart

sys.stdout = fout

print('done emcee in {0}s'.format(edur))

# In [19]
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

# In [20]
import plot_utils as pread
from prospect.io.read_results import results_from

# In [21]
# grab results, powell results, and our corresponding models
res, pr, mod = results_from("{}_mcmc.h5".format(outroot))

# In [22]
# To see how our MCMC samples look, we can examine a few traces
choice = np.random.choice
tracefig = pread.param_evol(res, figsize=(20,10), 
                            chains=choice(128, size=10, replace=False))

# In [23]
# Show samples in a triangle or corner plot
theta_truth = np.array([run_params[i] 
                        for i in ['mass','logzsol','tau','tage','dust2']])
theta_truth[0] = np.log10(theta_truth[0])
cornerfig = pread.subtriangle(res, start=0, thin=5, truths=theta_truth, fig=subplots(5,5,figsize=(27,27))[0])

# In [24]

# randomly chosen parameters from chain
randint = np.random.randint
nwalkers, niter = run_params['nwalkers'], run_params['niter']
theta = res['chain'][randint(nwalkers), randint(niter)]
# generate model
mspec, mphot, mextra = model.mean_model(theta, obs, sps=sps)

# establish bounds
xmin, xmax = wphot.min()*0.8, wphot.max()/0.8
temp = np.interp(np.linspace(xmin,xmax,10000), wspec * a, mspec)
ymin, ymax = temp.min()*0.8, temp.max()/0.8
figure(figsize=(16,8))

fig = plt.figure()

# plot data and model
plt.loglog(wspec * a, mspec, label='Model spectrum',
       lw=0.7, color='navy', alpha=0.7)
plt.errorbar(wphot, mphot, label='Model photometry',
         marker='s', markersize=10, alpha=0.8, ls='', lw=3, 
         markerfacecolor='none', markeredgecolor='blue', 
         markeredgewidth=3)
plt.errorbar(wphot, obs['maggies'], yerr=obs['maggies_unc'], 
         label='Observed photometry', ecolor='red', 
         marker='o', markersize=10, ls='', lw=3, alpha=0.8, 
         markerfacecolor='none', markeredgecolor='red', 
         markeredgewidth=3)

# plot transmission curves
for f in obs['filters']:
    w, t = f.wavelength.copy(), f.transmission.copy()
    while t.max() > 1:
        t /= 10.
    t = 0.1*(ymax-ymin)*t + ymin
    plt.loglog(w, t, lw=3, color='gray', alpha=0.7)

plt.xlabel('Wavelength [A]')
plt.ylabel('Flux Density [maggies]')
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.legend(loc='best', fontsize=20)
plt.tight_layout()

plt.savefig('ha.png')

