import numpy as np
import matplotlib.pyplot as plt
import spot_utils as pread
from prospect.io.read_results import results_from
from prospect.models import model_setup


paramfile = 'ad_params.py'
clargs = {'param_file':paramfile}
run_params = model_setup.get_run_params(argv=paramfile, **clargs)
print(run_params)

obs = model_setup.load_obs(**run_params)
sps = model_setup.load_sps(**run_params)
model = model_setup.load_model(**run_params)

wspec = sps.csp.wavelengths # *restframe* spectral wavelengths
a = 1.0 + model.params.get('zred', 0.6) # cosmological redshifting
wphot = np.array([f.wave_effective for f in obs['filters']])

# In [21]
# grab results, powell results, and our corresponding models
#res, pr, mod = results_from("{}_mcmc.h5".format(outroot))
res, pr, mod = results_from("demo_galphot_1508433060_mcmc.h5")

# In [22]
# To see how our MCMC samples look, we can examine a few traces
choice = np.random.choice
tracefig = pread.param_evol(res, figsize=(20,10), 
                            chains=choice(128, size=10, replace=False))

# In [23]
# Show samples in a triangle or corner plot
#theta_truth = np.array([pr[i] 
#                        for i in ['mass','logzsol','tau','tage','dust2']])
theta_truth = pr[0]['x']
theta_truth[0] = np.log10(theta_truth[0])
#cornerfig = pread.subtriangle(res, start=0, thin=5, truths=theta_truth, showpars='mass')#, fig=plt.subplots(5,5,figsize=(27,27))[0])
truths = theta_truth

import corner as triangle
sample_results = res
#parnames = np.array(sample_results['theta_labels'])
parnames = res['theta_labels']
flatchain = sample_results['chain']
flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])

fig = plt.figure()

# logify mass
#if 'mass' in parnames:
#    midx = [l=='mass' for l in parnames]
#    flatchain[:,midx] = np.log10(flatchain[:,midx])
#    parnames[midx] = 'logmass'
#flatchain = np.ma.masked_invalid(flatchain)
# restrict to parameters you want to show
#if showpars is not None:
#    ind_show = np.array([p in showpars for p in parnames], dtype=bool)
#    flatchain = flatchain[:, ind_show]
#    #truths = truths[ind_show]
#    parnames = parnames[ind_show]
trim_outliers=None
if trim_outliers is not None:
    trim_outliers = len(parnames) * [trim_outliers]
try:
    fig = triangle.corner(flatchain, labels=parnames, truths=truths[0],  verbose=False, quantiles=[0.16, 0.5, 0.84])#, range=trim_outliers)#, **kwargs)
except:
    fig = triangle.corner(flatchain, labels=parnames, truths=truths,  verbose=False,
                              quantiles=[0.16, 0.5, 0.84])#, range=trim_outliers)#, **kwargs)

plt.savefig('yo.png')
    
# In [24]

# randomly chosen parameters from chain
run_params = res['run_params']
randint = np.random.randint
nwalkers, niter = run_params['nwalkers'], run_params['niter']
theta = res['chain'][randint(nwalkers), randint(niter)]
# generate model
model = mod
mspec, mphot, mextra = model.mean_model(theta, obs, sps=sps)

# establish bounds
xmin, xmax = wphot.min()*0.8, wphot.max()/0.8
temp = np.interp(np.linspace(xmin,xmax,10000), wspec * a, mspec)
ymin, ymax = temp.min()*0.8, temp.max()/0.8
#figure(figsize=(16,8))

fig = plt.figure(figsize=(16,8))

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
