import numpy as np
import matplotlib.pyplot as plt

# In [20]
import plot_utils as pread
from prospect.io.read_results import results_from

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
    

