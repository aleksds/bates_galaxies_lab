import prospect.io.read_results as pread
import numpy as np
import os
import corner
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob


#h5name_arr = glob.glob('j*dtau_ssp*.h5')
h5name_arr = glob.glob('j1558_dust_ssp*.h5')
#h5name_arr = glob.glob('j*ssp*.h5')
#h5name_arr = glob.glob('j*dust_ssp*.h5')


def smpd(h5name):
    print(h5name)
    res, obs, mod = pread.results_from(h5name)#("demo_obj_<fitter>_<timestamp>_mcmc.h5")



    # Maximum posterior probability sample
    imax = np.argmax(res['lnprobability'])
    csz = res["chain"].shape
    if res["chain"].ndim > 2:
        # emcee
        i, j = np.unravel_index(imax, res['lnprobability'].shape)
        theta_max = res['chain'][i, j, :].copy()
        flatchain = res["chain"].reshape(csz[0] * csz[1], csz[2])
    else:
        # dynesty
        theta_max = res['chain'][imax, :].copy()
        flatchain = res["chain"]

    # 16th, 50th, and 84th percentiles of the posterior
    from prospect.utils.plotting import quantile
    post_pcts = [quantile(flatchain[:, i], percents=[16, 50, 84],
                                        weights=res.get("weights", None))
                          for i in range(mod.ndim)]

    # We need the correct sps object to generate models
    sps = pread.get_sps(res)

    # Choose the walker and iteration number,
    walker, iteration = 0, -1
    if res["chain"].ndim > 2:
        # if you used emcee for the inference
        theta = res['chain'][walker, iteration, :]
    else:
        # if you used dynesty
        theta = res['chain'][iteration, :]

    # Get the modeled spectra and photometry.
    # These have the same shape as the obs['spectrum'] and obs['maggies'] arrays.
    spec, phot, mfrac = mod.mean_model(theta, obs=res['obs'], sps=sps)
    # mfrac is the ratio of the surviving stellar mass to the formed mass (the ``"mass"`` parameter).

    # Plot the model SED
    filename = 'plot_demo_'+h5name+'.pdf'
    with PdfPages(filename) as pdf:

        fig = plt.figure()
        #res, obs, mod = pread.results_from('demo_galphot_1548731319_mcmc.h5')
        # Trace plots
        tfig = pread.traceplot(res)
        pdf.savefig()
        plt.close()
        # Corner figure of posterior PDFs
        cfig = pread.subcorner(res, logify=res['theta_labels'])#['mass', 'tau', 'tage'])
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        wave = [f.wave_effective for f in res['obs']['filters']]
        sedfig, sedax = plt.subplots()
        sedax.plot(wave, res['obs']['maggies'], '-o', label='Observations')
        sedax.plot(wave, phot, '-o', label='Model at {},{}'.format(walker, iteration))
        sedax.set_ylabel("Maggies")
        sedax.set_xlabel("wavelength")
        sedax.set_xscale('log')
        # plot filter transmission curves
        sedax.set_yscale('log')
        ymin = 5e-9
        ymax = 1e-7
        sedax.set_ylim([ymin, ymax])
        for f in obs['filters']:
            w, t = f.wavelength.copy(), f.transmission.copy()
            while t.max() > 1:
                t /= 10.
            t = 0.1*(ymax-ymin)*t + ymin
            sedax.plot(w, t, lw=3, color='gray', alpha=0.7)
        plt.legend(loc='upper left', prop={'size': 10})
        pdf.savefig()
        plt.close()

        # Plot residuals for this walker and iteration
        chifig, chiax = plt.subplots()
        chi = (res['obs']['maggies'] - phot) / res['obs']['maggies_unc']
        chiax.plot(wave, chi, 'o')
        chiax.set_ylabel("Chi")
        chiax.set_xlabel("wavelength")
        chiax.set_xscale('log')
        pdf.savefig()
        plt.close()

    
        os.system('open %s &' % filename)

for i in h5name_arr:
    print(i)
    smpd(i)