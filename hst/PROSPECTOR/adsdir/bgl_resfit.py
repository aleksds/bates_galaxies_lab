#!/usr/bin/env python
"""
Fit the SED using Prospector.

run bgl_resfit.py --nproc=4 --priors='delayed-tau' --galaxy='J0901' --sedfit  --verbose
run bgl_resfit.py --nproc=4 --priors='delayed-tau' --galaxy='J0901' --qaplots --verbose

"""
import os, time, argparse, pdb
import numpy as np
import multiprocessing
import prospect

ang2micron = 1e-4 # Angstrom --> micron
maggies2mJy = 10**(0.4*16.4) # maggies --> mJy

def _niceparnames(parnames):
    """Replace parameter names with nice names."""

    old = list(['tau',
           'tage',
           'mass',
           'logmass',
           'logzsol',
           'dust2'])
    new = list([r'$\tau$ (Gyr)',
           'Age (Gyr)',
           r'$M / M_{\odot}$',
           r'$\log_{10}\,(M / M_{\odot})$',
           r'$\log_{10}\, (Z / Z_{\odot})$',
           r'$\tau_{diffuse}$'])

    niceparnames = list(parnames).copy()
    for oo, nn in zip( old, new ):
        this = np.where(np.in1d(parnames, oo))[0]
        if len(this) > 0:
            niceparnames[this[0]] = nn
            
    return np.array(niceparnames)

def _galaxyphot(obs):
    """Get the galaxy photometry and inverse variances (converted to mJy) and filter
    effective wavelengths (converted to microns).

    """
    weff = np.array([f.wave_effective for f in obs['filters']]) * ang2micron
    fwhm = np.array([f.effective_width for f in obs['filters']]) * ang2micron

    if False:
        galphot = obs['maggies'] * maggies2mJy
        galphoterr = obs['maggies_unc'] * maggies2mJy
    else:
        galphot = -2.5 * np.log10(obs['maggies'])
        galphoterr = 2.5 * obs['maggies_unc'] / obs['maggies'] / np.log(10)

    return weff, fwhm, galphot, galphoterr

def _sed(model, theta, obs, sps):
    """Construct the SED for a given set of parameters.  Divide by mextra to account
    for the *current* mass in stars (rather than the integrated stellar mass
    based on the SFH.

    Also convert wavelengths from Angstroms to microns and fluxes from maggies
    to mJy.

    """
    modelwave = sps.wavelengths * (1 + obs['redshift']) # [observed-frame wavelengths]
    modelwave *= ang2micron
    
    modelspec, modelphot, mextra = model.mean_model(theta, obs, sps=sps)
    if False:
        modelspec *= maggies2mJy
        modelphot *= maggies2mJy
    else:
        modelspec = -2.5 * np.log10(modelspec)
        modelphot = -2.5 * np.log10(modelphot)
    #print(modelphot)
    
    return modelwave, modelspec, modelphot

def bestfit_sed(obs, chain=None, lnprobability=None, theta=None, sps=None,
                model=None, seed=None, nrand=100, png=None, priors=None):
    """Plot the (photometric) best-fitting SED.

    Either pass chain and lnprobability (to visualize the emcee fitting results)
    *or* theta (to visualize just a single SED fit).

    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, ScalarFormatter, FuncFormatter, FormatStrFormatter, NullFormatter
    
    import seaborn as sns
    sns.set(style='ticks', font_scale=1.6, palette='Set2')
    
    rand = np.random.RandomState(seed)

    # Get the galaxy photometry and filter info.
    weff, fwhm, galphot, galphoterr = _galaxyphot(obs)

    # Build the maximum likelihood model fit and also grab a random sampling of
    # the chains with weight equal to the posterior probability.    
    if chain is not None:
        if chain.ndim == 3: # emcee
            nwalkers, niter, nparams = chain.shape
            ntot = nwalkers * niter
            flatchain = chain.reshape(ntot, nparams)
            lnp = lnprobability.reshape(ntot)
        else: # dynesty
            ntot, nparams = chain.shape
            flatchain = chain
            lnp = lnprobability
            
        theta = flatchain[lnp.argmax(), :] # maximum likelihood values
        print('Maximum likelihood values: ', theta)

        prob = np.exp(lnp - lnp.max())
        prob /= prob.sum()
        rand_indx = rand.choice(ntot, size=nrand, replace=False, p=prob)
        theta_rand = flatchain[rand_indx, :]

    print('Rendering the maximum-likelihood model...', end='')
    t0 = time.time()
    modelwave, modelspec, modelphot = _sed(model=model, theta=theta, obs=obs, sps=sps)
    print('...took {:.2f} sec'.format(time.time()-t0))
    #print(modelspec.min(), modelspec.max())

    # Establish the wavelength and flux limits.
    minwave, maxwave = 0.3, 2.0 #0.1, 1000.0
    #minwave, maxwave = np.min(weff - 5*fwhm), np.max(weff + fwhm)

    inrange = (modelwave > minwave) * (modelwave < maxwave)
    #maxflux = np.hstack( (galphot + 5*galphoterr, modelspec[inrange]) ).max() * 1.2
    #minflux = -0.05 * maxflux
    minflux, maxflux = (17, 24)#(9, 24)

    fig, ax = plt.subplots(figsize=(8, 6))
    if chain is not None and nrand > 0:
        for ii in range(nrand):
            _, r_modelspec, _ = _sed(model=model, theta=theta_rand[ii, :], obs=obs, sps=sps)
            ax.plot(modelwave, r_modelspec, alpha=0.8, color='gray')
    ax.plot(modelwave, modelspec, alpha=1.0, label='Model spectrum', color='k')
    
    ax.errorbar(weff, modelphot, marker='s', ls='', lw=3, markersize=15, markerfacecolor='none',
                markeredgewidth=3, alpha=0.6, label='Model photometry')
    ax.errorbar(weff, galphot, yerr=galphoterr, marker='o', ls='', lw=2, markersize=10,
                markeredgewidth=2, alpha=0.8, label='Observed photometry',
                elinewidth=2, capsize=5)
                
    ax.set_xlabel(r'Observed-Frame Wavelength (${}$m)'.format('\mu'))
    ax.set_ylabel('Flux (AB mag)')
    #ax.set_ylabel('Flux Density (mJy)')
    ax.set_xlim(minwave, maxwave)
    ax.set_ylim(minflux, maxflux)
    ax.set_xscale('log')
    ax.invert_yaxis()
    #locs, labels = plt.xticks()
    #print(locs, labels)
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #ax.xaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks([0.3, 0.4, 0.6,1.0,2.0])
    ax.xaxis.set_ticklabels(['0.3', '0.4', '0.6', '1.0', '2.0'])
    #ax.xaxis.set_ticklabels([r'$0.3$',r'$0.4$','',r'$0.6$','','','',r'$1.0$', 
    #    r'$2.0$'],minor=True)
    #ax.set_yscale('log')
    #ax.legend(loc='upper right', fontsize=16, frameon=True)
    # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    #ax.get_xaxis().set_major_formatter(FuncFormatter(lambda y, _: '{:.16g}'.format(y)))
    #ax.get_xaxis().set_major_formatter(FuncFormatter(lambda y, _: '{:.1f}'.format(y)))
    #ax.get_xaxis().set_minor_formatter(FuncFormatter(lambda y, _: '{:.1f}'.format(y)))
    #ax.get_xaxis().set_major_formatter(ScalarFormatter())
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)

    loc = png.find('J')
    name = png[loc:loc+5]
    model = png[loc+6:loc+9]
    print('model for plot: ', model)
    if model=='ssp':
        name = name+' SSP'
    if model=='bur':
        name = name+' bursty'
    print(png)
    plt.text(1,16,name)

    # Add an inset with the posterior probability distribution.
    ax1 = fig.add_axes([0.23, 0.68, 0.22, 0.22])
    if priors=='ssp':
        ax1.hist(chain[:, 2], bins=40, histtype='step', linewidth=2, 
             edgecolor='k',fill=True)
    else:
        ax1.hist(chain[:, 3], bins=50, histtype='step', linewidth=2, 
             edgecolor='k',fill=True)
    ax1.set_xlim(9.3, 11.3)
    ax1.set_yticklabels([])
    ax1.set_xlabel(r'$\log_{10}(\mathcal{M}/\mathcal{M}_{\odot})$')
    ax1.set_ylabel(r'$P(\mathcal{M})$')
    ax1.xaxis.set_major_locator(MultipleLocator(0.5))
    #ax1.set_ylim([0,1.0])
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(16)

    if png:
        print('Writing {}'.format(png))
        fig.savefig(png)

def subtriangle(results, showpars=None, truths=None, start=0, thin=2,
                chains=slice(None), logify=None, extents=None, png=None,
                **kwargs):
    """Make a triangle plot of the (thinned, latter) samples of the posterior
    parameter space.  Optionally make the plot only for a supplied subset of
    the parameters.

    :param start:
        The iteration number to start with when drawing samples to plot.

    :param thin:
        The thinning of each chain to perform when drawing samples to plot.

    :param showpars:
        List of string names of parameters to include in the corner plot.

    :param truths:
        List of truth values for the chosen parameters
    
    """
    import corner as triangle

    # Get the ull out the parameter names and flatten the thinned chains.
    parnames = np.array(results['theta_labels'])
    print(parnames)

    # Restrict to a particular set of parameters.
    if showpars:
        ind_show = np.array([parnames.tolist().index(p) for p in showpars])
        parnames = parnames[ind_show]
    else:
        ind_show = slice(None)

    # Get the arrays we need (trace, wghts)
    trace = results['chain'][..., ind_show]
    if trace.ndim == 2:
        trace = trace[None, :]
    trace = trace[chains, start::thin, :]
    wghts = results.get('weights', None)
    if wghts is not None:
        wghts = wghts[start::thin]
    samples = trace.reshape(trace.shape[0] * trace.shape[1], trace.shape[2])

    # logify some parameters
    xx = samples.copy()
    if truths is not None:
        xx_truth = np.array(truths).copy()
    else:
        xx_truth = None
    if logify:
        for p in logify:
            if p in parnames:
                idx = parnames.tolist().index(p)
                xx[:, idx] = np.log10(xx[:,idx])
                parnames[idx] = "log({})".format(parnames[idx])
                if truths is not None:
                    xx_truth[idx] = np.log10(xx_truth[idx])

    # Make nice labels.
    niceparnames = _niceparnames(parnames)
        
    # mess with corner defaults
    corner_kwargs = {"plot_datapoints": False, "plot_density": False,
                     "fill_contours": True, "show_titles": True}
    corner_kwargs.update(kwargs)
    
    fig = triangle.corner(xx, labels=niceparnames, truths=xx_truth,
                          #quantiles=[0.25, 0.5, 0.75], weights=wghts,
                          quantiles=[0.16, 0.5, 0.84], weights=wghts,
                          color='k', **corner_kwargs)

    if png:
        print('Writing {}'.format(png))
        fig.savefig(png)

def logmass2mass(logmass=11.0, **extras):
    return 10**logmass

def load_obs(seed=1, nproc=1, nmin=10, verbose=False, sps=None, galaxy=None):
    """Load the photometry    
    """
    import sedpy
    from prospect.utils.obsutils import fix_obs
    from astropy.io import ascii

    # function to return a flux in maggies given an AB magnitude
    def flux(mag):
        flux = 10. ** (mag / (-2.5))
        return flux

    # function to return an inverse variance in maggies given a magnitude and uncertainty
    def ivar(mag, unc):
        flux = 10. ** (mag / (-2.5))
        func = flux / 1.086 * unc
        ivar = 1 / func ** 2
        return ivar
    
    data = ascii.read('../../autogalfit/bgl_phot.dat')
    dtot = ascii.read('../../autogalfit/bgl_tot.dat')

    match = data.field('Galaxy') == galaxy
    
    nuc_flux_f475w = flux(data['m475'][match][0] - data['ebv'][match][0] * 3.248)
    nuc_flux_f814w = flux(data['m814'][match][0] - data['ebv'][match][0] * 1.536)
    nuc_flux_f160w = flux(data['m160'][match][0] - data['ebv'][match][0] * 0.512)

    nuc_unc_f475w = nuc_flux_f475w / 1.086 * data['u475'][match][0]
    nuc_unc_f814w = nuc_flux_f814w / 1.086 * data['u814'][match][0]
    nuc_unc_f160w = nuc_flux_f160w / 1.086 * data['u160'][match][0]
    
    tot_flux_f475w = flux(dtot['m475'][match][0] - dtot['ebv'][match][0] * 3.248)
    tot_flux_f814w = flux(dtot['m814'][match][0] - dtot['ebv'][match][0] * 1.536)
    tot_flux_f160w = flux(dtot['m160'][match][0] - dtot['ebv'][match][0] * 0.512)

    tot_unc_f475w = tot_flux_f475w / 1.086 * np.sqrt((dtot['u475'][match][0])**2 + 0.03**2)
    tot_unc_f814w = tot_flux_f814w / 1.086 * np.sqrt((dtot['u814'][match][0])**2 + 0.03**2)
    tot_unc_f160w = tot_flux_f160w / 1.086 * np.sqrt((dtot['u160'][match][0])**2 + 0.03**2)
    
    res_mag_f475w = -2.5*np.log10(tot_flux_f475w-nuc_flux_f475w)
    res_mag_f814w = -2.5*np.log10(tot_flux_f814w-nuc_flux_f814w)
    res_mag_f160w = -2.5*np.log10(tot_flux_f160w-nuc_flux_f160w)

    #res_unc_f475w = 1.086 * np.sqrt(nuc_unc_f475w**2+tot_unc_f475w**2) / (tot_flux_f475w-nuc_flux_f475w)    
    #res_unc_f814w = 1.086 * np.sqrt(nuc_unc_f814w**2+tot_unc_f814w**2) / (tot_flux_f814w-nuc_flux_f814w)
    #res_unc_f160w = 1.086 * np.sqrt(nuc_unc_f160w**2+tot_unc_f160w**2) / (tot_flux_f160w-nuc_flux_f160w)
    res_unc_f475w = 1.086 * tot_unc_f475w / (tot_flux_f475w-nuc_flux_f475w)
    res_unc_f814w = 1.086 * tot_unc_f814w / (tot_flux_f814w-nuc_flux_f814w)
    res_unc_f160w = 1.086 * tot_unc_f160w / (tot_flux_f160w-nuc_flux_f160w)    

    
    print('res mags: ', res_mag_f475w, res_unc_f475w, res_mag_f814w, res_unc_f814w, res_mag_f160w, res_unc_f160w)
    print('nuc frac: ', nuc_flux_f475w / tot_flux_f475w, nuc_flux_f814w / tot_flux_f814w, nuc_flux_f160w / tot_flux_f160w)
    
    #unc_mag_f475w = np.sqrt((data['u475'][match][0])**2+(dtot['u475'][match][0])**2)
    #unc_mag_f814w = np.sqrt((data['u814'][match][0])**2+(dtot['u814'][match][0])**2)
    #unc_mag_f160w = np.sqrt((data['u160'][match][0])**2+(dtot['u160'][match][0])**2)

    #phot = dict(
    #    uvis_f475w=(flux(data['m475'][match][0] - data['ebv'][match][0] * 3.248),
    #                ivar(data['m475'][match][0] - data['ebv'][match][0] * 3.248, data['u475'][match][0])), 
    #    uvis_f814w=(flux(data['m814'][match][0] - data['ebv'][match][0] * 1.536),
    #                ivar(data['m814'][match][0] - data['ebv'][match][0] * 1.536, data['u814'][match][0])),
    #    ir_f160w=(flux(data['m160'][match][0] - data['ebv'][match][0] * 0.512),
    #                ivar(data['m160'][match][0] - data['ebv'][match][0] * 0.512, data['u160'][match][0])))
    
    phot = dict(
        uvis_f475w=(flux(res_mag_f475w),
                    ivar(res_mag_f475w, res_unc_f475w)), 
        uvis_f814w=(flux(res_mag_f814w),
                    ivar(res_mag_f160w, res_unc_f814w)),
        ir_f160w=(flux(res_mag_f160w),
                    ivar(res_mag_f160w, res_unc_f160w)))
    
    filternames = (['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])

    obs = {}
    obs['redshift'] = data['z'][match][0] 
    obs["filters"] = sedpy.observate.load_filters(filternames)

    obs["maggies"] = np.array([phot[filt][0] for filt in phot.keys()])
    obs["maggies_unc"] = np.array([1/np.sqrt(phot[filt][1]) for filt in phot.keys()])

    # mask out W4
    #obs["phot_mask"] = np.array(['w4' in f.name for f in obs["filters"]])    
    
    # Create a handy vector of effective wavelengths (optional) 
    obs["phot_wave"] = [f.wave_effective for f in obs["filters"]]
    obs["wavelength"] = None # spectral wavelength
    obs["spectrum"] = None
    obs['unc'] = None  # spectral uncertainties are given here
    #obs['mask'] = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    obs = fix_obs(obs)

    run_params = {}
    run_params['redshift'] = obs['redshift']
    run_params['verbose'] = verbose
    run_params['debug'] = False
    run_params['seed'] = seed
    run_params['nproc'] = nproc
    run_params['param_file'] = '' # no parameter file

    run_params['min_method'] = 'lm'
    run_params['nmin'] = 1
    
    if sps:
        run_params['sps_libraries'] = sps.ssp.libraries

    # dynesty Fitter parameters
    dyn_params = {
        'nested_bound': 'multi',  # bounding method
        'nested_sample': 'unif', # 'unif', 'slice' # sampling method
        'nested_nlive_init': 100,
        'nested_nlive_batch': 100,
        'nested_bootstrap': 0,
        'nested_dlogz_init': 0.05,
        'nested_weight_kwargs': {"pfrac": 1.0},
        #'nested_stop_kwargs': {"post_thresh": 0.05}
        }
    run_params.update(dyn_params)
    
    return obs, run_params

def load_sps(zcontinuous=1, verbose=False):
    """zcontinuous - interpolate between metallicity values.

    """
    from prospect.sources import CSPSpecBasis

    if verbose:
        print('Loading SPS models...', end='')
    t0 = time.time()
    sps = CSPSpecBasis(zcontinuous=zcontinuous)
    if verbose:
        print('...took {:.2f} sec'.format(time.time()-t0))
    return sps

def load_model(obs, template_library='delayed-tau', verbose=False):
    """
    http://dfm.io/python-fsps/current/stellarpop_api/#api-reference
    https://github.com/moustakas/siena-astrophysics/blob/master/research/redmapper/redmapper-stellar-mass.py#L125-L197    
    
    """
    from prospect.models import priors
    from prospect.models.sedmodel import SedModel
    from prospect.models.templates import TemplateLibrary
    from prospect.models.transforms import dustratio_to_dust1

    def base_delayed_tau():
        model_params = TemplateLibrary['parametric_sfh']

        # Initialize with sensible numbers.
        model_params['tau']['init'] = 10.
        model_params['tage']['init'] = 2.0 # --> 0.1
        #model_params['logzsol']['init'] = 0.2

        model_params['logzsol'] = {"N": 1, "isfree": False, "init": 0}

        
        # optimize log-stellar mass, not linear stellar mass
        model_params['logmass'] = {'N': 1, 'isfree': True, 'init': 10.5, # 11 --> 10
                                   'prior': priors.TopHat(mini=8.5, maxi=12.5), # 10 --> 8
                                   'units': '$M_{\odot}$'}

        model_params['mass']['isfree'] = False
        model_params['mass']['init'] = 10**model_params['logmass']['init']
        model_params['mass']['prior'] = None
        model_params['mass']['depends_on'] = logmass2mass
        
        # Adjust the prior ranges.
        model_params['tau']['prior'] = priors.LogUniform(mini=0.1, maxi=30.0) # 0.01 --> 0.001
        model_params['tage']['prior'] = priors.LogUniform(mini=1.0, maxi=7.0) # 0.01 --> 0.001, 10.0 --> 0.1
        #model_params['logzsol']['prior'] = priors.TopHat(mini=-0.5, maxi=0.3)

        add_duste = {"N": 1, "isfree": False, "init": False}
        model_params["add_dust_emission"] = add_duste

        
        #print('HACK!!!!!!!!!!!!!')
        #model_params['tau']['isfree'] = False
        #model_params['tage']['isfree'] = False
        #model_params['logzsol']['isfree'] = False
        #model_params['dust2']['isfree'] = False

        return model_params

    def base_ssp():
        model_params = TemplateLibrary['ssp']
        model_params['tage']['init'] = 2.0
        #model_params['logzsol']['init'] = 0.2

        model_params['logzsol'] = {"N": 1, "isfree": False, "init": 0}
        
        model_params['logmass'] = {'N': 1, 'isfree': True, 'init': 10.5, # 11 --> 10
                                   'prior': priors.TopHat(mini=8.5, maxi=12.5), # 10 --> 8
                                   'units': '$M_{\odot}$'}
        model_params['mass']['isfree'] = False
        model_params['mass']['init'] = 10**model_params['logmass']['init']
        model_params['mass']['prior'] = None
        model_params['mass']['depends_on'] = logmass2mass

        # Adjust prior ranges
        model_params['tage']['prior'] = priors.LogUniform(mini=1.0, maxi=7.0) # 0.01 --> 0.001, 10.0 --> 0.1
        #model_params['logzsol']['prior'] = priors.TopHat(mini=-0.5, maxi=0.3)

        add_duste = {"N": 1, "isfree": False, "init": False}
        model_params["add_dust_emission"] = add_duste
        
        return model_params
        
    if template_library == 'ssp':
        model_params = base_ssp()
        #print(model_params)
    
    if template_library == 'delayed-tau':
        # Underlying delayed tau model.
        model_params = base_delayed_tau()

    if template_library == 'bursty':
        # Underlying delayed tau model.
        model_params = base_delayed_tau()
        
        # Add bursts
        model_params.update(TemplateLibrary['burst_sfh'])
        
        model_params['fburst']['isfree'] = True
        model_params['fburst']['init'] = 0.1
        model_params['fburst']['prior'] = priors.TopHat(mini=0.0, maxi=0.5)

        model_params['fage_burst']['isfree'] = True
        model_params['fage_burst']['init'] = 0.9
        model_params['fage_burst']['prior'] = priors.TopHat(mini=0.5, maxi=1.0)

    # Add dust emission (with fixed dust SED parameters).
    #model_params.update(TemplateLibrary['dust_emission'])

    model_params['dust2']['init'] = 1.0 # diffuse dust
    model_params['dust2']['prior'] = priors.TopHat(mini=0.0, maxi=4.0)

    # Add more dust flexibility.
    model_params['dust_type'] = {'N': 1, 'isfree': False, 'init': 0, 'units': 'dust model'}
    model_params['dust_index'] = {'N': 1, 'isfree': False, 'init': -0.7,
                                  'units': 'power-law index', 'prior': None}

    model_params['dust1'] = {'N': 1, 'isfree': True, 'init': 1.0, 'prior': priors.TopHat(mini=0.0, maxi=4.0),
                             'units': 'optical depth towards young stars'}#,
    #                         'depends_on': dustratio_to_dust1}
    
    #model_params['dust1'] = {'N': 1, 'isfree': False, 'init': 0.0, 'prior': None,
    #                         'units': 'optical depth towards young stars',
    #                         'depends_on': dustratio_to_dust1}
    #model_params['dust_ratio'] = {'N': 1, 'isfree': True, 'init': 1.0,
    #                              'prior': priors.TopHat(mini=1.0, maxi=10.0),
    #                              'units': 'dust1/dust2 ratio (optical depth to young stars vs diffuse)'}

    ## Add nebular emission.
    #model_params.update(TemplateLibrary['nebular'])
    ##model_params['add_neb_continuum']['init'] = False
    #model_params['gas_logu']['init'] = -1.0 # harder radiation field [default is -2.0]

    # Fixed redshift.
    model_params['zred']['init'] = obs['redshift']
    model_params['zred']['isfree'] = False 

    # Change the IMF from Kroupa to Salpeter.
    model_params['imf_type']['init'] = 0
        
    # Now instantiate the model using this new dictionary of parameter specifications
    model = SedModel(model_params)
    if verbose:
        print(model)

    return model

def main():
    """
    Main wrapper script.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--priors', default='delayed-tau', type=str, choices=['delayed-tau', 'bursty', 'ssp'],
                        help='Choose the model priors.')
    parser.add_argument('--galaxy', default='j0826', type=str, help='Output file galaxy.')
    parser.add_argument('--seed', default=1, type=int, help='Seed for random number generation.')
    parser.add_argument('--nproc', default=1, type=int, help='Number of cores to use.')
    parser.add_argument('--sedfit', action='store_true', help='Do the SED fit.')
    parser.add_argument('--qaplots', action='store_true', help='Make pretty plots.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose.')
    args = parser.parse_args()

    datadir = os.getcwd() #os.path.join(os.getenv('HIZEA_PROJECT'), 'j2118-nebula')
    hfile = os.path.join(datadir, 'res-'+'{}-{}.h5'.format(args.galaxy, args.priors))

    if args.sedfit:
        import prospect.io
        import prospect.fitting

        # Initialize the SPS library (takes a bit), the photometry, the "run
        # parameters" dictionary, and the model priors.
        sps = load_sps(verbose=args.verbose)
        obs, rp = load_obs(seed=args.seed, nproc=args.nproc, verbose=args.verbose, sps=sps, galaxy=args.galaxy)
        model = load_model(obs, args.priors, verbose=args.verbose)
        
        #with multiprocessing.Pool(args.nproc) as P:
        output = prospect.fitting.fit_model(obs, model, sps, noise=(None, None),
                                            optimize=False, dynesty=True, emcee=False,
                                            #nested_posterior_thresh=0.05,
                                            pool=None, **rp)
        #output = prospect.fitting.run_emcee_sampler(obs, model, sps, noise=(None, None),
         #                                           optimize = False, dynesty = True, emcee = False,
          #                                          # nested_posterior_thresh=0.05,
           #                                         pool = None, ** rp)

        if os.path.isfile(hfile):
            os.remove(hfile)
        print('Writing {}'.format(hfile))
        prospect.io.write_results.write_hdf5(
            hfile, rp, model, obs, output['sampling'][0],
            output['optimization'][0], tsample=output['sampling'][1],
            toptimize=output['optimization'][1])

    if args.qaplots:
        from prospect.io import read_results as reader
        
        print('Reading {}...'.format(hfile), end='')
        t0 = time.time()
        result, obs, _ = reader.results_from(hfile, dangerous=False)
        print('...took {:.2f} sec'.format(time.time()-t0))

        png = os.path.join(datadir, 'res-'+'{}-{}-corner.png'.format(args.galaxy, args.priors))
        if args.priors=='ssp':
            subtriangle(result, showpars=['logmass', 'tage', 'dust2', 'dust1'], png=png)
                    #logify=['tage'], png=png)

        if args.priors=='delayed-tau':
            subtriangle(result, showpars=['logmass', 'tage', 'tau', 'dust2', 'dust1'],
                    logify=['tau'], png=png)

        if args.priors=='bursty':
            subtriangle(result, showpars=['logmass', 'tage', 'tau', 'dust2', 'dust1', 'fburst', 'fage_burst'], png=png)

        #reader.subcorner(result, start=0, thin=1, fig=plt.subplots(5,5,figsize=(27,27))[0])

        pdb.set_trace()

        # SED
        sps = load_sps(verbose=args.verbose)
        model = load_model(obs, args.priors, verbose=args.verbose)

        png = os.path.join(datadir, 'res-'+'{}-{}-sed.png'.format(args.galaxy, args.priors))
        bestfit_sed(obs, chain=result['chain'], lnprobability=result['lnprobability'], 
                    sps=sps, model=model, seed=1, nrand=100, png=png, priors=args.priors)
    
if __name__ == '__main__':
    main()
