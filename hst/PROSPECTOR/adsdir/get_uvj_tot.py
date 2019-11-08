#!/usr/bin/env python
"""
Fit the SED using Prospector.

run get_uvj_nuc.py --galaxy='J0901' 

"""
import os, time, argparse, pdb
import numpy as np
import multiprocessing
from sedpy.observate import getSED
import sedpy

ang2micron = 1e-4 # Angstrom --> micron
maggies2mJy = 10**(0.4*16.4) # maggies --> mJy

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
                model=None, seed=None, nrand=100, png=None):
    """Plot the (photometric) best-fitting SED.
    Either pass chain and lnprobability (to visualize the emcee fitting results)
    *or* theta (to visualize just a single SED fit).
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, ScalarFormatter, FuncFormatter
    
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

    # new things to calculate UVJ
    tobs = obs
    filternames = (['bessell_U','bessell_V','twomass_J'])
    tobs["filters"] = sedpy.observate.load_filters(filternames)
    
    #filternames = (['bessell_U','bessell_V','twomass_J'])
    filters = tobs["filters"]
    fnames = [f.name for f in filters]
    unique_names, uinds, filter_ind = np.unique(fnames, return_index=True, return_inverse=True)
    unique_filters = np.array(filters)[uinds]

    modelspec_cgs = 10.**(modelspec/(-2.5))*3631 * 1e-23 
    c = 2.998e18 # Angstrom / s
    modelwave_rest_aa = modelwave * 1e4 / (1+obs['redshift'])
    nu = c / (modelwave_rest_aa)
    modelspec_fl = modelspec_cgs * nu**2/c

    mags = getSED(modelwave_rest_aa, modelspec_fl, unique_filters)

    wteff = np.array([f.wave_effective for f in tobs['filters']]) * ang2micron * (1. + obs['redshift'])

    #print(phot)
    print(mags)
    print(wteff)
    print(obs['redshift'])

    pdb.set_trace()

    
    # Establish the wavelength and flux limits.
    minwave, maxwave = 0.1, 40
    #minwave, maxwave = np.min(weff - 5*fwhm), np.max(weff + fwhm)

    inrange = (modelwave > minwave) * (modelwave < maxwave)
    #maxflux = np.hstack( (galphot + 5*galphoterr, modelspec[inrange]) ).max() * 1.2
    #minflux = -0.05 * maxflux
    minflux, maxflux = (11, 24)

    fig, ax = plt.subplots(figsize=(8, 6))

    r_mags = np.zeros([nrand, len(filternames)])
    if chain is not None and nrand > 0:
        for ii in range(nrand):
            _, r_modelspec, _ = _sed(model=model, theta=theta_rand[ii, :], obs=obs, sps=sps)
            ax.plot(modelwave, r_modelspec, alpha=0.8, color='gray')
            r_modelspec_cgs = 10.**(r_modelspec/(-2.5))*3631 * 1e-23 
            r_modelspec_fl = r_modelspec_cgs * nu**2/c
            r_mags[ii] = getSED(modelwave_rest_aa, r_modelspec_fl, unique_filters)
            ax.scatter(wteff, r_mags[ii])
            #print(r_mags[ii])

    umag = np.median(r_mags[:nrand,0])
    ustd = np.std(r_mags[:nrand,0])
    vmag = np.median(r_mags[:nrand,1])
    vstd = np.std(r_mags[:nrand,1])
    jmag = np.median(r_mags[:nrand,2])
    jstd = np.std(r_mags[:nrand,2])

    uvj_phot = np.array([umag, vmag, jmag])
    uvj_err = np.array([ustd, vstd, jstd])

    print(umag, ustd, vmag, vstd, jmag, jstd)
    ax.errorbar(wteff, uvj_phot, yerr=uvj_err, marker='o', ls='', lw=2, markersize=10,
                markeredgewidth=2, alpha=0.8, label='UVJ photometry',
                elinewidth=2, capsize=5)
    
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
    #ax.set_yscale('log')
    #ax.legend(loc='upper right', fontsize=16, frameon=True)
    # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax.get_xaxis().set_major_formatter(FuncFormatter(lambda y, _: '{:.16g}'.format(y)))
    #ax.get_xaxis().set_major_formatter(ScalarFormatter())
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)

    # Add an inset with the posterior probability distribution.
    ax1 = fig.add_axes([0.23, 0.68, 0.22, 0.22])
    ax1.hist(chain[:, 4], bins=50, histtype='step', linewidth=2, 
             edgecolor='k',fill=True)    
    ax1.set_xlim(10.5, 11.5)
    ax1.set_yticklabels([])
    ax1.set_xlabel(r'$\log_{10}(\mathcal{M}/\mathcal{M}_{\odot})$')
    ax1.set_ylabel(r'$P(\mathcal{M})$')
    ax1.xaxis.set_major_locator(MultipleLocator(0.5))
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(16)

    if png:
        print('Writing {}'.format(png))
        fig.savefig(png)

def logmass2mass(logmass=11.0, **extras):
    return 10**logmass


def load_obs(seed=1, nproc=1, nmin=10, verbose=False, sps=None, galaxy=None):
    #Load the photometry

    # import relevant modules
    import sedpy
    from prospect.utils.obsutils import fix_obs
    import numpy as np
    import argparse
    from astropy.io import ascii

    # Here is photometric information for one galaxy.
    # This includes the following:
    # [0] flux in units of maggies
    # [1] inverse variance in units of maggies
    # [2] effective wavelength in units of microns

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

    # read in data from a table
    table = ascii.read('../../umeh_table.dat')
    print("now printing table")
    print(table)
    print("done printing table ")

    # match to the galaxy you want
    match = table.field('Galaxy') == galaxy
    print("now printing match")
    print(match)
    print("done printing match")

    # create a photometry dictionary
    phot = dict(
        FUV=(flux(table.field('fuv_mag')[match][0] - table.field('ebv')[match][0] * 6.783),
             ivar(table.field('fuv_mag')[match][0] - table.field('ebv')[match][0] * 6.783, table.field('fuv_unc')[match][0]),
             0.1528),
        NUV=(flux(table.field('nuv_mag')[match][0] - table.field('ebv')[match][0] * 6.620),
             ivar(table.field('nuv_mag')[match][0] - table.field('ebv')[match][0] * 6.620, table.field('nuv_unc')[match][0]),
             0.2271),
        u=(flux(table.field('u_mag')[match][0] - table.field('ebv')[match][0] * 4.214),
           ivar(table.field('u_mag')[match][0] - table.field('ebv')[match][0] * 4.214, table.field('u_unc')[match][0]),
           0.3543),
        g=(flux(table.field('g_mag')[match][0] - table.field('ebv')[match][0] * 3.269),
           ivar(table.field('g_mag')[match][0] - table.field('ebv')[match][0] * 3.269, table.field('u_unc')[match][0]),
           0.4770),
        r=(flux(table.field('r_mag')[match][0] - table.field('ebv')[match][0] * 2.270),
           ivar(table.field('r_mag')[match][0] - table.field('ebv')[match][0] * 2.270, table.field('u_unc')[match][0]),
           0.6231),
        i=(flux(table.field('i_mag')[match][0] - table.field('ebv')[match][0] * 1.689),
           ivar(table.field('i_mag')[match][0] - table.field('ebv')[match][0] * 1.689, table.field('u_unc')[match][0]),
           0.7625),
        z=(flux(table.field('z_mag')[match][0] - table.field('ebv')[match][0] * 1.261),
           ivar(table.field('z_mag')[match][0] - table.field('ebv')[match][0] * 1.261, table.field('u_unc')[match][0]),
           0.9134),
        w1=(flux(table.field('w1_mag')[match][0] - table.field('ebv')[match][0] * 0.186) * 306.681 / 3631,
            ivar(table.field('w1_mag')[match][0] - table.field('ebv')[match][0] * 0.186, table.field('w1_unc')[match][0]) * (3631 / 306.681) ** 2,
            3.368),
        w2=(flux(table.field('w2_mag')[match][0] - table.field('ebv')[match][0] * 0.123) * 170.663 / 3631,
            ivar(table.field('w2_mag')[match][0] - table.field('ebv')[match][0] * 0.123, table.field('w2_unc')[match][0]) * (3631 / 170.663) ** 2,
            4.618))
        #w3=(flux(table.field('w3_mag')[match][0]) * 29.0448 / 3631,
        #    ivar(table.field('w3_mag')[match][0], table.field('w3_unc')[match][0]) * (3631 / 29.0448) ** 2,
        #    12.082),
        #w4=(flux(table.field('w4_mag')[match][0]) * 8.2839 / 3631,
        #    ivar(table.field('w4_mag')[match][0], table.field('w4_unc')[match][0]) * (3631 / 8.2839) ** 2,
        #    22.194))

    galex = ['galex_FUV', 'galex_NUV']
    sdss = ['sdss_{}0'.format(b) for b in ['u', 'g', 'r', 'i', 'z']]
    #spitzer = ['spitzer_irac_ch{}'.format(n) for n in ['1', '2']]
    #wise = ['wise_w{}'.format(n) for n in ['1', '2', '3', '4']]
    wise = ['wise_w{}'.format(n) for n in ['1', '2']]
    filternames = galex + sdss + wise

    #filternames = (['wfc3_uvis_f475w', 'wfc3_uvis_f814w', 'wfc3_ir_f160w'])

    obs = {}
    obs['redshift'] = table['z'][match][0]
    obs["filters"] = sedpy.observate.load_filters(filternames)

    obs["maggies"] = np.array([phot[filt][0] for filt in phot.keys()])
    obs["maggies_unc"] = np.array([1 / np.sqrt(phot[filt][1]) for filt in phot.keys()])

    # mask out W4
    # obs["phot_mask"] = np.array(['w4' in f.name for f in obs["filters"]])

    # Create a handy vector of effective wavelengths (optional)
    obs["phot_wave"] = [f.wave_effective for f in obs["filters"]]
    obs["wavelength"] = None  # spectral wavelength
    obs["spectrum"] = None
    obs['unc'] = None  # spectral uncertainties are given here
    # obs['mask'] = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    obs = fix_obs(obs)

    run_params = {}
    run_params['redshift'] = obs['redshift']
    run_params['verbose'] = verbose
    run_params['debug'] = False
    run_params['seed'] = seed
    run_params['nproc'] = nproc
    run_params['param_file'] = ''  # no parameter file

    run_params['min_method'] = 'lm'
    run_params['nmin'] = 1

    if sps:
        run_params['sps_libraries'] = sps.ssp.libraries

    # dynesty Fitter parameters
    dyn_params = {
        'nested_bound': 'multi',  # bounding method
        'nested_sample': 'unif',  # 'unif', 'slice' # sampling method
        'nested_nlive_init': 100,
        'nested_nlive_batch': 100,
        'nested_bootstrap': 0,
        'nested_dlogz_init': 0.05,
        'nested_weight_kwargs': {"pfrac": 1.0},
        # 'nested_stop_kwargs': {"post_thresh": 0.05}
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
        model_params['tau']['init'] = 10.0
        model_params['tage']['init'] = 2.0
        model_params['logzsol']['init'] = 0.2

        # optimize log-stellar mass, not linear stellar mass
        model_params['logmass'] = {'N': 1, 'isfree': True, 'init': 11.0,
                                   'prior': priors.TopHat(mini=10.0, maxi=12.0),
                                   'units': '$M_{\odot}$'}

        model_params['mass']['isfree'] = False
        model_params['mass']['init'] = 10**model_params['logmass']['init']
        model_params['mass']['prior'] = None
        model_params['mass']['depends_on'] = logmass2mass
        
        # Adjust the prior ranges.
        model_params['tau']['prior'] = priors.LogUniform(mini=0.1, maxi=30.0)
        model_params['tage']['prior'] = priors.LogUniform(mini=1.0, maxi=7.0)
        model_params['logzsol']['prior'] = priors.TopHat(mini=-0.5, maxi=0.3)

        #print('HACK!!!!!!!!!!!!!')
        #model_params['tau']['isfree'] = False
        #model_params['tage']['isfree'] = False
        #model_params['logzsol']['isfree'] = False
        #model_params['dust2']['isfree'] = False

        return model_params

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
    add_duste = {"N": 1, "isfree": False, "init": False}
    model_params["add_dust_emission"] = add_duste
        
    model_params['dust2']['init'] = 1.0 # diffuse dust
    model_params['dust2']['prior'] = priors.TopHat(mini=0.0, maxi=4.0)

    # Add more dust flexibility.
    model_params['dust_type'] = {'N': 1, 'isfree': False, 'init': 0, 'units': 'dust model'}
    model_params['dust_index'] = {'N': 1, 'isfree': False, 'init': -0.7,
                                  'units': 'power-law index', 'prior': None}
    
    model_params['dust1'] = {'N': 1, 'isfree': False, 'init': 0.0, 'prior': None,
                             'units': 'optical depth towards young stars',
                             'depends_on': dustratio_to_dust1}
    model_params['dust_ratio'] = {'N': 1, 'isfree': True, 'init': 1.0,
                                  'prior': priors.TopHat(mini=1.0, maxi=10.0),
                                  'units': 'dust1/dust2 ratio (optical depth to young stars vs diffuse)'}

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
    parser.add_argument('--priors', default='bursty', type=str, choices=['delayed-tau', 'bursty'],
                        help='Choose the model priors.')
    parser.add_argument('--galaxy', default='J0901', type=str, help='Galaxy name.')
    parser.add_argument('--verbose', action='store_true', help='Be verbose.')
    args = parser.parse_args()

    j2118dir = os.path.join(os.getenv('HIZEA_PROJECT'), 'j2118-nebula')
    hfile = os.path.join(j2118dir, 'tot-{}-{}.h5'.format(args.galaxy, args.priors))

    from prospect.io import read_results as reader
    
    print('Reading {}...'.format(hfile), end='')
    t0 = time.time()
    result, obs, _ = reader.results_from(hfile, dangerous=False)
    print('...took {:.2f} sec'.format(time.time()-t0))

    # SED
    sps = load_sps(verbose=args.verbose)
    model = load_model(obs, args.priors, verbose=args.verbose)

    png = os.path.join(j2118dir, 'uvj-tot-{}-{}-sed.png'.format(args.galaxy, args.priors))
    bestfit_sed(obs, chain=result['chain'], lnprobability=result['lnprobability'], 
                sps=sps, model=model, seed=1, nrand=100, png=png)
    
if __name__ == '__main__':
    main()
