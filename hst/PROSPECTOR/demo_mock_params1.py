from copy import deepcopy
import numpy as np
from prospect.models import priors, sedmodel
from prospect.sources import CSPSpecBasis

from sedpy.observate import load_filters

tophat = None

# --------------
# RUN_PARAMS
# --------------

run_params = {'verbose':True,
              'debug':False,
              'outfile':'output/demo_mock',
              # Optimization parameters
              'do_powell': False,
              'ftol':0.5e-5, 'maxfev':5000,
              'initial_disp':0.1,
              'do_levenburg': True,
              'nmin': 10,
              # emcee Fitter parameters
              'nwalkers':128,
              'nburn':[32, 32, 64], 'niter':512,
              # nestle Fitter parameters
              'nestle_method': 'single',
              'nestle_npoints': 200,
              'nestle_maxcall': int(1e6),
              # Mock data parameters
              'snr': 20.0,
              'add_noise': False,
              # Input mock model parameters
              'mass': 5e10,
              'logzsol': 0.19,
              'tage': 7.9,
              'tau': 9.6,
              'dust2': 0.18,
              'zred': 0.6,
              # Data manipulation parameters
              'logify_spectrum':False,
              'normalize_spectrum':False,
              # SPS parameters
              'zcontinuous': 1,
              }

# --------------
# OBS
# --------------

# Here we are going to put together some filter names
# All these filters are available in sedpy.  If you want to use other filters,
# add their transmission profiles to sedpy/sedpy/data/filters/ with appropriate
# names (and format).  See sedpy documentation
galex = ['galex_FUV', 'galex_NUV']
sdss = ['sdss_{0}0'.format(b) for b in ['u','g','r','i','z']]
twomass = ['twomass_{}'.format(b) for b in ['J', 'H', 'Ks']]
spitzer = ['spitzer_irac_ch'+n for n in ['1','2','3','4']]

filterset = ['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w']

def load_obs(snr=10.0, add_noise=True, **kwargs):
    """Make a mock dataset.  Feel free to add more complicated kwargs, and put
    other things in the run_params dictionary to control how the mock is
    generated.

    :param snr:
        The S/N of the phock photometry.  This can also be a vector of same
        lngth as the number of filters.

    :param add_noise: (optional, boolean, default: True)
        If True, add a realization of the noise to the mock spectrum
    """
    # We'll put the mock data in this dictionary, just as we would for real
    # data.  But we need to know which bands (and wavelengths if doing
    # spectroscopy) in which to generate mock data.
    mock = {}
    mock['wavelength'] = None # No spectrum
    #filterset = galex + sdss + twomass + spitzer[:2] # only warm spitzer
    #filterset = ['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w']
    mock['filters'] = load_filters(filterset)

    # We need the models to make a mock
    sps = load_sps(**kwargs)
    mod = load_model(**kwargs)

    # Now we get the mock params from the kwargs dict
    params = {}
    for p in mod.params.keys():
        if p in kwargs:
            params[p] = np.atleast_1d(kwargs[p])

    # And build the mock
    mod.params.update(params)
    spec, phot, _ = mod.mean_model(mod.theta, mock, sps=sps)
    # Now store some output
    mock['true_spectrum'] = spec.copy()
    mock['true_maggies'] = phot.copy()
    mock['mock_params'] = deepcopy(mod.params)
    # And add noise
    pnoise_sigma = phot / snr
    if add_noise:
        pnoise = np.random.normal(0, 1, len(phot)) * pnoise_sigma
        mock['maggies'] = phot + pnoise
    else:
        mock['maggies'] = phot.copy()
    mock['maggies_unc'] = pnoise_sigma
    mock['mock_snr'] = snr
    mock['phot_mask'] = np.ones(len(phot), dtype=bool)

    return mock

# --------------
# SPS Object
# --------------

def load_sps(zcontinuous=1, compute_vega_mags=False, **extras):
    sps = CSPSpecBasis(zcontinuous=zcontinuous,
                       compute_vega_mags=compute_vega_mags)
    return sps

# -----------------
# Noise Model
# ------------------

def load_gp(**extras):
    return None, None

# --------------
# MODEL_PARAMS
# --------------

# You'll note below that we have 5 free parameters:
# mass, logzsol, tage, tau, dust2
# They are all scalars.
# mass and tau have logUniform priors (i.e. TopHat priors in log(mass) and
# log(tau)), the rest have TopHat priors.
# You should adjust the prior ranges (particularly in mass) to suit your objects.
#
# The other parameters are all fixed, but we may want to explicitly set their
# values, which can be done here, to override any defaults in python-FSPS


model_params = []

# --- Distance ---
model_params.append({'name': 'zred', 'N': 1,
                        'isfree': False,
                        'init': 0.1,
                        'units': '',
                        'prior': priors.TopHat(mini=0.0, maxi=4.0)})

# --- SFH --------
# FSPS parameter
model_params.append({'name': 'sfh', 'N': 1,
                        'isfree': False,
                        'init': 4,  # This is delay-tau
                        'units': 'type',
                        'prior': None})

model_params.append({'name': 'mass', 'N': 1,
                        'isfree': True,
                        'init': 1e10,
                        'init_disp': 1e9,
                        'units': r'M_\odot',
                        'prior': priors.LogUniform(mini=1e8, maxi=1e12)})

model_params.append({'name': 'logzsol', 'N': 1,
                        'isfree': True,
                        'init': -0.3,
                        'init_disp': 0.3,
                        'units': r'$\log (Z/Z_\odot)$',
                        'prior': priors.TopHat(mini=-1.5, maxi=0.19)})

# If zcontinuous > 1, use 3-pt smoothing
model_params.append({'name': 'pmetals', 'N': 1,
                        'isfree': False,
                        'init': -99,
                        'prior': None})
                        
# FSPS parameter
model_params.append({'name': 'tau', 'N': 1,
                        'isfree': True,
                        'init': 1.0,
                        'init_disp': 0.5,
                        'units': 'Gyr',
                        'prior':priors.LogUniform(mini=0.101, maxi=100)})

# FSPS parameter
model_params.append({'name': 'tage', 'N': 1,
                        'isfree': True,
                        'init': 5.0,
                        'init_disp': 3.0,
                        'units': 'Gyr',
                        'prior': priors.TopHat(mini=0.101, maxi=14.0)})

model_params.append({'name': 'fage_burst', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': 'time at wich burst happens, as a fraction of `tage`',
                        'prior': priors.TopHat(mini=0.9, maxi=1.0)})

# This function transfroms from a fractional age of a burst to an absolute age.
# With this transformation one can sample in ``fage_burst`` without worry about
# the case tburst > tage.
def tburst_fage(tage=0.0, fage_burst=0.0, **extras):
    return tage * fage_burst

# FSPS parameter
model_params.append({'name': 'tburst', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': '',
                        'prior': None,})
                        #'depends_on': tburst_fage})  # uncomment if using bursts.

# FSPS parameter
model_params.append({'name': 'fburst', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': '',
                        'prior': priors.TopHat(mini=0.0, maxi=0.5)})

# --- Dust ---------
# FSPS parameter
model_params.append({'name': 'dust2', 'N': 1,
                        'isfree': True,
                        'init': 0.35,
                        'reinit': True,
                        'init_disp': 0.3,
                        'units': 'Diffuse dust optical depth towards all stars at 5500AA',
                        'prior': priors.TopHat(mini=0.0, maxi=2.0)})

# FSPS parameter
model_params.append({'name': 'dust1', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': 'Extra optical depth towards young stars at 5500AA',
                        'prior': None,})

# FSPS parameter
model_params.append({'name': 'dust_tesc', 'N': 1,
                        'isfree': False,
                        'init': 7.0,
                        'units': 'definition of young stars for the purposes of the CF00 dust model, log(Gyr)',
                        'prior': None,})

# FSPS parameter
model_params.append({'name': 'dust_index', 'N': 1,
                        'isfree': False,
                        'init': -0.7,
                        'units': 'power law slope of the attenuation curve for diffuse dust',
                        'prior': None,})

# FSPS parameter
model_params.append({'name': 'dust1_index', 'N': 1,
                        'isfree': False,
                        'init': -1.0,
                        'units': 'power law slope of the attenuation curve for young-star dust',
                        'prior': None,})

# FSPS parameter
model_params.append({'name': 'dust_type', 'N': 1,
                        'isfree': False,
                        'init': 0,  # power-laws
                        'units': 'index',
                        'prior': None})

# FSPS parameter
model_params.append({'name': 'add_dust_emission', 'N': 1,
                        'isfree': False,
                        'init': True,
                        'units': 'index',
                        'prior': None})

# An example of the parameters controlling the dust emission SED.  There are others!
model_params.append({'name': 'duste_umin', 'N': 1,
                        'isfree': False,
                        'init': 1.0,
                        'units': 'MMP83 local MW intensity',
                        'prior': None})

# --- Stellar Pops ------------
# One could imagine changing these, though doing so *during* the fitting will
# be dramatically slower.
# FSPS parameter
model_params.append({'name': 'tpagb_norm_type', 'N': 1,
                        'isfree': False,
                        'init': 2,
                        'units': 'index',
                        'prior': None})

# FSPS parameter
model_params.append({'name': 'add_agb_dust_model', 'N': 1,
                        'isfree': False,
                        'init': True,
                        'units': 'index',
                        'prior': None})

# FSPS parameter
model_params.append({'name': 'agb_dust', 'N': 1,
                        'isfree': False,
                        'init': 1,
                        'units': 'index',
                        'prior': None})

# --- Nebular Emission ------

# For speed we turn off nebular emission in the demo
model_params.append({'name': 'add_neb_emission', 'N': 1,
                        'isfree': False,
                        'init': False,
                        'prior': None})

# Here is a really simple function that takes a **dict argument, picks out the
# `logzsol` key, and returns the value.  This way, we can have gas_logz find
# the value of logzsol and use it, if we uncomment the 'depends_on' line in the
# `gas_logz` parameter definition.
#
# One can use this kind of thing to transform parameters as well (like making
# them linear instead of log, or divide everything by 10, or whatever.) You can
# have one parameter depend on several others (or vice versa).  Just remember
# that a parameter with `depends_on` must always be fixed.

def stellar_logzsol(logzsol=0.0, **extras):
    return logzsol

# FSPS parameter
model_params.append({'name': 'gas_logz', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': r'log Z/Z_\odot',
#                        'depends_on': stellar_logzsol,
                        'prior': priors.TopHat(mini=-2.0, maxi=0.5)})

# FSPS parameter
model_params.append({'name': 'gas_logu', 'N': 1,
                        'isfree': False,
                        'init': -2.0,
                        'units': '',
                        'prior': priors.TopHat(mini=-4, maxi=-1)})

# --- Calibration ---------
# Only important if using a NoiseModel
model_params.append({'name': 'phot_jitter', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': 'mags',
                        'prior': priors.TopHat(mini=0.0, maxi=0.2)})

def load_model(zred=0.0, **extras):
    # In principle (and we've done it) you could have the model depend on
    # command line arguments (or anything in run_params) by making changes to
    # `model_params` here before instantiation the SedModel object.  Up to you.

    # Here we are going to set the intial value (and the only value, since it
    # is not a free parameter) of the redshift parameter to whatever was used
    # to generate the mock, listed in the run_params dictionary.
    pn = [p['name'] for p in model_params]
    zind = pn.index('zred')
    model_params[zind]['init'] = zred
    
    return sedmodel.SedModel(model_params)

