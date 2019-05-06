import numpy as np
from prospect.models import priors, sedmodel
from prospect.sources import CSPSpecBasis
from sedpy.observate import load_filters
from astropy.io import ascii

# --------------
# RUN_PARAMS
# --------------

run_params = {'verbose':True,
              'debug':False,
              'outfile':'J2140_ssp',
              # Fitter parameters
              'nwalkers':128,
              'nburn':[10, 10, 10], 'niter':512,
              'do_powell': False,
              'do_levenburg': False,
              'ftol':0.5e-5, 'maxfev':5000,
              'initial_disp':0.1,
              # Obs data parameters
              #'objid':1,
              'phottable': '../hst/autogalfit/flux.dat',
              'logify_spectrum':False,
              'normalize_spectrum':False,
              'wlo':3750., 'whi':7200.,
              # SPS parameters
              'zcontinuous': 1,
              'mass': 1e10,
              'logzsol': 0.19,
              'tage': 0.01,
              'tau': 0.01,
              'dust2': 0.18,
              'zred': 0.6,
              }

# --------------
# OBS
# --------------

# Here we are going to put together some filter names
galex = ['galex_FUV', 'galex_NUV']
spitzer = ['spitzer_irac_ch'+n for n in ['1','2','3','4']]
bessell = ['bessell_'+n for n in ['U', 'B', 'V', 'R', 'I']]
sdss = ['sdss_{0}0'.format(b) for b in ['u','g','r','i','z']]

# The first filter set is Johnson/Cousins, the second is SDSS. We will use a
# flag in the photometry table to tell us which set to use for each object
# (some were not in the SDSS footprint, and therefore have Johnson/Cousins
# photometry)
#
# All these filters are available in sedpy.  If you want to use other filters,
# add their transmission profiles to sedpy/sedpy/data/filters/ with appropriate
# names (and format)
#filtersets = (galex + bessell + spitzer,
#              galex + sdss + spitzer)
filtersets = (['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])
galaxy = 11

def load_obs(objid=galaxy, phottable='../hst/autogalfit/flux.dat', **kwargs):
    """Load photometry from an ascii file.  Assumes the following columns:
    `objid`, `filterset`, [`mag0`,....,`magN`] where N >= 11.  The User should
    modify this function (including adding keyword arguments) to read in their
    particular data format and put it in the required dictionary.

    :param objid:
        The object id for the row of the photomotery file to use.  Integer.
        Requires that there be an `objid` column in the ascii file.

    :param phottable:
        Name (and path) of the ascii file containing the photometry.

    :returns obs:
        Dictionary of observational data.
    """
    # Writes your code here to read data.  Can use FITS, h5py, astropy.table,
    # sqlite, whatever.
    # e.g.:
    # import astropy.io.fits as pyfits
    # catalog = pyfits.getdata(phottable)

    catalog = ascii.read(phottable)
    print(catalog)
    wavelengths = np.array([475,814,1600])
    #mags = np.array(catalog[0]['f_475'], catalog[0]['f_814'], catalog[0]['f_1600'])
    mags = np.array([catalog[objid]['Nanomaggies_{}'.format(i)] for i in wavelengths])
    print(mags)
    #ivars = np.array([catalog['ivar_{}'.format(i)] for i in wavelengths])
    ids = catalog['Galaxy']
    zs = catalog['z']

    # Build output dictionary.
    obs = {}
    # This is a list of sedpy filter objects.    See the
    # sedpy.observate.load_filters command for more details on its syntax.
    obs['filters'] = load_filters(filtersets)
    # This is a list of maggies, converted from mags.  It should have the same
    # order as `filters` above.
    obs['maggies'] = np.squeeze(mags*1e-9)
    # HACK.  You should use real flux uncertainties
    obs['maggies_unc'] = obs['maggies'] * 0.07
    # Here we mask out any NaNs or infs
    obs['phot_mask'] = np.isfinite(np.squeeze(mags))
    # We have no spectrum.
    obs['wavelength'] = None

    # Add unessential bonus info.  This will be stored in output
    #obs['dmod'] = catalog[ind]['dmod']
    obs['objid'] = ids[objid]
    obs['z'] = zs[objid]
    obs['spectrum'] = None
    obs['logify_spectrum'] = False

    return obs

# --------------
# SPS Object
# --------------

def load_sps(zcontinuous=1, compute_vega_mags=False, **extras):
    sps = CSPSpecBasis(zcontinuous=zcontinuous,
                       compute_vega_mags=compute_vega_mags)
    return sps

# -----------------
# Gaussian Process
# ------------------

def load_gp(**extras):
    return None, None

# --------------
# MODEL_PARAMS
# --------------

# You'll note below that we have 5 free parameters:
# mass, logzsol, tage, tau, dust2
# They are all scalars.
#
# The other parameters are all fixed, but we want to explicitly set their
# values, possibly from something differnt than the FSPS defaults

model_params = []

# --- Distance ---
# This is the redshift.  Because we are not separately supplying a ``lumdist``
# parameter, the distance will be determined from the redshift using a WMAP9
# cosmology, unless the redshift is 0, in which case the distance is assumed to
# be 10pc (i.e. for absolute magnitudes)
model_params.append({'name': 'zred', 'N': 1,
                        'isfree': False,
                        'init': 0.6,
                        'units': '',
                        'prior':priors.TopHat(mini=0.0, maxi=4.0)})

# --- SFH --------
# FSPS parameter.  sfh=4 is a delayed-tau SFH
model_params.append({'name': 'sfh', 'N': 1,
                        'isfree': False,
                        'init': 0,
                        'units': 'type'
                    })

# Normalization of the SFH.  If the ``mass_units`` parameter is not supplied,
# this will be in surviving stellar mass.  Otherwise it is in the total stellar
# mass formed.
model_params.append({'name': 'mass', 'N': 1,
                        'isfree': True,
                        'init': 5e9,
                        'init_disp': 3e9,
                        'units': r'M_\odot',
                        'prior':priors.TopHat(mini=1e7, maxi=1e12)})

# Since we have zcontinuous=1 above, the metallicity is controlled by the
# ``logzsol`` parameter.
model_params.append({'name': 'logzsol', 'N': 1,
                        'isfree': False,
                        'init': 0,
                        'init_disp': 0.1,
                        'units': r'$\log (Z/Z_\odot)$',
                        'prior': priors.TopHat(mini=-1, maxi=0.19)})

# FSPS parameter
#model_params.append({'name': 'tau', 'N': 1,
#                        'isfree': True,
#                        'init': 1.0,
#                        'units': 'Gyr',
#                        'prior':priors.LogUniform(mini=0.001, maxi=10)})

# FSPS parameter
model_params.append({'name': 'tage', 'N': 1,
                        'isfree': True,
                        'init': 0.01,
                        'init_disp': 0.01,
                        'units': 'Gyr',
                        #'prior':priors.TopHat(mini=0.01, maxi=10.0)})
                        'prior':priors.LogUniform(mini=0.003, maxi=0.1)})

# FSPS parameter
#model_params.append({'name': 'tburst', 'N': 1,
#                        'isfree': False,
#                        'init': 0.0,
#                        'units': '',
#                        'prior':priors.TopHat(mini=0.0, maxi=1.3)})

# FSPS parameter
#model_params.append({'name': 'fburst', 'N': 1,
#                        'isfree': False,
#                        'init': 0.0,
#                        'units': '',
#                        'prior':priors.TopHat(mini=0.0, maxi=0.5)})

# --- Dust ---------
# FSPS parameter
model_params.append({'name': 'dust1', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': '',
                        'prior':priors.TopHat(mini=0.1, maxi=2.0)})

# FSPS parameter
model_params.append({'name': 'dust2', 'N': 1,
                        'isfree': True,
                        'init': 0.35,
                        'reinit': True,
                        'init_disp': 0.3,
                        'units': '',
                        'prior':priors.TopHat(mini=0.0, maxi=2.0)})

# FSPS parameter
model_params.append({'name': 'dust_index', 'N': 1,
                        'isfree': False,
                        'init': -0.7,
                        'units': '',
                        'prior':priors.TopHat(mini=-1.5, maxi=-0.5)})

# FSPS parameter
model_params.append({'name': 'dust1_index', 'N': 1,
                        'isfree': False,
                        'init': -1.0,
                        'units': '',
                        'prior':priors.TopHat(mini=-1.5, maxi=-0.5)})

# FSPS parameter
model_params.append({'name': 'dust_tesc', 'N': 1,
                        'isfree': False,
                        'init': 7.0,
                        'units': 'log(Gyr)',
                        'prior': None})

# FSPS parameter
model_params.append({'name': 'dust_type', 'N': 1,
                        'isfree': False,
                        'init': 0,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'add_dust_emission', 'N': 1,
                        'isfree': False,
                        'init': True,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'duste_umin', 'N': 1,
                        'isfree': False,
                        'init': 1.0,
                        'units': 'MMP83 local MW intensity'})

# --- Stellar Pops ------------
# FSPS parameter
model_params.append({'name': 'tpagb_norm_type', 'N': 1,
                        'isfree': False,
                        'init': 2,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'add_agb_dust_model', 'N': 1,
                        'isfree': False,
                        'init': True,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'agb_dust', 'N': 1,
                        'isfree': False,
                        'init': 1,
                        'units': 'index'})

# --- Nebular Emission ------

# FSPS parameter
model_params.append({'name': 'add_neb_emission', 'N': 1,
                     'isfree': False,
                     'init': False})

# Here is a really simple function that takes a **dict argument, picks out the
# `logzsol` key, and returns the value.  This way, we can have gas_logz find
# the value of logzsol and use it, if we uncomment the 'depends_on' line in the
# `gas_logz` parameter definition.
#
# One can use this kind of thing to transform parameters as well (like making
# them linear instead of log, or divide everything by 10, or whatever.) You can
# have one parameter depend on several others (or vice versa).  Just remember
# that a parameter with `depends_on` must always be fixed.  It's also not a
# good idea to have one parameter depend on another parameter that *also*
# depends on something, since dependency resolution order is arbitrary.

def stellar_logzsol(logzsol=0.0, **extras):
    return logzsol


# FSPS parameter
model_params.append({'name': 'gas_logz', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': r'log Z/Z_\odot',
#                        'depends_on': stellar_logzsol,
                        'prior':priors.TopHat(mini=-2.0, maxi=0.5)})

# FSPS parameter
model_params.append({'name': 'gas_logu', 'N': 1,
                        'isfree': False,
                        'init': -2.0,
                        'units': '',
                        'prior':priors.TopHat(mini=-4, maxi=-1)})

# --- Calibration ---------
model_params.append({'name': 'phot_jitter', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': 'mags',
                        'prior':priors.TopHat(mini=0.0, maxi=0.2)})

def load_model(**extras):
    # In principle (and we've done it) you could have the model depend on
    # command line arguments (or anything in run_params) by making changes to
    # `model_params` here before instantiation the SedModel object.  Up to you.
    return sedmodel.SedModel(model_params)

