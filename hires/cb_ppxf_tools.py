from time import perf_counter as clock
from os import path
import importlib


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from uncertainties import ufloat

import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib

import importlib
from os import listdir
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

#Vars
cbfits = importlib.import_module('cb_galaxy_fits_sfr_analysis_tools')


#Functions
##############################################################################

def get_z(file):
    hdu = fits.open(file)
    full_gal = cbfits.get_plot_title(file)
    short_n = full_gal[full_gal.find('J'):full_gal.find('J')+5]
    problematic = ['J0944', 'J1341']
    redshifts = {'J0944': 0.51383984, 'J1341': 0.66115844}

    if short_n in problematic:
        z = redshifts.get(short_n)
    elif short_n not in problematic:
        z = float(hdu[2].data['Z']) #replacing 'SPECOBJ' with index 2 since galaxy J1506 doesn't have name SPECOBJ

    return z



def ppxf_example_population_gas_sdss(file,tie_balmer, limit_doublets):
    ppxf_dir = path.dirname(path.realpath(ppxf_package.__file__))

    # Read SDSS DR8 galaxy spectrum taken from here http://www.sdss3.org/dr8/
    # The spectrum is *already* log rebinned by the SDSS DR8
    # pipeline and log_rebin should not be used in this case.
    #
    # file = ppxf_dir + '/spectra/NGC3522_SDSS_DR8.fits'

    hdu = fits.open(file)
    t = hdu[1].data
    # z = float(hdu[1].header["Z"]) # SDSS redshift estimate

    # Used to check if SDSS has correct redshift
    z = get_z(file)
    print('Z ----------- ',z)


    # Only use the wavelength range in common between galaxy and stellar library.
    coadd = hdu['COADD'].data
    fluxdata = coadd['flux']
    wavelength = 10 ** coadd['loglam']
    mask = (wavelength > 3540 * (1 + z)) & (wavelength < 7409 * (1 + z))
    flux = fluxdata[mask]
    galaxy = flux / np.median(flux)  # Normalize spectrum to avoid numerical issues
    wave = wavelength[mask]

    # Get an array with the uncertainty
    uncert_arr = np.sqrt(1/coadd['ivar'])[mask]
    uncert = np.array([ufloat(flux[x], uncert_arr[x]) for x in range(0, len(uncert_arr))])

    # The SDSS wavelengths are in vacuum, while the MILES ones are in air.
    # For a rigorous treatment, the SDSS vacuum wavelengths should be
    # converted into air wavelengths and the spectra should be resampled.
    # To avoid resampling, given that the wavelength dependence of the
    # correction is very weak, I approximate it with a constant factor.
    #
    wave *= np.median(util.vac_to_air(wave) / wave)

    # The noise level is chosen to give Chi^2/DOF=1 without regularization (REGUL=0).
    # A constant noise is not a bad approximation in the fitted wavelength
    # range and reduces the noise in the fit.
    #
    noise = np.full_like(galaxy, 0.01635)  # Assume constant noise per pixel here

    # The velocity step was already chosen by the SDSS pipeline
    # and we convert it below to km/s
    #
    c = 299792.458  # speed of light in km/s
    velscale = c * np.log(wave[1] / wave[0])  # eq.(8) of Cappellari (2017)
    FWHM_gal = 2.76  # SDSS has an approximate instrumental resolution FWHM of 2.76A.

    # ------------------- Setup templates -----------------------

    pathname = ppxf_dir + '/miles_models/Mun1.30*.fits'

    # The templates are normalized to mean=1 within the FWHM of the V-band.
    # In this way the weights and mean values are light-weighted quantities
    miles = lib.miles(pathname, velscale, FWHM_gal)

    # The stellar templates are reshaped below into a 2-dim array with each
    # spectrum as a column, however we save the original array dimensions,
    # which are needed to specify the regularization dimensions
    #
    reg_dim = miles.templates.shape[1:]
    stars_templates = miles.templates.reshape(miles.templates.shape[0], -1)

    # See the pPXF documentation for the keyword REGUL,
    regul_err = 0.013  # Desired regularization error

    # Estimate the wavelength fitted range in the rest frame.
    lam_range_gal = np.array([np.min(wave), np.max(wave)]) / (1 + z)

    # Construct a set of Gaussian emission line templates.
    # The `emission_lines` function defines the most common lines, but additional
    # lines can be included by editing the function in the file ppxf_util.py.
    gas_templates, gas_names, line_wave = util.emission_lines(
        miles.log_lam_temp, lam_range_gal, FWHM_gal,
        tie_balmer=tie_balmer, limit_doublets=limit_doublets)

    # Combines the stellar and gaseous templates into a single array.
    # During the PPXF fit they will be assigned a different kinematic
    # COMPONENT value
    #
    templates = np.column_stack([stars_templates, gas_templates])

    # -----------------------------------------------------------

    # The galaxy and the template spectra do not have the same starting wavelength.
    # For this reason an extra velocity shift DV has to be applied to the template
    # to fit the galaxy spectrum. We remove this artificial shift by using the
    # keyword VSYST in the call to PPXF below, so that all velocities are
    # measured with respect to DV. This assume the redshift is negligible.
    # In the case of a high-redshift galaxy one should de-redshift its
    # wavelength to the rest frame before using the line below as described
    # in PPXF_EXAMPLE_KINEMATICS_SAURON and Sec.2.4 of Cappellari (2017)
    #
    c = 299792.458
    dv = c * (miles.log_lam_temp[0] - np.log(wave[0]))  # eq.(8) of Cappellari (2017)
    vel = c * np.log(1 + z)  # eq.(8) of Cappellari (2017)
    start = [vel, 180.]  # (km/s), starting guess for [V, sigma]

    n_temps = stars_templates.shape[1]
    n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
    n_balmer = len(gas_names) - n_forbidden

    # Assign component=0 to the stellar templates, component=1 to the Balmer
    # gas emission lines templates and component=2 to the forbidden lines.
    component = [0] * n_temps + [1] * n_balmer + [2] * n_forbidden
    gas_component = np.array(component) > 0  # gas_component=True for gas templates

    # Fit (V, sig, h3, h4) moments=4 for the stars
    # and (V, sig) moments=2 for the two gas kinematic components
    moments = [4, 2, 2]

    # Adopt the same starting value for the stars and the two gas components
    start = [start, start, start]

    # If the Balmer lines are tied one should allow for gas reddeining.
    # The gas_reddening can be different from the stellar one, if both are fitted.
    gas_reddening = 0 if tie_balmer else None

    # Here the actual fit starts.
    #
    # IMPORTANT: Ideally one would like not to use any polynomial in the fit
    # as the continuum shape contains important information on the population.
    # Unfortunately this is often not feasible, due to small calibration
    # uncertainties in the spectral shape. To avoid affecting the line strength of
    # the spectral features, we exclude additive polynomials (DEGREE=-1) and only use
    # multiplicative ones (MDEGREE=10). This is only recommended for population, not
    # for kinematic extraction, where additive polynomials are always recommended.
    #
    t = clock()
    pp = ppxf(templates, galaxy, noise, velscale, start,
              plot=False, moments=moments, degree=-1, mdegree=10, vsyst=dv,
              lam=wave, clean=False, regul=1. / regul_err, reg_dim=reg_dim,
              component=component, gas_component=gas_component,
              gas_names=gas_names, gas_reddening=gas_reddening)

    # When the two Delta Chi^2 below are the same, the solution
    # is the smoothest consistent with the observed spectrum.
    #
    print('Desired Delta Chi^2: %.4g' % np.sqrt(2 * galaxy.size))
    print('Current Delta Chi^2: %.4g' % ((pp.chi2 - 1) * galaxy.size))
    print('Elapsed time in PPXF: %.2f s' % (clock() - t))

    weights = pp.weights[~gas_component]  # Exclude weights of the gas templates
    weights = weights.reshape(reg_dim) / weights.sum()  # Normalized

    miles.mean_age_metal(weights)
    miles.mass_to_light(weights, band="r")


    def plot_pp():
        # Plot fit results for stars and gas.
        plt.clf()
        plt.subplot(211)
        pp.plot()

        # Plot stellar population mass fraction distribution
        plt.subplot(212)
        miles.plot(weights)
        plt.tight_layout()
        plt.pause(1)
    plotpp = False
    if plotpp:
        plot_pp()

    return [pp,wave,galaxy,flux,uncert]


#gets flux value (does the sums and stuff). Parameters: file of galaxy, wavelength of emission line from which flux
#is desired, spread in same units as wavelength of the emission line, data generated by fits
def pp_get_flux(wavel,wspread,fits_data):
    print('----WSPREAD ',wspread)
    cbfits = importlib.import_module('cb_galaxy_fits_sfr_analysis_tools')

    pp,wave = fits_data[0],fits_data[1]
    galaxy = fits_data[2]
    fluxdata = fits_data[3]
    #make galaxy not normalized
    # galaxy = fits_data[2]*np.median(fits_data[2])
    gas_fit = pp.gas_bestfit
    gal_fit = pp.bestfit
    continuum = (gal_fit - gas_fit)*np.median(fluxdata)
    contin_u = np.array([ufloat(x, 0.1*x) for x in continuum])
    flux = galaxy
    lam=wave
    print('--------------------WAVEL', wavel)
    print('--------------------LAM', lam, 'TYPE ', type(lam))
    #next we get indeces to mark regions over which we can find the peak of the fit
    em_center_ex = cbfits.cb_match_lambda(lam,wavel)
    em_lowlim_ex = cbfits.cb_match_lambda(lam,wavel - wspread)
    em_maxlim_ex = cbfits.cb_match_lambda(lam,wavel + wspread)
    cont_high_ex = cbfits.cb_match_lambda(lam,wavel + 4*wspread)
    cont_low_ex = cbfits.cb_match_lambda(lam,wavel - 4*wspread)

    #next we find the peak of the fit
    print('GASFIT ',gas_fit,' emlowlim ', em_lowlim_ex, ' emmaxlim ', em_maxlim_ex)
    linefit_peak = np.amax(gas_fit[em_lowlim_ex:em_maxlim_ex])
    wpeak = lam[np.where(gas_fit == linefit_peak)]

    #next we get indeces to mark regions over which we can sum
    em_center = cbfits.cb_match_lambda(lam,wpeak)
    em_lowlim = cbfits.cb_match_lambda(lam,wpeak - wspread)
    em_maxlim = cbfits.cb_match_lambda(lam,wpeak + wspread)
    cont_high = cbfits.cb_match_lambda(lam,wpeak + 4*wspread)
    cont_low = cbfits.cb_match_lambda(lam,wpeak - 4*wspread)

    #continuum = np.mean([np.mean(flux[cont_low:em_lowlim]),np.mean(flux[em_maxlim:cont_high])])
    #continuum_const = np.median(np.concatenate((flux[cont_low:em_lowlim],flux[em_maxlim:cont_high])))
    dlambda = np.array([lam[x+1]-lam[x] for x in range(em_lowlim,em_maxlim+1)])
    #continuum_area = np.sum(float(continuum_const)*dlambda)
    #prelim_flux_integral = np.sum(dlambda*flux[em_lowlim:em_maxlim+1])
    flux_area = np.sum(dlambda*fluxdata[em_lowlim:em_maxlim+1])
    cont_area = np.sum(dlambda * continuum[em_lowlim:em_maxlim+1])
    flux_val = flux_area - cont_area
    #flux_val = prelim_flux_integral-continuum_area

    ##### Uncertainty Calculations #### (separated as they came after the completion of the first stage of this code)
    flux_uncert = fits_data[4]

    # It would be good if we had uncertainty for the ppxf fits to take into account
    # UPDATE take 10% of continuum as the uncertainty
    flux_area_unc = np.sum(dlambda * flux_uncert[em_lowlim:em_maxlim + 1])
    cont_area_unc = np.sum(dlambda * contin_u[em_lowlim:em_maxlim+1])

    flux_val_unc = flux_area_unc - cont_area_unc

    #### /Uncertainty Calculations ####

    ##print tests
    #print(lam[cont_low:cont_high])
    #since we don't deal with continuum constant, we change it to be useful data we can plot
    continuum_const = [pp.bestfit,pp.gas_bestfit,fluxdata]

    plot_data = [em_lowlim,em_maxlim,cont_low,cont_high,flux,lam,continuum_const]
    #flux_vals = [flux_val, flux_val_unc]


    return flux_val_unc, plot_data



##############################################################################
#
# #Variables
# filepath = 'galaxies_fits/spec-0761-54524-0409.fits'
#
# #Execution sequence
#
# cbfits = importlib.import_module('cb_galaxy_fits_sfr_analysis_tools')
#
# data = cbfits.get_quantities(filepath,'H_beta')
# hb_spread, hb_wav = cbfits.wspread(filepath, 'H_beta', 400)
# fits_data = ppxf_example_population_gas_sdss(filepath,tie_balmer=True, limit_doublets=True)
# hb_flux, plot_data = pp_get_flux(hb_wav, hb_spread,fits_data)
# hb_lum = cbfits.get_lum(filepath, hb_flux)
# sfr = cbfits.get_sfr_hbeta(hb_lum)
# with PdfPages('Emissionlines_ppxf.pdf') as pdf:
#     cbfits.plot_area_of_interest(plot_data, [hb_flux, hb_lum, sfr], filepath, ppxf_data=True)
#     pdf.savefig(bbox_inches="tight")
#     plt.close('all')



