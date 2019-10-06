import importlib
from os import listdir
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt


#Variables
filepath = 'galaxies_fits/spec-0761-54524-0409.fits'
galfolpath = 'C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hires/galaxies_fits/'
fitsfiles = [f for f in listdir(galfolpath)]
cbfits = importlib.import_module('cb_galaxy_fits_sfr_analysis_tools')
ppcbfits = importlib.import_module('cb_ppxf_tools')
velocity_ranges = [60, 820, 20]
vel_space = np.array(np.linspace(60, 800, 38))

#Execution sequence


def get_sfr(filepath, fits_data, vel_range):
    min, max, step = vel_range[0], vel_range[1], vel_range[2]
    sfr_array = np.array([])
    for vel in range(min, max, step):
        # data = cbfits.get_quantities(filepath,'H_beta')
        hb_spread, hb_wav = cbfits.wspread(filepath, 'H_beta', vel)
        hb_flux, plot_data = ppcbfits.pp_get_flux(hb_wav, hb_spread,fits_data)
        hb_lum = cbfits.get_lum(filepath, hb_flux)
        sfr = cbfits.get_sfr_hbeta(hb_lum)
        sfr_array = np.append(sfr_array, sfr)

    return sfr_array

def plot_sfr_vel(sfr_arr, filepath):
    # Vars
    plot_title = cbfits.get_plot_title(filepath)

    # Init
    fig = plt.figure()
    fig.suptitle(plot_title)
    plt.ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")
    plt.xlabel("Velocity " + r"$\frac{km}{s}$")


    # Plotting
    sfrs = np.array([x.n for x in sfr_arr])
    uncert = np.array([x.s for x in sfr_arr])

    plt.scatter(vel_space, sfrs)
    plt.errorbar(vel_space, sfrs, yerr=uncert)




with PdfPages('Vel_vs_SFR_U.pdf') as pdf:
    for file in fitsfiles:
        print('\nGetting SFR Values For: ', file, '..........')
        filepath = galfolpath + file
        fits_data = ppcbfits.ppxf_example_population_gas_sdss(filepath, tie_balmer=True, limit_doublets=True)
        gal_sfrs = get_sfr(filepath, fits_data, velocity_ranges)
        plot_sfr_vel(gal_sfrs, filepath)
        print(gal_sfrs)
        pdf.savefig()
        plt.close('all')



print('-----------Finished-----------!')

