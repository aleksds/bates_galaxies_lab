import importlib
from os import listdir
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


#Variables
filepath = 'galaxies_fits/spec-0761-54524-0409.fits'
galfolpath = 'C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hires/galaxies_fits_DR15/'
fitsfiles = [f for f in listdir(galfolpath)]
cbfits = importlib.import_module('cb_galaxy_fits_sfr_analysis_tools')
ppcbfits = importlib.import_module('cb_ppxf_tools')

#Execution sequence


with PdfPages('Emissionlines_PPXF_U_DR15_AAS.pdf') as pdf:
    for file in fitsfiles:
        print('\nPlotting file: ', file, '..........')
        filepath = galfolpath + file
        data = cbfits.get_quantities(filepath,'H_beta')
        hb_spread, hb_wav = cbfits.wspread(filepath, 'H_beta', 400)
        fits_data = ppcbfits.ppxf_example_population_gas_sdss(filepath, tie_balmer=True, limit_doublets=True) # gets [pp,wave,norm_flux,flux,flux&uncert]
        hb_flux, plot_data = ppcbfits.pp_get_flux(hb_wav, hb_spread, fits_data)
        hb_lum = cbfits.get_lum(filepath, hb_flux)
        sfr = cbfits.get_sfr_hbeta(hb_lum)
        cbfits.plot_area_of_interestAAS(plot_data, [hb_flux, hb_lum, sfr], filepath, ppxf_data=True)
        pdf.savefig(bbox_inches="tight")
        plt.close('all')
        print('Success')

print('-----------Finished-----------!')

