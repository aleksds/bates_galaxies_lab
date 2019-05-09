import importlib
from os import listdir
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


cbfits = importlib.import_module('cb_galaxy_fits_sfr_analysis_tools')

galfolpath = 'C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hires/galaxies_fits/'

fitsfiles = [f for f in listdir(galfolpath)]

with PdfPages('Emissionlines.pdf') as pdf:
    for file in fitsfiles:
        print('\nPlotting file: ', file,'..........')
        filepath = galfolpath+file
        data = cbfits.get_quantities(filepath)
        hb_spread, hb_wav = cbfits.wspread(filepath,'H_beta',400)
        hb_flux, plot_data = cbfits.get_flux(filepath,hb_wav,hb_spread)
        hb_lum = cbfits.get_lum(filepath,hb_flux)
        sfr = cbfits.get_sfr_hbeta(hb_lum)
        cbfits.plot_area_of_interest(plot_data, [hb_flux, hb_lum, sfr],filepath)
        pdf.savefig(bbox_inches="tight")
        plt.close('all')
        print('Success')

print('-----------Finished-----------!')