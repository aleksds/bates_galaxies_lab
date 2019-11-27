import importlib
from os import listdir
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#Variables
filepath = 'galaxies_fits/spec-0761-54524-0409.fits'
galfolpath = 'C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hires/galaxies_fits_DR15/'
fitsfiles = [f for f in listdir(galfolpath)]
cbfits = importlib.import_module('cb_galaxy_fits_sfr_analysis_tools')
ppcbfits = importlib.import_module('cb_ppxf_tools')

# Main Function

def get_ratios():
    columns = np.array([])
    hbfluxes = np.array([])
    hgfluxes = np.array([])
    for file in fitsfiles:
        print('\nAnalyzing: ', file, '..........')
        filepath = galfolpath + file
        hg_spread, hg_wav = cbfits.wspread(filepath, 'H_gamma', 400)
        hb_spread, hb_wav = cbfits.wspread(filepath, 'H_beta', 400)
        fits_data = ppcbfits.ppxf_example_population_gas_sdss(filepath, tie_balmer=True, limit_doublets=True) # gets [pp,wave,norm_flux,flux,flux&uncert]
        hg_flux, plot_datab = ppcbfits.pp_get_flux(hg_wav, hg_spread, fits_data)
        hb_flux, plot_datag = ppcbfits.pp_get_flux(hb_wav, hb_spread, fits_data)
        galname = cbfits.get_plot_title(filepath)

        hbfluxes = np.append(hbfluxes, hb_flux)
        hgfluxes = np.append(hgfluxes, hg_flux)
        columns = np.append(columns, galname)
    print('Success')
    return columns, hbfluxes, hgfluxes

print('-----------Finished-----------!')

# Organize data in pandas, send to csv

index, hbfluxes, hgfluxes = get_ratios()
ratios = hgfluxes/hbfluxes
data = [[hbfluxes[i], hgfluxes[i], ratios[i]] for i in range(0, len(ratios))]
columns = ['H_beta Flux', 'H_gamma Flux', 'Hg/Hb Ratio']
df = pd.DataFrame(data, index, columns)
df.to_csv('Hb-Hg-ratios.csv')
