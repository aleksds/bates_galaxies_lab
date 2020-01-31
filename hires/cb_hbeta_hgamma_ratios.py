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


def plot_area_of_interestAAS(plot_data,plottext,filepath,ppxf_data=False):

    #Get relevant stuff for generating a figure
    em_lowlim,em_maxlim,cont_low,cont_high,flux,lam,cont_const = plot_data[0],plot_data[1],plot_data[2],plot_data[3],\
                                                                 plot_data[4],plot_data[5],plot_data[6]
    lam_c = lam[cont_low:cont_high+1]
    flux_c = flux[cont_low:cont_high+1]
    plot_fulln = cbfits.get_plot_title(filepath)
    plot_title = plot_fulln[21:26]
    typetag = plottext[0]
    fluxvalue = plottext[1]

    fluxdata = cont_const[2]

    #Generating that figure
    #with PdfPages('Emissionline.pdf') as pdf:
    fig = plt.figure()
    fig.suptitle(plot_title + ' - ' + typetag + ' - ' + fluxvalue)

    ax = fig.add_subplot(1, 2, 1)
    ax.axvspan(lam[em_lowlim], lam[em_maxlim], alpha=0.5, color='cyan', label='Area of Emission')
    # ax.axvspan(lam[em_lowlim], lam[em_maxlim], alpha=0.2, color='red')
    ax.plot(lam_c, flux_c, color='black', linewidth=0.3)

    if ppxf_data:
        galfit_c = cont_const[0][cont_low:cont_high+1]
        gasfit_c = cont_const[1][cont_low:cont_high+1]
        tallstarfit = (galfit_c - gasfit_c) * np.median(fluxdata)
        ax.plot(lam_c, galfit_c, color='orange', linewidth=0.3)
        ax.plot(lam_c, (gasfit_c*np.median(fluxdata))+tallstarfit[0], color='red', linewidth=0.3)
        ax.plot(lam_c, galfit_c-gasfit_c, color='blue', linewidth=0.3)
        ax.plot(lam_c, tallstarfit, color='blue', linewidth=0.3)
        ax.plot(lam_c, fluxdata[cont_low:cont_high+1], color='grey', linewidth=0.3)
    elif not ppxf_data:
        plt.axhline(y=cont_const, color='blue', label='Continuum', alpha=0.3)


    ax.set_xlabel("$\AA ngstr \ddot{o} ms$")
    ax.set_ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")
    plt.legend(loc=2)
    ax.grid(True)
    plt.savefig(plot_title + ' - ' + typetag + '.png')



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
        hg_flux, plot_datag = ppcbfits.pp_get_flux(hg_wav, hg_spread, fits_data)
        hb_flux, plot_datab = ppcbfits.pp_get_flux(hb_wav, hb_spread, fits_data)
        galname = cbfits.get_plot_title(filepath)
        plot_area_of_interestAAS(plot_datab, ['Hbeta', str(hb_flux)], filepath, ppxf_data=True)
        plot_area_of_interestAAS(plot_datag, ['Hgamma', str(hg_flux)], filepath, ppxf_data=True)

        hbfluxes = np.append(hbfluxes, hb_flux)
        hgfluxes = np.append(hgfluxes, hg_flux)
        columns = np.append(columns, galname)
    print('Success')
    return columns, hbfluxes, hgfluxes

# Organize data in pandas, send to csv

index, hbfluxes, hgfluxes = get_ratios()
ratios = hgfluxes/hbfluxes
data = [[hbfluxes[i], hgfluxes[i], ratios[i]] for i in range(0, len(ratios))]
columns = ['H_beta Flux', 'H_gamma Flux', 'Hg/Hb Ratio']
df = pd.DataFrame(data, index, columns)
df.to_csv('Hb-Hg-ratios.csv')
print('-----------Finished-----------!')

