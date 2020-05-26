import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sedpy import observate
import pysynphot as S

########################################################################################################################
# Preliminary variables and tables
########################################################################################################################

templates = ['Composite1', 'Composite2', 'Composite3', 'AGN1', 'AGN2', 'AGN3', 'AGN4', 'SFG1', 'SFG2', 'SFG3',
             'IR_COLOR1', 'IR_COLOR2', 'IR_COLOR3', 'IR_COLOR4', 'IR_COLOR5', 'IR_COLOR6', 'IR_COLOR7', 'IR_COLOR8']

# ----- Opening table
table = fits.open('hizea_photo_galex_wise_v1.0.fit')

# ----- Setting up variables
cols = ['fuv', 'nuv', 'u', 'g', 'r', 'i', 'z', 'w1', 'w2', 'w3', 'w4',
        'fuv_unc', 'nuv_unc', 'u_unc', 'g_unc', 'r_unc', 'i_unc', 'z_unc', 'w1_unc', 'w2_unc', 'w3_unc', 'w4_unc', 'Z']
filt_waves = table[1].data['PHOT_WAVE'][0].byteswap().newbyteorder()*(10**-4)
gals_redshifts = np.array([[i] for i in table[1].data['Z']])
np.array([[i] for i in table[1].data['Z']])
gal_names = table[1].data['SHORT_NAME'].byteswap().newbyteorder()
# ----- Making a table containing flux data and error - converting FLUX_FLAM to fnu
def flam_to_fnu(flux):
    spec = S.ArraySpectrum(table[1].data['PHOT_WAVE'][0], flux, fluxunits='Flam')
    spec.convert('Fnu')
    return spec.flux

flam = table[1].data['FLUX_FLAM'].byteswap().newbyteorder()
fnu = np.array([flam_to_fnu(flammie) for flammie in flam])
flamu = table[1].data['FLUX_FLAM_ERR'].byteswap().newbyteorder()
fnuu = np.array([flam_to_fnu(flammieu) for flammieu in flamu])
flux_w_err = np.concatenate((fnu, fnuu, gals_redshifts), axis=1)
gals_flux = pd.DataFrame(data=flux_w_err,
                        index=gal_names,
                        columns=cols)

# ----- Making table containing mags data and error
mags = table[1].data['AB_MAG'].byteswap().newbyteorder()
magsu = table[1].data['AB_MAG_ERR'].byteswap().newbyteorder()
mags_w_err = np.concatenate((mags, magsu, gals_redshifts), axis=1)
gals_mag = pd.DataFrame(data=mags_w_err,
                        index=gal_names,
                        columns=cols)

templates = ['Composite1', 'Composite2', 'Composite3', 'AGN1', 'AGN2', 'AGN3', 'AGN4', 'SFG1', 'SFG2', 'SFG3',
             'IR_COLOR1', 'IR_COLOR2', 'IR_COLOR3', 'IR_COLOR4', 'IR_COLOR5', 'IR_COLOR6', 'IR_COLOR7', 'IR_COLOR8']

fits_cols = ['galaxy', 'template_name', 'filter', 'mags', 'wave_eff', 'model_phot', 'model_phot_mags']

colors = {'Composite1': 'silver', 'Composite2': 'rosybrown', 'Composite3': 'darksalmon',
        'AGN1': 'cornflowerblue', 'AGN2': 'blue', 'AGN3': 'slateblue', 'AGN4': 'paleturquoise',
        'SFG1': 'blueviolet', 'SFG2': 'plum', 'SFG3': 'mediumorchid',
        'IR_COLOR1': 'olive', 'IR_COLOR2': 'olivedrab', 'IR_COLOR3': 'yellowgreen', 'IR_COLOR4': 'greenyellow',
        'IR_COLOR5': 'lawngreen', 'IR_COLOR6': 'lightgreen', 'IR_COLOR7': 'darkgreen', 'IR_COLOR8': 'aquamarine'}

########################################################################################################################
# Funtions
########################################################################################################################

# ----- Read the template into a data frame. Name columns something convenient.
def read_template(name):
    temps_path = 'kirkpatrick/'
    temp = pd.read_csv(temps_path+name+'.txt',
                       names=['rest_wavelength','luminosity','DLnu'],
                       skiprows=[0, 1, 2, 3],
                       sep=r'\s{1,}',
                       engine='python')
    return temp

# ----- Function to match a target wavelength's position.
def mask_wave(temp_wavel, target_wave):
    return np.abs(temp_wavel-target_wave) == np.amin(np.abs(temp_wavel-target_wave))

# ----- Designed to take in flux and return mags
def mag(flux):
    mag = -2.5*np.log10(flux)
    return mag

# ----- Fit templates, return column to feed into a new table
def sed_fitting(gal_name, template_name):
    z = gals_mag.loc[gal_name, 'Z']
    template = tempsdf[tempsdf.template_name == template_name]

    # ----- Organizing wavelength and luminosity
    z_temp_wavel = template.rest_wavelength * (1 + z)
    gal_fluxes = gals_flux.loc[gal_name, :][:11].values
    W3_wavelength = filt_waves[9]
    # Figure out where the template lines up with W1
    mask = mask_wave(z_temp_wavel, W3_wavelength)
    # Scale template to match value at W1
    factor = gal_fluxes[9] / float(template.luminosity[mask].values[0])
    luminosity = template.luminosity * factor  # Scale

    # ----- Readying wavelength and flux for sedpy
    wave_aa = np.array(z_temp_wavel[0:-1]) * 1e4
    flux = np.array(luminosity[0:-1])
    fnu = flux * 3631. * 1e-23
    flambda = fnu * 2.998e18 / (wave_aa) ** 2

    # ----- Using sedpy to get wise band photometry based on templates
    filternames = ['wise_w{}'.format(n) for n in ['1', '2', '3', '4']]
    wise_filters = observate.load_filters(filternames)
    model_mags = observate.getSED(wave_aa, flambda, filterlist=wise_filters)
    wave_eff = [f.wave_effective for f in wise_filters]
    model_phot = 10. ** (model_mags / (-2.5))

    rows = pd.DataFrame([[gal_name, template_name, filternames[i],
                          model_mags[i], wave_eff[i], model_phot[i],
                          mag(model_phot[i])] for i in range(len(filternames))],
                        columns=fits_cols)
    return rows

# ----- Plots templates for comparison and calculates lowest chi
def template_comparison(gal_name, template_name):
    z = gals_mag.loc[gal_name, 'Z']
    template = tempsdf[tempsdf.template_name == template_name].reset_index(drop=True)[['rest_wavelength', 'luminosity']]
    z_temp_wavel = template.rest_wavelength * (1 + z)
    gal_fluxes = gals_flux.loc[gal_name, :][:11].values
    W3_wavelength = filt_waves[9]
    # Figure out where the template lines up with W1
    mask = mask_wave(z_temp_wavel, W3_wavelength)
    # Scale template to match value at W1
    factor = gal_fluxes[9]/float(template.luminosity[mask].values[0])
    luminosity = template.luminosity*factor # Scale
    model_phot = sed_fits[(sed_fits.galaxy == gal_name) & (sed_fits.template_name == template_name)].model_phot.array

    gal_unc = gals_flux.iloc[:,11:-1].loc[gal_name].values
    chi = np.sum(np.array([((gal_fluxes[i + 8] - model_phot[i + 1]) / gal_unc[i + 8]) ** 2 for i in range(3)])) / 3

    # Plot
    plot = True

    if plot:
        plot_color = colors[template_name]
        title = gal_name + '-' + template_name

        g = sns.lineplot(x=z_temp_wavel, y=luminosity, color=plot_color, label=template_name, alpha=0.6, ax=ax)
        h = sns.scatterplot(x=filt_waves, y=gal_fluxes, ax=ax, color='blue')
        # ax.scatter(x=np.array(wave_eff)/1e4, y=model_phot, color='red', s=11)
        ax.errorbar(filt_waves, gal_fluxes, yerr=gal_unc, color='blue', ls='none')
        ax.set_ylim([1e-14, 1e-7])
        ax.set_xlim([0.1, 1000.])
        ax.loglog()
        ax.legend()
        ax.set_title(title)
        plt.ioff()
        # plt.savefig(title+'.png')
        # plt.clf()
        # plt.close()

    return (chi)

# ----- Make general plots showing all galaxies (different colors and different filters)
def general_color_plots():
    # Different filters
    fig = plt.figure(figsize=(17, 10))
    ax = fig.add_subplot()
    w2 = sed_fits[sed_fits['filter'] == 'wise_w2'].model_phot_mags.array
    w3 = sed_fits[sed_fits['filter'] == 'wise_w3'].model_phot_mags.array
    w4 = sed_fits[sed_fits['filter'] == 'wise_w4'].model_phot_mags.array
    sns.scatterplot(x=w3 - w4,
                    y=w2 - w3,
                    data=sed_fits.iloc[::4, :].reset_index(), hue='template_name', ax=ax)
    ax.set_ylabel('W3-W4')
    ax.set_xlabel('W2-W3')
    ax.set_title('Color vs color - filters')
    plt.savefig('galex_color_byfilt.png')
    plt.close('all')
    # Different galaxies
    fig = plt.figure(figsize=(17, 10))
    ax = fig.add_subplot()
    w2 = sed_fits[sed_fits['filter'] == 'wise_w2'].model_phot.array
    w3 = sed_fits[sed_fits['filter'] == 'wise_w3'].model_phot.array
    w4 = sed_fits[sed_fits['filter'] == 'wise_w4'].model_phot.array
    sns.scatterplot(x=mag(w3) - mag(w4),
                    y=mag(w2) - mag(w3),
                    data=sed_fits.iloc[::4, :].reset_index(), hue='galaxy', palette='Paired')
    ax.set_ylabel('W3-W4')
    ax.set_xlabel('W2-W3')
    ax.set_title('Color vs color - galaxies')
    plt.savefig('galex_color_bygal.png')
    plt.close('all')

# ----- Making color plots - 12 plots
def color_plots():
    print(' ---------> Making color plots showing different filters...')
    with PdfPages('galex_sed_fitting_colorplots.pdf') as pdf:
        for galaxy in gal_names:
            print("---Making color plot: ", galaxy)
            fig = plt.figure(figsize=(17, 10))
            ax = fig.add_subplot()
            sedf_g = sed_fits.copy()[sed_fits.galaxy == galaxy].reset_index(drop=True)

            w2 = sedf_g[sedf_g['filter'] == 'wise_w2'].model_phot_mags.array
            w3 = sedf_g[sedf_g['filter'] == 'wise_w3'].model_phot_mags.array
            w4 = sedf_g[sedf_g['filter'] == 'wise_w4'].model_phot_mags.array
            w2_gal = gals_mag['w2'][galaxy]
            w3_gal = gals_mag['w3'][galaxy]
            w4_gal = gals_mag['w4'][galaxy]

            sns.scatterplot(x=w3 - w4, y=w2 - w3,
                            data=sedf_g[sedf_g['filter'] == 'wise_w2'], hue='template_name', ax=ax)

            ax.plot(w3_gal - w4_gal, w2_gal - w3_gal, marker='*', markersize=14, label=galaxy)
            ax.set_title(galaxy)
            ax.set_ylabel('W3-W4')
            ax.set_xlabel('W2-W3')
            pdf.savefig(bbox_inches="tight")
            plt.close('all')
        print('-----Finished color plots pdf')

########################################################################################################################
# Running stuff
########################################################################################################################

# ----- Making a table with the templates
print('Reading templates into data frame...')
tempsdf = pd.DataFrame([],columns=['rest_wavelength','luminosity','DLnu'])
for temp_name in templates:
    newdf = read_template(temp_name)
    newdf['template_name'] = [temp_name for i in range(newdf.shape[0])]
    tempsdf = tempsdf.append(newdf)
print('Templates read.')

# ----- Generating fits for galaxies
sed_fits = pd.DataFrame([], columns=fits_cols)
print(' ---------> Fitting data to templates...')
for gal in gal_names:
    print("---Fitting ", gal)
    for tem in templates:
        sed_fits = sed_fits.append(sed_fitting(gal, tem))
print('---Finished fitting data to templates.\n')
sed_fits.reset_index(inplace=True, drop=True)

color_plots()
print('-----Generated galex_sed_fitting_colorplots.pdf with color vs color plost.')
general_color_plots()
print('----Generated galex_color_byfilt.png and galex_color_bygal.png showing color vs color for all gals.')

with PdfPages('galexsed_fitting.pdf') as pdf:
    print('\n ---------> Plotting templates and calculating chi values...\n')
    for gal in gal_names:
        fig = plt.figure(figsize=(17,10))
        ax = fig.add_subplot()
        print("---Fitting ", gal)
        chis = []
        for tem in templates:
            new_chi = template_comparison(gal, tem)
            chis.append(new_chi)
        bestchi_pos = chis.index(min(chis))
        result_string = "Galaxy - " + gal + " - lowest chi tempalte: " + str(templates[bestchi_pos])
        print(result_string)
        plt.text(10, 10**-6.5, result_string, ha='center')
        pdf.savefig(bbox_inches="tight")
        plt.close('all')
    print('Finished!')