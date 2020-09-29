"""
Author: Sofía Edgar
Date: 09/14/20
Purpose: ...
"""

import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sedpy import observate
import pysynphot as S
from uncertainties import ufloat
from numpy import math


# ----- Read the template into a data frame. Name columns something convenient.
def read_template(name):
    temps_path = 'kirkpatrick/'
    temp = pd.read_csv(temps_path + name + '.txt',
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
    mag = -2.5 * np.log10(flux)
    return mag

#----- Function to convert between f_lamda and f_nu
def flam_to_fnu(flux):
    spec = S.ArraySpectrum(table[1].data['PHOT_WAVE'][0], flux, fluxunits='Flam')
    spec.convert('Fnu')
    return spec.flux

"""
----- This is a function that approximates the w4 value by using the uncertainty as an upper limit on 
         the flux and the negative value as the flux."""
def w4_new(replace, galaxy):
    w4_nanomaggies = replace[1].data['w4_nanomaggies'][0]
    w4_nanomaggies_ivar = replace[1].data['w4_nanomaggies_ivar'][0]
    w4_unc = 1. / np.sqrt(w4_nanomaggies_ivar)
    w4_snr = w4_nanomaggies / w4_unc

    # making new numbers
    w4_upper_1sigma = w4_unc
    w4_mag_upper_1sigma = -2.5 * np.log10(w4_upper_1sigma * 1e-9)  # W4 unc as an upper limit on the flux
    w4_mag_perhaps = -2.5 * np.log10(np.abs(w4_nanomaggies) * 1e-9)  # treat the negative value as the flux
    w4_err = w4_mag_perhaps - w4_mag_upper_1sigma  # difference between upper limit and "value

    return w4_mag_perhaps, w4_err


#----- Replaces all the NaN galaxies
def replace_galaxies(gals_flux):
    for galaxy in gal_names:

        if np.isnan(gals_flux.at[galaxy, 'w4']) == True:
            print("---- Replacing galaxy" + galaxy)
            # Assign file name and open file using SDSS-version DR13
            galaxy_file = galaxy + "sdsswise.fits"
            replace = fits.open(galaxy_file)

            # Create table from data
            replace[1].data

            # Find w4 and w4_unc data from file

            new_w4, new_w4_unc = w4_new(replace, galaxy)

            # new_w4 = galaxy_file.data['w4_mag'][0](# These are placeholder values, not the actual replacement data for certain
            # new_w4_unc = w4_new(galaxy)

            print("---- Replaced galaxy" + galaxy)
            # Replace data in gals_flux
            gals_flux.at[galaxy, 'w4'] = new_w4
            gals_flux.at[galaxy, 'w4_unc'] = new_w4_unc

    # What if we run through all the templates here instead of in generating
    # Plots. Will that save time?

#-----fit templates, return a new column to fit into a new table of all 50 galaxies
def sed_fitting(gal_name, template_name):

    z = gals_mag.loc[gal_name, 'Z']  # Get galaxy redshift
    template = tempsdf[tempsdf.template_name == template_name]#extract template using param template_name

    # ----- Organizing wavelength and luminosity
    z_temp_wavel = template.rest_wavelength * (1 + z)  # moves wavelength to redshift
    gal_fluxes = gals_flux.loc[gal_name, :][:11].values
    W2_wavelength = filt_waves[8]

    # Figure out where the template lines up with W2
    mask = mask_wave(z_temp_wavel, W2_wavelength)

    # Scale template to match value at W2
    factor = gal_fluxes[8] / float(template.luminosity[mask].values[0])
    luminosity = template.luminosity * factor  # Scale

    # ----- Readying wavelength and flux for sedpy observate
    #Note, these numbers should not be hard coded in like this..., we need to change it at some point
    wave_aa = np.array(z_temp_wavel[0:-1]) * 1e4  # changes microns to angstroms
    flux = np.array(luminosity[0:-1])
    fnu = flux * 3631. * 1e-23  # converts from maggies to ergs/s/cm^2/Hz<---standard unit people use
    flambda = fnu * 2.998e18 / (wave_aa) ** 2  # see flux.doc in the Drive for math explanation

    # ----- Using sedpy to get wise band photometry based on templates
    filternames = ['wise_w{}'.format(n) for n in ['1', '2', '3', '4']]#Creates 4 element array, wise_Wn
    wise_filters = observate.load_filters(filternames)
    model_mags = observate.getSED(wave_aa, flambda, filterlist=wise_filters)
    wave_eff = [f.wave_effective for f in wise_filters]
    model_phot = 10. ** (model_mags / (-2.5))  # converts from magnitudes to flux??

    # diff = model_phot[1]/float(template.luminosity[mask].values[0])
    #What is diff?
    #
    diff = model_phot[1] / gal_fluxes[8]

    # apply the correction to the model photometry
    model_phot = model_phot / diff

    # how can we caculate the rest-frame 24-micron luminosity
    # apply similar logic to what's done above for the model_phot calculation
    # the difference is that we'll be using rest-wavelength rather than observed wavelength
    # Becca's good idea: we should store this information rather then continually repeating calculations
    wave_rest = np.array(template.rest_wavelength[0:-1])
    wave_rest_aa = wave_rest * 1e4  # microns 10^-6 m to Angstroms 10^-10 m
    flux_me = np.array(luminosity[0:-1]) / diff

    fnu_me = flux_me * 3631. * 1e-23  # converts from "maggies" to erg/s/cm^2/Hz
    flam_me = fnu_me * 2.998e18 / (wave_rest_aa) ** 2
    filter_24 = ['spitzer_mips_24']
    new_filter = observate.load_filters(filter_24)
    model_mag_24 = observate.getSED(wave_rest_aa, flam_me, filterlist=new_filter)
    wave_eff_24 = [f.wave_effective for f in new_filter]
    model_phot24 = 10. ** (model_mag_24 / (-2.5))

    # add chi^2 calculation
    gal_unc = gals_flux.iloc[:, 11:-1].loc[gal_name].values
    chi = np.sum(np.array([((gal_fluxes[i + 8] - model_phot[i + 1]) / gal_unc[i + 8]) ** 2 for i in range(3)])) / 3


    rows = pd.DataFrame([[gal_name, template_name, filternames[i],
                          model_mags[i], wave_eff[i], model_phot[i],
                          mag(model_phot[i]), diff, chi] for i in range(len(filternames))],
                        columns=fits_cols)
    return rows

#Calculates chi for each filter, may not need this later, this is the innermost function
#of the sed_fitting_comparison function. This function calculates all the information for each filter given a template
#and a galaxy.

def calculate_template_chis(gal_fluxes,template,gal_name):

    # create empty lists for information to be stored
    luminosities = []
    uncertainties = []
    #filter_chis = []
    # ----- Organizing wavelength and luminosity
    z = gals_mag.loc[gal_name, 'Z']
    z_temp_wavel = template.rest_wavelength * (1 + z)
    wave_index = 7


    wavelength = filt_waves[wave_index]  # Assigns
    # Figure out where the template lines up with W2
    mask = mask_wave(z_temp_wavel, wavelength)

    # Scale template to match value at W2
    factor = gal_fluxes[wave_index] / float(template.luminosity[mask].values[0])
    luminosity = template.luminosity * factor  # Scale

    # ----- Readying wavelength and flux for sedpy
    wave_aa = np.array(z_temp_wavel[0:-1]) * 1e4  # changes angstroms to microns
    flux = np.array(luminosity[0:-1])
    fnu = flux * 3631. * 1e-23  # converts from angstroms to microns?
    flambda = fnu * 2.998e18 / (wave_aa) ** 2  # see flux.doc in the Drive for math explanation

    # ----- Using sedpy to get wise band photometry based on templates
    filternames = ['wise_w{}'.format(n) for n in ['1', '2', '3', '4']]
    wise_filters = observate.load_filters(filternames)
    model_mags = observate.getSED(wave_aa, flambda, filterlist=wise_filters)
    wave_eff = [f.wave_effective for f in wise_filters]
    model_phot = 10. ** (model_mags / (-2.5))  # converts from magnitudes to flux??

    diff = model_phot[1] / gal_fluxes[wave_index]

    # apply the correction to the model photometry
    model_phot = model_phot / diff

    # ------ Calculate Chi
    gal_unc = gals_flux.iloc[:, 11:-1].loc[gal_name].values

    # add chi^2 calculation
    chi = np.sum(np.array(
        [((gal_fluxes[i + 8] - model_phot[i + 1]) / gal_unc[i + 8]) ** 2 for i in range(3)])) / 3
    filter_chis.append(chi)
    wave_index += 1
    luminosities.append(luminosity)
    uncertainties.append(gal_unc)

    return [filter_chis, luminosities, uncertainties,z_temp_wavel]

#For each template extract best chi value, called in sed_fitting_comparison, part of refactorization
#This function extracts all of the optimal information for each template in the set of 19
def template_fitting_comparison(gal_name,template_name):
    flux_24 = []
    best_chis = []
    template_info = []
    # This is supposed to fix the reassignment of 'rows' every time, not sure how it was affecting the end result
    #Next thing is to figure out what rows does and if we need to store each row_element, which is galaxy to template specific
    rows = []



    template = tempsdf[tempsdf.template_name == template_name]

    # print("for template "+ template_name+" z_temp_wave1 is:"+str(z_temp_wave1))
    gal_fluxes = gals_flux.loc[gal_name, :][:11].values

    # Call the template_chis_comparison to

    filter_chis_info = template_fitting_comparison(gal_fluxes, template, gal_name)
    filter_chis = filter_chis_info[0]
    luminosities = filter_chis_info[1]
    uncertainties = filter_chis_info[2]
    z_temp_wavel = filter_chis_info[3]

    chi = min(filter_chis)  # finds best chi of the four bandpass filters

    chi_index = filter_chis.index(chi)

    best_luminosity = luminosities[chi_index]
    best_gal_unc = uncertainties[chi_index]
    plotting_inf = [z_temp_wavel, best_luminosity, gal_fluxes, best_gal_unc, chi_index]

    best_chis.append(chi)  # Adds best template chi to list of template chis
    template_info.append(plotting_inf)  # Adds info for best temp chi to list of accompanying  info for best chis

    # how can we caculate the rest-frame 24-micron luminosity
    # apply similar logic to what's done above for the model_phot calculation
    # the difference is that we'll be using rest-wavelength rather than observed wavelength
    # Becca's good idea: we should store this information rather then continually repeating calculations
    wave_rest = np.array(template.rest_wavelength[0:-1])
    wave_rest_aa = wave_rest * 1e4  # microns 10^-6 m to Angstroms 10^-10 m
    flux_me = np.array(luminosity[0:-1]) / diff

    fnu_me = flux_me * 3631. * 1e-23  # converts from "maggies" to erg/s/cm^2/Hz
    flam_me = fnu_me * 2.998e18 / (wave_rest_aa) ** 2  #

    filter_24 = observate.load_filters(['spitzer_mips_24'])

    model_mag_24 = observate.getSED(wave_rest_aa, flam_me, filterlist=filter_24)
    wave_eff_24 = [f.wave_effective for f in filter_24]
    model_phot24 = 10. ** (model_mag_24 / (-2.5))

    flux_24.append(model_phot24)

    row_element = pd.DataFrame([[gal_name, template_name, filternames[i],
                          model_mags[i], wave_eff[i], model_phot[i],
                          mag(model_phot[i]), diff, chi]
                          for i in range(len(filternames))],
                          columns=fits_cols)
    rows.append[row_element]
    return [rows, best_chis,template_info]

'''
def sed_fitting_comparison(gal_name):

    for template_name in templates: # Goes through each of the templates
        for filt in filt_waves[7:]: #Goes through each of the filters
           chi_information_for_a_galaxy =  template_fitting_comparison(gal_name, template_name)

    flux_24_array = np.array(model_phot24)

    mean = np.mean(flux_24_array)
    std = np.std(flux_24_array)
    minimum = flux_24_array.min()
    maximum = flux_24_array.max()
    flux_stats = [mean, std, minimum, maximum]
    # print(mean,std,mean/std, maximum, minimum)


#Why do I return a list of best_chis and template info instead of finding the minimum chi here?
    return [rows, best_chis, template_info, flux_stats]
'''

# What if we run through all the templates here instead of in generating
# Plots. Will that save time?

#Original sed_fitting_comparison
'''
def sed_fitting_comparison(gal_name):
    flux_24 = []
    best_chis = []
    template_info = []
    for template_name in templates:

        z = gals_mag.loc[gal_name, 'Z']
        template = tempsdf[tempsdf.template_name == template_name]
        filter_chis = []
        filter_info = []
        luminosities = []
        uncertainties = []
        # ----- Organizing wavelength and luminosity
        z_temp_wavel = template.rest_wavelength * (1 + z)

        # print("for template "+ template_name+" z_temp_wave1 is:"+str(z_temp_wave1))
        gal_fluxes = gals_flux.loc[gal_name, :][:11].values
        wave_index = 7

        for filt in filt_waves[7:]:
            wavelength = filt_waves[wave_index]  # Assigns
            # Figure out where the template lines up with W2
            mask = mask_wave(z_temp_wavel, wavelength)

            # Scale template to match value at W2
            factor = gal_fluxes[wave_index] / float(template.luminosity[mask].values[0])
            luminosity = template.luminosity * factor  # Scale

            # ----- Readying wavelength and flux for sedpy
            wave_aa = np.array(z_temp_wavel[0:-1]) * 1e4  # changes angstroms to microns
            flux = np.array(luminosity[0:-1])
            fnu = flux * 3631. * 1e-23  # converts from angstroms to microns?
            flambda = fnu * 2.998e18 / (wave_aa) ** 2  # see flux.doc in the Drive for math explanation

            # ----- Using sedpy to get wise band photometry based on templates
            filternames = ['wise_w{}'.format(n) for n in ['1', '2', '3', '4']]
            wise_filters = observate.load_filters(filternames)
            model_mags = observate.getSED(wave_aa, flambda, filterlist=wise_filters)
            wave_eff = [f.wave_effective for f in wise_filters]
            model_phot = 10. ** (model_mags / (-2.5))  # converts from magnitudes to flux??

            # diff = model_phot[1]/float(template.luminosity[mask].values[0])
            diff = model_phot[1] / gal_fluxes[wave_index]

            # apply the correction to the model photometry
            model_phot = model_phot / diff

            # ------ Calculate Chi

            gal_unc = gals_flux.iloc[:, 11:-1].loc[gal_name].values
            # print(gal_fluxes)
            # print(model_phot)
            # print(gal_unc)
            # add chi^2 calculation
            chi = np.sum(
                np.array([((gal_fluxes[i + 8] - model_phot[i + 1]) / gal_unc[i + 8]) ** 2 for i in range(3)])) / 3
            filter_chis.append(chi)
            wave_index += 1
            luminosities.append(luminosity)
            uncertainties.append(gal_unc)
            # info = [z_temp_wavel,luminosity,gal_fluxes]
            # filter_info.append(info)

        # print("filter chis" + str(filter_chis) )

        chi = min(filter_chis)  # finds best chi of the four bandpass filters

        chi_index = filter_chis.index(chi)

        best_luminosity = luminosities[chi_index]
        best_gal_unc = uncertainties[chi_index]
        # plotting_inf = filter_info[chi_index]# Gets the accompanying info for best filter chi
        plotting_inf = [z_temp_wavel, best_luminosity, gal_fluxes, best_gal_unc, chi_index]

        best_chis.append(chi)  # Adds best template chi to list of template chis
        template_info.append(plotting_inf)  # Adds info for best temp chi to list of accompanying  info for best chis

        # how can we caculate the rest-frame 24-micron luminosity
        # apply similar logic to what's done above for the model_phot calculation
        # the difference is that we'll be using rest-wavelength rather than observed wavelength
        # Becca's good idea: we should store this information rather then continually repeating calculations
        wave_rest = np.array(template.rest_wavelength[0:-1])
        wave_rest_aa = wave_rest * 1e4  # microns 10^-6 m to Angstroms 10^-10 m
        flux_me = np.array(luminosity[0:-1]) / diff

        fnu_me = flux_me * 3631. * 1e-23  # converts from "maggies" to erg/s/cm^2/Hz
        flam_me = fnu_me * 2.998e18 / (wave_rest_aa) ** 2  #

        filter_24 = observate.load_filters(['spitzer_mips_24'])

        model_mag_24 = observate.getSED(wave_rest_aa, flam_me, filterlist=filter_24)
        wave_eff_24 = [f.wave_effective for f in filter_24]
        model_phot24 = 10. ** (model_mag_24 / (-2.5))

        flux_24.append(model_phot24)

        rows = pd.DataFrame([[gal_name, template_name, filternames[i],
                              model_mags[i], wave_eff[i], model_phot[i],
                              mag(model_phot[i]), diff, chi]
                             for i in range(len(filternames))],
                            columns=fits_cols)
        print(rows)

    flux_24_array = np.array(flux_24)

    mean = np.mean(flux_24_array)
    std = np.std(flux_24_array)
    minimum = flux_24_array.min()
    maximum = flux_24_array.max()
    flux_stats = [mean, std, minimum, maximum]

    return [rows, best_chis, template_info, flux_stats]
'''
#-----This function fits each template to each of the 50 galaxies and returns the template with the smallest chi value

#Plots templates for comparison and calculate the lowest chi

def template_comparison(gal_name, template_name):
    z = gals_mag.loc[gal_name, 'Z']
    template = tempsdf[tempsdf.template_name == template_name].reset_index(drop=True)[['rest_wavelength', 'luminosity']]
    z_temp_wavel = template.rest_wavelength * (1 + z)
    #print("z_temp_wavel" + str(z_temp_wavel))
    gal_fluxes = gals_flux.loc[gal_name, :][:11].values
    W2_wavelength = filt_waves[8]

    # Figure out where the template lines up with W2
    mask = mask_wave(z_temp_wavel, W2_wavelength)

    # Scale template to match value at W2
    factor = gal_fluxes[8] / float(template.luminosity[mask].values[0])
    luminosity = template.luminosity * factor

    # Scale and match luminosity of rest wavelength with template
    model_phot = sed_fits[(sed_fits.galaxy == gal_name) & (sed_fits.template_name == template_name)].model_phot.array

    gal_unc = gals_flux.iloc[:, 11:-1].loc[gal_name].values
    chi = np.sum(np.array([((gal_fluxes[i + 8] - model_phot[i + 1]) / gal_unc[i + 8]) ** 2 for i in range(3)])) / 3

    # After this the function can stop right? What is the rest for? Why do we plot each template?
    # Plot
    return (chi)

def plot_template_comparison(gal_name,template_name,luminosity,gal_fluxes,z_temp_wavel,gal_unc):

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
    # plt.close() '''


'''
def template_comparison(gal_name, template_name):
    z = gals_mag.loc[gal_name, 'Z']
    template = tempsdf[tempsdf.template_name == template_name].reset_index(drop=True)[['rest_wavelength', 'luminosity']]
    z_temp_wavel = template.rest_wavelength * (1 + z)
    print("z_temp_wavel" + str(z_temp_wavel))
    gal_fluxes = gals_flux.loc[gal_name, :][:11].values

    W2_wavelength = filt_waves[8]
    # Figure out where the template lines up with W2
    mask = mask_wave(z_temp_wavel, W2_wavelength)
    # Scale template to match value at W2
    factor = gal_fluxes[8] / float(template.luminosity[mask].values[0])
    luminosity = template.luminosity * factor
    # Scale and match luminosity of rest wavelength with template
    model_phot = sed_fits[(sed_fits.galaxy == gal_name) & (sed_fits.template_name == template_name)].model_phot.array

    gal_unc = gals_flux.iloc[:, 11:-1].loc[gal_name].values
    chi = np.sum(np.array([((gal_fluxes[i + 8] - model_phot[i + 1]) / gal_unc[i + 8]) ** 2 for i in range(3)])) / 3

    # After this the function can stop right? What is the rest for? Why do we plot each template?
    # Plot



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

    return (chi)'''

# Whichever template has the smallest chi value is assigned to the galaxy
def template_comparison_optimal(gal_name):


        returned_values = sed_fitting_comparison(gal_name)  # creates list/array? of values returned from function

        chis = returned_values[1]  # Extracts the list of chi values from returned chis-list element in list

        plotting_values = returned_values[2]  # Extracts the accompanying info for best chis

        minimum = min(chis)
        bestchi_pos = chis.index(minimum)

        true_plot_info = plotting_values[bestchi_pos]  # Creates list of associated all info for the best chi value
        # print("minimum chi is" +str(minimum)+ "at index "+str(bestchi_pos))
        # print("true plot info is: "+ str(true_plot_info)+" and length is "+str(len(true_plot_info)))

        # Note that the values of true_plot_info are assigned as follows
        # z_temp_wave1, luminosity, gal_fluxes, gal_unc,wave_index

        # Plot
        plot = True

        if plot:
            plot_color = colors[str(templates[bestchi_pos])]
            title = gal + '-' + str(templates[bestchi_pos]) + ' using W' + str(true_plot_info[4] + 1)

            g = sns.lineplot(x=true_plot_info[0], y=true_plot_info[1], color=plot_color,
                             label=str(templates[bestchi_pos]), alpha=0.6, ax=ax)
            h = sns.scatterplot(x=filt_waves, y=true_plot_info[2], ax=ax, color='blue')
            # ax.scatter(x=np.array(wave_eff)/1e4, y=model_phot, color='red', s=11)
            ax.errorbar(filt_waves, true_plot_info[2], yerr=true_plot_info[3], color='blue', ls='none')
            ax.set_ylim([1e-14, 1e-7])
            ax.set_xlim([0.1, 1000.])
            ax.loglog()
            ax.legend()
            ax.set_title(title)
            plt.ioff()

        return (minimum)
    # Whichever template has the smallest chi value is assigned to the galaxy

"""-----This function ..."""
def template_comparison_full(gal_name):
    returned_values = sed_fitting_comparison(gal_name)  # creates list/array? of values returned from function

    chis = returned_values[1]  # Extracts the list of chi values from returned chis-list element in list

    plotting_values = returned_values[2]

    for template_name in templates:
        # Extracts the accompanying info for best chis
        fig = plt.figure(figsize=(17, 10))
        ax = fig.add_subplot(1, 1, 1)
        this_index = templates.index(template_name)
        temp_chi = chis[this_index]
        minimum = min(chis)
        # bestchi_pos = chis.index(minimum)

        true_plot_info = plotting_values[this_index]  # Creates list of associated all info for the best chi value
        # print("minimum chi is" +str(minimum)+ "at index "+str(bestchi_pos))
        # print("true plot info is: "+ str(true_plot_info)+" and length is "+str(len(true_plot_info)))

        # Note that the values of true_plot_info are assigned as follows
        # z_temp_wave1, luminosity, gal_fluxes, gal_unc,wave_index

        # Plot
        plot = True

        if plot:
            plot_color = colors[template_name]
            title = gal + '-' + template_name + ' using W' + str(true_plot_info[4] + 1)

            g = sns.lineplot(x = true_plot_info[0], y = true_plot_info[1], color = plot_color, label= template_name, alpha = 0.6,
                             ax = ax)
            h = sns.scatterplot(x = filt_waves, y = true_plot_info[2], ax=ax, color='blue')
            # ax.scatter(x=np.array(wave_eff)/1e4, y=model_phot, color='red', s=11)
            ax.errorbar(filt_waves, true_plot_info[2], yerr = true_plot_info[3], color = 'blue', ls ='none')
            ax.set_ylim([1e-14, 1e-7])
            ax.set_xlim([0.1, 1000.])
            ax.loglog()
            ax.legend()
            ax.set_title(title)
            plt.ioff()

            result_string = "Galaxy - " + gal + " - template: " + template_name + " lowest chi:" + str(temp_chi)
            #print(result_string)
            plt.text(10, 10 ** -6.5, result_string, ha='center')
            pdf.savefig(bbox_inches ="tight")
            plt.close('all')

    return (minimum)


# Whichever template has the smallest chi value is assigned to the galaxy

"""
----- This function makes general color plots showing all the galaxies 
----- with different colors and different filters"""
def general_color_plots():
    # Different filters
    fig = plt.figure(figsize=(17, 10))
    ax = fig.add_subplot(1,1,1)
    w2 = sed_fits[sed_fits['filter'] == 'wise_w2'].model_phot_mags.array
    w3 = sed_fits[sed_fits['filter'] == 'wise_w3'].model_phot_mags.array
    w4 = sed_fits[sed_fits['filter'] == 'wise_w4'].model_phot_mags.array
    sns.scatterplot(x=w3 - w4,
                    y=w2 - w3,
                    data=sed_fits.iloc[::4, :].reset_index(), hue='template_name', ax=ax)
    ax.set_ylabel('W3-W4')
    ax.set_xlabel('W2-W3')
    ax.set_title('Color vs color - filters')
    plt.savefig('se_galex_color_byfilt.png')
    plt.close('all')
    # Different galaxies
    fig = plt.figure(figsize=(17, 10))
    ax = fig.add_subplot(1,1,1)
    w2 = sed_fits[sed_fits['filter'] == 'wise_w2'].model_phot.array
    w3 = sed_fits[sed_fits['filter'] == 'wise_w3'].model_phot.array
    w4 = sed_fits[sed_fits['filter'] == 'wise_w4'].model_phot.array
    sns.scatterplot(x=mag(w3) - mag(w4),
                    y=mag(w2) - mag(w3),
                    data=sed_fits.iloc[::4, :].reset_index(), hue='galaxy', palette='Paired')
    ax.set_ylabel('W3-W4')
    ax.set_xlabel('W2-W3')
    ax.set_title('Color vs color - galaxies')
    plt.savefig('se_galex_color_bygal.png')
    plt.close('all')


"""
----- This function makes color v. color plots for each individual galaxy in the data set
"""
def color_plots():
    print(' ---------> Making color plots showing different filters...')
    with PdfPages('se_galex_sed_fitting_colorplots.pdf') as pdf:
        for galaxy in gal_names:
            print("---Making color plot: ", galaxy)
            fig = plt.figure(figsize=(17, 10))
            ax = fig.add_subplot(1, 1, 1)
            sedf_g = sed_fits.copy()[sed_fits.galaxy == galaxy].reset_index(drop=True)

            w2 = sedf_g[sedf_g['filter'] == 'wise_w2'].model_phot_mags.array
            w3 = sedf_g[sedf_g['filter'] == 'wise_w3'].model_phot_mags.array
            w4 = sedf_g[sedf_g['filter'] == 'wise_w4'].model_phot_mags.array
            w2_gal = gals_mag['w2'][galaxy]
            w3_gal = gals_mag['w3'][galaxy]
            w4_gal = gals_mag['w4'][galaxy]

            w2_unc = ufloat(gals_mag['w2'][galaxy], gals_mag['w2_unc'][galaxy])
            w3_unc = ufloat(gals_mag['w3'][galaxy], gals_mag['w3_unc'][galaxy])
            w4_unc = ufloat(gals_mag['w4'][galaxy], gals_mag['w4_unc'][galaxy])

            sns.scatterplot(x=w3 - w4, y=w2 - w3,
                            data=sedf_g[sedf_g['filter'] == 'wise_w2'], hue='template_name', ax=ax)

            ax.plot(w3_gal - w4_gal, w2_gal - w3_gal, marker='*', markersize=14, label=galaxy)
            ax.errorbar((w3_unc - w4_unc).n, (w2_unc - w3_unc).n, xerr=(w3_unc - w4_unc).s, yerr=(w2_unc - w3_unc).s,
                        color='blue', ls='none')
            # plots both vertical and horizontal errorbars
            ax.set_title(galaxy)
            ax.set_ylabel('W2-W3')
            ax.set_xlabel('W3-W4')
            pdf.savefig(bbox_inches="tight")
            plt.close('all')
        print('-----Finished color plots pdf')


#def main():
#---------------------------------------------------------------------------------------------------------------------
#Non-function code

# Declaring Variables
templates = ['Composite1', 'Composite2', 'Composite3', 'Composite4', 'AGN1', 'AGN2', 'AGN3', 'AGN4', 'SFG1', 'SFG2',
                 'SFG3',
                 'IR_COLOR1', 'IR_COLOR2', 'IR_COLOR3', 'IR_COLOR4', 'IR_COLOR5', 'IR_COLOR6', 'IR_COLOR7', 'IR_COLOR8']


# This is a combination of the Libraries for High Z Dusty Galaxies and AGN
# It uses all templates from the Color-Based Library and all templates from Comprehensive Library


fits_cols = ['galaxy', 'template_name', 'filter', 'mags', 'wave_eff', 'model_phot', 'model_phot_mags', 'diff', 'chi']

colors = {'Composite1': 'silver', 'Composite2': 'rosybrown', 'Composite3': 'darksalmon', 'Composite4': 'deeppink',
          'AGN1': 'cornflowerblue', 'AGN2': 'blue', 'AGN3': 'slateblue', 'AGN4': 'paleturquoise',
          'SFG1': 'blueviolet', 'SFG2': 'plum', 'SFG3': 'mediumorchid',
          'IR_COLOR1': 'olive', 'IR_COLOR2': 'olivedrab', 'IR_COLOR3': 'yellowgreen', 'IR_COLOR4': 'greenyellow',
          'IR_COLOR5': 'lawngreen', 'IR_COLOR6': 'lightgreen', 'IR_COLOR7': 'darkgreen', 'IR_COLOR8': 'aquamarine'}

# Open the table
table = fits.open('hizea_photo_galex_wise_v1.0.fit')

cols = ['fuv', 'nuv', 'u', 'g', 'r', 'i', 'z', 'w1', 'w2', 'w3', 'w4',
        'fuv_unc', 'nuv_unc', 'u_unc', 'g_unc', 'r_unc', 'i_unc', 'z_unc', 'w1_unc', 'w2_unc', 'w3_unc', 'w4_unc', 'Z']
filt_waves = table[1].data['PHOT_WAVE'][0].byteswap().newbyteorder()*(10**-4)

# PHOT_WAVE finds central wavelength of each photometric data point given in angstroms

gals_redshifts = np.array([[i] for i in table[1].data['Z']])
np.array([[i] for i in table[1].data['Z']]) # What does this line actually do?
gal_names = table[1].data['SHORT_NAME'].byteswap().newbyteorder()


print("Success we made it through the imports!")

# Make a table with fluxes and errors
flam = table[1].data['FLUX_FLAM'].byteswap().newbyteorder()
fnu = np.array([flam_to_fnu(flammie) for flammie in flam])
flamu = table[1].data['FLUX_FLAM_ERR'].byteswap().newbyteorder()
fnuu = np.array([flam_to_fnu(flammieu) for flammieu in flamu])
flux_w_err = np.concatenate((fnu, fnuu, gals_redshifts), axis=1) # Creates an array of fluxes and errors from smaller arrays
gals_flux = pd.DataFrame(data=flux_w_err,# Creates pandas data frame from array
                    index=gal_names,
                    columns=cols)
replace_galaxies(gals_flux)
# Really it is unecessary to write a function. This can just be done in a script w/o it



mags = table[1].data['AB_MAG'].byteswap().newbyteorder()
magsu = table[1].data['AB_MAG_ERR'].byteswap().newbyteorder()
mags_w_err = np.concatenate((mags, magsu, gals_redshifts), axis=1)
gals_mag = pd.DataFrame(data=mags_w_err,
                        index=gal_names,
                        columns=cols)

# Making a table with the templates
print('Reading templates into data frame...')
tempsdf = pd.DataFrame([], columns=['rest_wavelength', 'luminosity', 'DLnu'])

for temp_name in templates:
    newdf = read_template(temp_name)
    newdf['template_name'] = [temp_name for i in range(newdf.shape[0])]
    print("Read " + temp_name)
    tempsdf = tempsdf.append(newdf)
print('Templates read.')
# tempsdf

#Generating fits for galaxies
#Full galaxy sedfit
sed_fits = pd.DataFrame([], columns=fits_cols)
print(' ---------> Fitting templates to data...')
for gal in gal_names:  # [0:5]
    print("---Fitting ", gal)
    for tem in templates:
        sed_fits = sed_fits.append(sed_fitting(gal, tem))
print('---Finished fitting templates to data.\n')
sed_fits.reset_index(inplace=True, drop=True)'''


#optimized sedfit
sed_fits = pd.DataFrame([], columns=fits_cols)
print(' ---------> Fitting templates to data...')
for gal in gal_names:
    print("---Fitting ", gal)
    output_list = sed_fitting_comparison(gal)
    extract_rows = output_list[0]
    sed_fits = sed_fits.append(extract_rows)

print('---Finished fitting templates to data.\n')
sed_fits.reset_index(inplace = True, drop = True)

#Make color v color plots
color_plots()
print('-----Generated se_galex_sed_fitting_colorplots.pdf with color vs color plots.')

general_color_plots()
print('----Generated se_galex_color_byfilt.png and se_galex_color_bygal.png showing color vs color for all gals.')


#PRODUCE A PDF FILE SHOWING TEMPLATE FITS

#Plots the best fit templates
# this takes some time to run
with PdfPages('se_galexsed_fitting_optimize.pdf') as pdf:  # prints out best fit
    print('\n ---------> Plotting templates and calculating chi values...\n')
    for gal in gal_names:
        fig = plt.figure(figsize=(17, 10))
        ax = fig.add_subplot(1, 1, 1)
        print("---Fitting ", gal)
        best_fit = template_comparison_optimal(gal)

        result_string = "Galaxy - " + gal + " - lowest chi: " + str(best_fit)
        #print(result_string)
        plt.text(10, 10 ** -6.5, result_string, ha='center')
        pdf.savefig(bbox_inches="tight")
        plt.close('all')
    print('Finished!')

#Plots all 19 templates for 1 galaxy on a separate plot
#This takes some time to run (about a minute per galaxy?)
'''with PdfPages('se_galexsed_fitting_individuals.pdf') as pdf:  # prints out all templates
    print('\n ---------> Plotting templates and calculating chi values...\n')
    for gal in gal_names:  # [0:5]
        fig = plt.figure(figsize=(17, 10))
        ax = fig.add_subplot(1, 1, 1)
        result_string = "Galaxy - " + gal + " - lowest chi template: " + str(templates[bestchi_pos])
        print("---Fitting ", gal, "w/ all templates")
        best_chi = template_comparison_full(gal)
        print("Smallest chi:", best_chi)
        plt.text(10, 10 ** -6.5, result_string, ha='center')
        pdf.savefig(bbox_inches="tight")
        plt.close('all')

    print('Finished!')'''

# Plots all 19 templates for 1 galaxy on a single plot
# this takes some time to run (about a minute per galaxy?)
'''with PdfPages('se_galexsed_fitting_full.pdf') as pdf:  # prints out all templates
    print('\n ---------> Plotting templates and calculating chi values...\n')
    for gal in gal_names[0:3]:  # [0:5]
        fig = plt.figure(figsize=(17, 10))
        ax = fig.add_subplot(1, 1, 1)
        print("---Fitting ", gal)
        chis = []
        for tem in templates:
            new_chi = template_comparison(gal, tem)
            chis.append(new_chi)
        bestchi_pos = chis.index(min(chis))
        result_string = "Galaxy - " + gal + " - lowest chi template: " + str(templates[bestchi_pos])
        print(result_string)
        plt.text(10, 10 ** -6.5, result_string, ha ='center')
        pdf.savefig(bbox_inches="tight")
        plt.close('all')
    print('Finished!')'''


#main()
