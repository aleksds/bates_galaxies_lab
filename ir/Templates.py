#
# Author: Grayson Petter
# Templates.py
# Functions to redshift, interpolate, and integrate IR templates in order to simulate observed colors & luminosities

import numpy as np
from astropy import units as u
import math
import pickle
import WISE
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
import glob
import importlib
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii


importlib.reload(WISE)
# reload(WISE)

# speed of light
c = 299792458.  # m/s

# path to project

projpath = os.getcwd() + '/'  # '../'#'/Users/graysonpetter/Desktop/IRSFRs/'
print("projpath:", projpath)
# set cosmology for calculating luminosity distance
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Read in all templates
templates = glob.glob(projpath + 'Comprehensive_library/SFG*.txt')
templates.extend(
    glob.glob(projpath + 'Comprehensive_library/Comp*.txt'))
templates.extend(
    glob.glob(projpath + 'Comprehensive_library/AGN*.txt'))

print('templates:', templates)

# read in the W3 & W4 bandpasses
wise_bandpasses_3_4 = sorted(glob.glob(projpath + 'bandpass/*.txt'))[2:4]


# ("wise_bandpasses1",wise_bandpasses_3_4)

# redshift a template SED (wavelengths only)
def redshift_spectrum(z, template, trim):
    # t = pd.read_csv(template, delim_whitespace=True, engine='python', header=None, skiprows=3)
    t = pd.read_csv(template, delim_whitespace=True, engine='python', header=None, skiprows=9)

    # read off wavelengths and luminosities from template
    wavelengths = np.array(t.iloc[:, 0])
    wavelengths = wavelengths.astype(float)
    # print('wavelengths: ', wavelengths)
    Lums = np.array(t.iloc[:, 1])

    # cut template down to 8-1000 microns (TIR) definition
    if trim:
        spec_range = np.where((wavelengths >= 8.) & (wavelengths <= 1000.))[0]
        wavelengths = wavelengths[spec_range]
        Lums = Lums[spec_range]

    # redshift wavelengths
    shifted_len = np.array(wavelengths) * (1 + z)

    # get luminosity at 12 & 22 micron in observed frame
    twelve_mu = (np.abs(shifted_len - 12)).argmin()
    twenty_two_mu = (np.abs(shifted_len - 22)).argmin()

    return wavelengths, Lums, Lums[twelve_mu], Lums[twenty_two_mu], shifted_len


# linearly interpolate the SED in frequency space to make integration simple
def interpolate_spec(shifted_spec, model):
    # convert wavelengths in microns to frequencies in Hz, 10**6 converts microns to meters
    nus = (10 ** 6) * c / (shifted_spec[0])

    # reverse lists so frequencies go from low to high for simplicity
    reversed_nus = np.flipud(nus).flatten()
    # also reverse luminosities
    reversed_lums = np.flipud(shifted_spec[1])

    # calculate constant frequency interval to interpolate on
    if model:
        step = reversed_nus[1] - reversed_nus[0]
        dx = round(step, -(len(str(int(step))) - 1))
    else:
        dx = 10000000000

    # find smallest factor of dx Hz greater than the smallest frequency in the list to start the interpolation
    start = (reversed_nus[0] + int(dx)) - (reversed_nus[0] % int(dx))

    # range of frequency across entire template
    span = reversed_nus[len(reversed_nus) - 1] - reversed_nus[0]
    # number of frequency intervals to interpolate on
    chunks = int(math.floor(span / dx))

    # lists for interpolated values
    new_nus, new_lums = [], []
    current_nu = start

    # linearly interpolate to frequencies in dx Hz steps
    for x in range(chunks):
        new_nus.append(current_nu)
        new_lums.append(np.interp(current_nu, reversed_nus, reversed_lums))
        current_nu += dx

    return new_nus, new_lums


# integrate spectrum using trapezoidal method
def integrate_spectrum(freqs, Ls):
    return np.trapz(y=Ls, x=freqs)


def simulate_wise_fluxes_for_colors(z, tems, bands, csv):
    tot_mag_list, template_names = [], []
    # iterate through templates
    for tem in tems:
        # redshift template
        red_spec = redshift_spectrum(z, tem, False)
        red_waves = np.array(red_spec[4])
        lumi = np.array(red_spec[1])

        normalized = []

        # iterate through WISE bands
        for y in range(len(bands)):
            if csv:
                band = pd.read_csv(bands[y], header=None, engine='python')
            else:
                band = pd.read_csv(bands[y], header=None, delim_whitespace=True, engine='python')
            bandwaves = np.array(band.iloc[:, 0])
            band_response = np.array(band.iloc[:, 1])

            # trim template to same wavelength range as WISE band
            cut = np.where((red_waves >= np.min(bandwaves)) & (red_waves <= np.max(bandwaves)))[0]
            trimmed_y = red_waves[cut]
            trimmed_L = lumi[cut]

            # interpolate template to band wavelengths, multiply by the response at that wavelength
            inter_lum = []
            for j in range(len(bandwaves)):
                inter_lum.append(band_response[j] * (np.interp(bandwaves[j], trimmed_y, trimmed_L)))

            # crude method
            """sum_lum = np.sum(np.array(inter_lum))
            sum_waves = np.sum(np.array(band_response))
            normalized.append(sum_lum/sum_waves)"""

            # integrate template multiplied by response function
            spectrum = [bandwaves, inter_lum]
            interped_again = interpolate_spec(spectrum, True)
            wise_lums = integrate_spectrum(interped_again[0], interped_again[1])

            # integrate wise band
            band_spectrum = [bandwaves, band_response]
            interped_band = interpolate_spec(band_spectrum, True)
            integrated_band = integrate_spectrum(interped_band[0], interped_band[1])
            # divide two
            normalized.append(wise_lums / integrated_band)

        tot_mag_list.append(normalized)

        template_names.append(tem.split('.txt')[0].split('/')[8])

    return tot_mag_list, template_names


# simulate observed WISE fluxes by integrating templates over WISE bandpasses
def simulate_wise_fluxes(z, tem, bands, csv):
    # redshift template
    red_spec = redshift_spectrum(z, tem, False)

    shifted_wavelengths = np.array(red_spec[4])
    lumi = np.array(red_spec[1])
    # print("lumi",lumi)
    normalized = []

    # iterate through WISE bands
    for y in range(len(bands)):
        if csv:
            band = pd.read_csv(bands[y], header=None, engine='python')
        else:
            band = pd.read_csv(bands[y], header=None, delim_whitespace=True, engine='python')
        # bandpass wavelength list
        bandwaves = np.array(band.iloc[:, 0])
        # print("bandwaves", np.max(bandwaves))
        # ("min shifted_wavelengths", np.min(shifted_wavelengths))
        # bandpass response function values
        band_response = np.array(band.iloc[:, 1])
        # convolve wavelength with response function per Greg Rudnick's suggestion to account for the fact
        # that the WISE detectors are photon counting devices, while the templates are energy templates
        band_convolved = np.multiply(bandwaves, band_response)

        # trim template to same wavelength range as WISE band
        cut = np.where((shifted_wavelengths >= np.min(bandwaves)) & (shifted_wavelengths <= np.max(bandwaves)))[0]
        # ('cut',cut)
        # print(np.min(bandwaves))

        trimmed_y = shifted_wavelengths[cut]
        trimmed_L = lumi[cut]
        # print('trimmed_y',trimmed_y)
        # print('trimmed_L',trimmed_L)

        # interpolate template to wavelengths in the bandpass list, multiply by the convolved response function
        # at that wavelength
        inter_lum = []
        # print('bandwaves',bandwaves)
        # print('len(bandwaves)',len(bandwaves))
        # print('band_convolved',band_convolved)
        # print('band_convolved', len(band_convolved))

        for i in range(len(bandwaves)):
            inter_lum.append(band_convolved[i] * (np.interp(bandwaves[i], trimmed_y, trimmed_L)))

        # integrate template multiplied by response function
        spectrum = [bandwaves, inter_lum]
        interped_again = interpolate_spec(spectrum, True)
        wise_lums = integrate_spectrum(interped_again[0], interped_again[1])

        # integrate wise response function to divide out
        band_spectrum = [bandwaves, band_convolved]
        interped_band = interpolate_spec(band_spectrum, True)
        integrated_band = integrate_spectrum(interped_band[0], interped_band[1])
        # divide two
        normalized.append(wise_lums / integrated_band)
        # print(normalized)
    return normalized


# integrate the IR templates and write out the total IR luminosity so it can be recalled quickly without doing an
# integration each time
def writetotals():
    totlist = []
    # for each template
    for x in range(len(templates)):
        # call redshift function, but don't actually redshift, just trim to 8-1000 microns
        shifted_spectrum = redshift_spectrum(0, templates[x], True)
        # interpolate the template
        interped_spectrum = interpolate_spec(shifted_spectrum, False)
        # integrate template from 8-1000 micron
        total_ir = integrate_spectrum(interped_spectrum[0], interped_spectrum[1])
        totlist.append(total_ir)
    # write out the integral totals in a file
    with open(projpath + 'integrations/kirk.txt', 'wb') as fb:
        pickle.dump(totlist, fb)


writetotals()


# calculate SFRs using calibration given in Murphy+11
def murphyIRSFR(L_IR):
    L_IR = L_IR.to('erg/s').value
    SFR = 3.88e-44 * L_IR
    return SFR


# calculate IR SFRs
def IR_SFRs(z, name, pdfEnder, calc_SFR = False, tems = templates):
    # luminosity distance
    d = cosmo.luminosity_distance(z)
    # convert WISE mag to flux in Janskys
    fluxes = WISE.mag_to_flux(name)
    w3_flux = fluxes[0] * u.Jy
    w3_flux_err = fluxes[2] * u.Jy
    # print('W3:', w3_flux, w3_flux_err)
    # assume no W4 data for now
    w_four_good = False

    # calculate luminosities with fluxes & distances
    w3_lum = (w3_flux * 4 * np.pi * d ** 2).to('W/Hz')
    w3_lum_err = ((4 * np.pi * d ** 2) * w3_flux_err).to('W/Hz')
    # print('W3:', w3_lum, w3_lum_err)

    # if there's data for W4
    if not np.isnan(fluxes[1]):
        w4_flux = fluxes[1] * u.Jy
        w4_flux_err = fluxes[3] * u.Jy
        w_four_good = True
        w4_lum = (w4_flux * 4 * np.pi * d ** 2).to('W/Hz')
        w4_lum_err = ((4 * np.pi * d ** 2) * w4_flux_err).to('W/Hz')
    # print('W4:', w4_lum, w4_lum_err)

    # lists for SFR results
    SFRs = []

    # read in template total IR luminosities previously calculated
    with open(projpath + 'integrations/kirk.txt', 'rb') as fb:
        total_ir = np.array(pickle.load(fb))

    #print('tems', tems)
    # for each template
    with PdfPages(name + pdfEnder + ".pdf") as pdf:
        # fig = plt.figure()
        # tems = np.sort(tems)
        SED_list = []
        x_vals = []
        y_vals = []
        for i, tem in enumerate(tems):
            start = tem.find('y')
            tname = tem[start + 2:len(tem)]

            # redshift wavelengths of template
            tem_lum = redshift_spectrum(z, tem, False)

            # if there is W4 data, do least squares fit of W3 & W4 points to the template curve
            if w_four_good:
                # join W3 & W4 observed luminosities
                measured_lums = np.array([float(w3_lum.value), float(w4_lum.value)])
                measured_lum_errs = np.array([float(w3_lum_err.value), float(w4_lum_err.value)])

                # simulate a WISE flux by integrating the template over the response curves
                simulated = np.array(simulate_wise_fluxes(z, tem, wise_bandpasses_3_4, False))

                # perform least squares fit of observed W3, W4 luminosities to the simulated W3, W4 luminosities
                # this gives a normalization parameter which can be multiplied by the template TIR luminosity to give an
                # estimate of the intrinsic luminosity of the source
                # note: this is equation 3 from this paper: https://aip.scitation.org/doi/pdf/10.1063/1.168428
                # can use flux instead of lum , need to figure out how simulated piece works
                l_ratio = (measured_lums[0] * simulated[0] / (measured_lum_errs[0]) ** 2 + measured_lums[1] * simulated[
                    1] / (measured_lum_errs[1]) ** 2) / ((simulated[0] / measured_lum_errs[0]) ** 2 + (
                        simulated[1] / measured_lum_errs[1]) ** 2)

            # if there is no W4 data, simply take ratio of template and observed luminosity at W3
            else:
                measured_lums = np.array([float(w3_lum.value), 0.])
                measured_lum_errs = np.array([float(w3_lum_err.value), 0.])

                simulated = np.array(simulate_wise_fluxes(z, tem, wise_bandpasses_3_4, False))

                # l_ratio = float(w3_lum.value/tem_lum[2])
                l_ratio = float(w3_lum.value / simulated[0])

            # the observed LIR is just the template TIR luminosity multiplied by the normalization factor determined
            # print("len total_ir",len(total_ir))

            if calc_SFR == True:
                L_ir_tot = total_ir[i] * l_ratio * u.W
                SFR = murphyIRSFR(L_ir_tot)
                print(tem, SFR)
                SFRs.append(SFR)

                # make a plot
                fig = plt.figure()

                # ax = plt.subplot(4,3,i+1)
                # print('tem_lum[0]', tem_lum[0])
                plt.plot(tem_lum[0], tem_lum[1] * l_ratio, linewidth=0.5, label="Template")
                wave = np.array([12, 22]) / (1 + z)
                plt.scatter(wave, simulated * l_ratio, marker='s', facecolors='none', edgecolors='green',
                            label="Model Photometry")
                # plt.scatter(wave, measured_lums, color='red', marker='*')
                plt.errorbar(wave, measured_lums, yerr=measured_lum_errs, color='red', marker='*', ls='none',
                             label="Photometry")
                plt.xlim(1, 1000)
                print(SFR)
                plt.title(name + ': ' + tname + ' - SFR: ' + str(round(SFR)))
                plt.xlabel('Wavelength [microns]')
                plt.ylabel('Luminosity [W/Hz]')
                plt.legend()
                plt.xscale('log')
                plt.yscale('log')

                pdf.savefig()
                plt.close()
            else:
                tem_break = tem.split("/")
                #print(len(tem_break)-1)
                tem_trunc = tem_break[len(tem_break)-1].split('_spec')[0]
                # make a plot
                fig = plt.figure()

                # ax = plt.subplot(4,3,i+1)
                # print('tem_lum[0]', tem_lum[0])
                plt.plot(tem_lum[0], tem_lum[1] * l_ratio, linewidth=0.5, label="Template")

                #x_vals.append(tem_lum[0])
                #y_vals.append(tem_lum[1] * l_ratio)

                SED_list.append([(tem_lum[0]).tolist(), (tem_lum[1] * l_ratio).tolist()])
                wave = np.array([12, 22]) / (1 + z)
                plt.scatter(wave, simulated * l_ratio, marker='s', facecolors='none', edgecolors='green',
                            label="Model Photometry")
                # plt.scatter(wave, measured_lums, color='red', marker='*')
                plt.errorbar(wave, measured_lums, yerr=measured_lum_errs, color='red', marker='*', ls='none',
                             label="Photometry")
                plt.xlim(1, 100)
                plt.ylim(10e21, 10e26)
                plt.title(name + ': ' + tem_trunc)
                plt.xlabel('Wavelength [microns]')
                plt.ylabel('Luminosity [W/Hz]')
                plt.legend()
                plt.xscale('log')
                plt.yscale('log')

                pdf.savefig()
                plt.close()

                '''
                new_name = name + "_plotInfo.dat"
                print(new_name)
                y_vals = np.array(y_vals)
                x_vals = np.array(x_vals)
                ascii.write([micro_wave, fnu], new_name,names=['Column 1: Rest Wavelength (microns)', 'Column 2: Flux (ergs/s/cm^2/Hz)'])
                ascii.write([y_vals,x_vals], new_name,names=['Wavelength [microns]','Luminosity [W/Hz]'])
                '''
        return np.average(SFRs), np.std(SFRs),SED_list
'''
def plot_all(name,SED_list):
    with PdfPages(name + ".pdf") as pdf:
        plt.figure()
        print(len(SED_list))
        for i
       # for s in SED_list():

         #plt.plot(s[0], s[1])

        #print(SED_list)

        #plt.xlabel('Wavelength [microns]')
        #plt.ylabel('Luminosity [W/Hz]')
        #plt.legend()
        #plt.xscale('log')
        #plt.yscale('log')

        pdf.savefig()
        plt.close()
'''



