import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import os
from astropy import units as u
import math
import pickle
import glob
from astropy.table import Table
import random
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as st
from scipy.stats import rv_discrete

#names = ['J0106', 'J0826', 'J0827', 'J0905', 'J0908', 'J0944', 'J1039', 'J1107', 'J1125', 'J1219', 'J1229', 'J1232',
#         'J1248', 'J1341', 'J1506', 'J1613', 'J2116', 'J2118', 'J2140', 'J2256']
#names = ['J0901', 'J1506', 'J1622', 'J1713']
#names = ['J1506B']
#names = ['J0106', 'J0826']
#redshift = [0.454, 0.603]
#names = ['J0106']
names = ['J1558']
#redshift = [0.454]
redshift = [0.402]
# redshift = [0.454, 0.603, 0.681, 0.711, 0.502, 0.514, 0.634, 0.467, 0.519, 0.451, 0.614, 0.401, 0.632, 0.661, 0.437,
#              0.449, 0.728, 0.459, 0.751, 0.727]
SFR_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
SFR_std_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# convert from WISE magnitudes to fluxes in Jy
flux3s = []


# convert from WISE magnitudes to fluxes in Jy
def mag_to_flux(name):
    WISEdir = projpath + 'unWISE/%s.fits' % name
    t = Table.read(WISEdir)
    print("name: ", name)
    # constants given at http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#example
    w_three_const = 31.674
    w_four_const = 8.363
    # calculate W3, W4 fluxes in Jy
    flux_three = w_three_const * (10 ** (-(float(t['w3_mag']) / 2.5)))
    flux_four = w_four_const * (10 ** (-(float(t['w4_mag']) / 2.5)))
    # errors computed by error propagation formula
    flux_three_err = w_three_const * np.log(10) / 2.5 * (10 ** (-(float(t['w3_mag']) / 2.5))) * float(t['w3_mag_err'])
    flux_four_err = w_four_const * np.log(10) / 2.5 * (10 ** (-(float(t['w4_mag']) / 2.5))) * float(t['w4_mag_err'])
    #print("flux three:", flux_three)
    #flux3s.append(flux_three)

    mu3, sigma3 = flux_three, flux_three_err
    s3 = np.random.normal(mu3, sigma3, 200) #remember to switch back to 100
    mu4, sigma4 = flux_four, flux_four_err
    s4 = np.random.normal(mu4, sigma4, 200) #remember to switch back to 100
    #np.random.choice(s3) #did this get used?
    #print("flux_three:",np.array([flux_three]))
    #print("s3:", s3)
    s3_awesome = np.concatenate((np.array([flux_three]),s3))
    s4_awesome = np.concatenate((np.array([flux_four]),s4))
    return s3_awesome, s4_awesome, flux_three_err, flux_four_err


# calculate WISE colors and errors
def colors(name):
    pat = projpath + 'unWISE/%s.fits' % name
    t = Table.read(pat)
    # W1-W2
    one_two = float(t['w1_mag']) - float(t['w2_mag'])
    one_two_err = np.sqrt((float(t['w1_mag_err'])) ** 2 + (float(t['w2_mag_err'])) ** 2)
    # W3-W4
    three_four = float(t['w3_mag']) - float(t['w4_mag'])
    three_four_err = np.sqrt((float(t['w3_mag_err'])) ** 2 + (float(t['w4_mag_err'])) ** 2)
    # W2-W3
    two_three = float(t['w2_mag']) - float(t['w3_mag'])
    two_three_err = np.sqrt((float(t['w2_mag_err'])) ** 2 + (float(t['w3_mag_err'])) ** 2)
    return one_two, three_four, two_three, one_two_err, three_four_err, two_three_err


# TEMPLATES50: _______________________________________________________________
# speed of light
c = 299792458.  # m/s
# path to project
projpath = os.getcwd() + '/ir/'  # '../'#'/Users/graysonpetter/Desktop/IRSFRs/'
# set cosmology for calculating luminosity distance
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# Read in all templates except for AGN (SFGs and Composite)
templates = glob.glob(projpath + 'Comprehensive_library/SFG*.txt')
templates.extend(
    glob.glob(projpath + 'Comprehensive_library/Comp*.txt'))
# print('templates:', templates)
# read in the W3 & W4 bandpasses
wise_bandpasses_3_4 = sorted(glob.glob(projpath + 'bandpass/*.txt'))[2:4]


# redshift a template SED (wavelengths only)
def redshift_spectrum(z, template, trim):
    t = pd.read_csv(template, delim_whitespace=True, engine='python', header=None, skiprows=3)

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

    normalized = []

    # iterate through WISE bands
    for y in range(len(bands)):
        if csv:
            band = pd.read_csv(bands[y], header=None, engine='python')
        else:
            band = pd.read_csv(bands[y], header=None, delim_whitespace=True, engine='python')
        # bandpass wavelength list
        bandwaves = np.array(band.iloc[:, 0])
        # bandpass response function values
        band_response = np.array(band.iloc[:, 1])
        # convolve wavelength with response function per Greg Rudnick's suggestion to account for the fact
        # that the WISE detectors are photon counting devices, while the templates are energy templates
        band_convolved = np.multiply(bandwaves, band_response)

        # trim template to same wavelength range as WISE band
        cut = np.where((shifted_wavelengths >= np.min(bandwaves)) & (shifted_wavelengths <= np.max(bandwaves)))[0]
        trimmed_y = shifted_wavelengths[cut]
        trimmed_L = lumi[cut]

        # interpolate template to wavelengths in the bandpass list, multiply by the convolved response function
        # at that wavelength
        inter_lum = []
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
def IR_SFRs(z, i, tems=templates):
    # luminosity distance
    d = cosmo.luminosity_distance(z)
    # convert WISE mag to flux in Janskys
    fluxes = mag_to_flux(i)
    #we could re-write mag_to_flux so that this gives us all the 100 flux values we want
    #Then we could start a for loop right here that runs the code below 100 times and calculate SFR each time
    #print(fluxes[0])
    SFRs=[]
    for j in range(0,len(fluxes[0])):
        w3_flux = random.choice(fluxes[0]) * u.Jy
        #w3_flux = fluxes[0]*u.Jy
        w3_flux_err = fluxes[2] * u.Jy
        print('W3:', w3_flux, w3_flux_err)
        # assume no W4 data for now
        w_four_good = False

        # calculate luminosities with fluxes & distances
        w3_lum = (w3_flux * 4 * np.pi * d ** 2).to('W/Hz')
        w3_lum_err = ((4 * np.pi * d ** 2) * w3_flux_err).to('W/Hz')

        # if there's data for W4
        if not np.isnan(fluxes[1][0]): #if not np.isnan(fluxes[1][0]):
            w4_flux = random.choice(fluxes[1]) * u.Jy
            #w4_flux = fluxes[1] * u.Jy
            w4_flux_err = fluxes[3] * u.Jy
            #print('W4:', w4_flux, w4_flux_err)
            w_four_good = True
            w4_lum = (w4_flux * 4 * np.pi * d ** 2).to('W/Hz')
            w4_lum_err = ((4 * np.pi * d ** 2) * w4_flux_err).to('W/Hz')
        #print("j,w3_flux,w4_flux",j,w3_flux,w4_flux)
        # lists for SFR results
        #SFRs = []

        # read in template total IR luminosities previously calculated
        with open(projpath + 'integrations/kirk.txt', 'rb') as fb:
            total_ir = np.array(pickle.load(fb))

        print('tems', tems)
        # for each template
        for i, tem in enumerate(tems):
            print('tem: ', tem)
            # redshift wavelengths of template
            tem_lum = redshift_spectrum(z, tem, False)

            # if there is W4 data, do least squares fit of W3 & W4 points to the template curve
            if w_four_good:
                # join W3 & W4 observed luminosities
                measured_lums = np.array([float(w3_lum.value), float(w4_lum.value)])
                measured_lum_errs = np.array([float(w3_lum_err.value), float(w4_lum_err.value)])

                # simulate a WISE flux by integrating the template over the response curves
                simulated = np.array(simulate_wise_fluxes(z, tem, wise_bandpasses_3_4, False))
                print('simulated WISE flux:', simulated)

                # perform least squares fit of observed W3, W4 luminosities to the simulated W3, W4 luminosities
                # this gives a normalization parameter which can be multiplied by the template TIR luminosity to give an
                # estimate of the intrinsic luminosity of the source
                l_ratio = (measured_lums[0] * simulated[0] / (measured_lum_errs[0]) ** 2 + measured_lums[1] * simulated[
                    1] / (measured_lum_errs[1]) ** 2) / ((simulated[0] / measured_lum_errs[0]) ** 2 + (
                        simulated[1] / measured_lum_errs[1]) ** 2)

            # if there is no W4 data, simply take ratio of template and observed luminosity at W3
            else:
                l_ratio = float(w3_lum.value / tem_lum[2])

            # the observed LIR is just the template TIR luminosity multiplied by the normalization factor determined
            L_ir_tot = total_ir[i] * l_ratio * u.W

            SFR = murphyIRSFR(L_ir_tot)
            #print(tem, SFR)
            SFRs.append(SFR)
            #print("SFRs", SFRs)
            # percent68=SFRs.interval(.68)
            # percent50=SFRs.interval(.50)
            # percent16=SFRs.interval(.16)
            # percent84=SFRs.interval(.84)

            # percent68 = st.t.interval(alpha=0.68, df=len(SFRs) - 1, loc=np.mean(SFRs), scale=st.sem(SFRs))
            # percent50 = st.t.interval(alpha=0.50, df=len(SFRs) - 1, loc=np.mean(SFRs), scale=st.sem(SFRs))
            # percent16 = st.t.interval(alpha=0.16, df=len(SFRs) - 1, loc=np.mean(SFRs), scale=st.sem(SFRs))
            # percent84 = st.t.interval(alpha=0.84, df=len(SFRs) - 1, loc=np.mean(SFRs), scale=st.sem(SFRs))

    percent68 = np.percentile(SFRs, 68)
    percent50 = np.percentile(SFRs, 50)
    percent16 = np.percentile(SFRs, 16)
    percent84 = np.percentile(SFRs, 84)
    #print("DELETE THIS LATER: ", SFRs)
    #print("REALLLLLLLLLY DELETE THIS: ", percent16, percent50, percent84)



    return np.average(SFRs), np.std(SFRs), SFRs, percent68, percent50, percent16, percent84


# READ IN NAMES, REDSHIFTS: _______________________________________________________________
SFR = []
SFR_std = []
dfSFR = pd.DataFrame()
SFR68=[]
SFR50=[]
SFR16=[]
SFR84=[]



for n, r, s, d in zip(names, redshift, SFR_list, SFR_std_list):
    def gal_sfr2(n, r):
        SFR2 = IR_SFRs(r, n)
        # percent68 = np.percentile(SFRs, 68)
        # percent50 = np.percentile(SFRs, 50)
        # percent16 = np.percentile(SFRs, 16)
        # percent84 = np.percentile(SFRs, 84)
        # print("DELETE THIS LATER: ", SFRs)
        # print("REALLLLLLLLLY DELETE THIS: ", percent16, percent50, percent84)
        return SFR2
    #call a function that returns 100 w3 flux values for one galaxy
    #start a for loop that runs through 100 times
    #at the end of each iteration save the sfr value
    #need to figure out how to pass the right fluxes to IR_SFRs
    #alternatively, we can set things up so that IR SFRs automatically runs through 100 different times
    value = gal_sfr2(n, r)
    # percent68 = np.percentile(SFRs, 68)
    # percent50 = np.percentile(SFRs, 50)
    # percent16 = np.percentile(SFRs, 16)
    # percent84 = np.percentile(SFRs, 84)
    # print("DELETE THIS LATER: ", SFRs)
    # print("REALLLLLLLLLY DELETE THIS: ", percent16, percent50, percent84)
    print("value: ", value)
    SFR.append(value[0])
    SFR_std.append(value[1])
    dfSFR[n]=value[2]
    SFR68.append(value[3])
    SFR50.append(value[4])
    SFR16.append(value[5])
    SFR84.append(value[6])


print(SFR)
print(SFR_std)
print(zip(names, redshift, SFR, SFR_std))
print(flux3s)

df = pd.DataFrame(list(zip(names, redshift, SFR, SFR_std,SFR68, SFR50, SFR16, SFR84)), columns=['name', 'redshift', 'SFR', 'SFRstd','SFR68', 'SFR50', 'SFR16', 'SFR84'])
print(df)
#print(dfSFR)

df.to_csv(r'/Users/natashajones/Desktop/PycharmProjects/HIZEA-IR/galSFRsMC_w4.csv', index=False)

dfSFR.to_csv(r'/Users/natashajones/Desktop/PycharmProjects/HIZEA-IR/Mc_Squiggles.csv', index=False)

table = fits.open('hizea_wind_ancillary.fit')
logsfr=table[1].data['SFR']
sfr=10**logsfr

short_name = table[1].data['SHORT_NAME']
sfr_names = ['J0106-1023', 'J0826+4305','J0827+2954','J0905+5759','J0908+1039','J0944+0930','J1039+4537','J1107+0417','J1125-0145','J1219+0336','J1229+3545','J1232+0723','J1248+0601','J1341-0321','J1506+6131','J1613+2834','J2116-0634','J2118+0017','J2140+1209','J2256+1504']
tremonti_sfrs = np.zeros(len(sfr_names))
for i in range (0,len(sfr_names)):
    match=np.where(short_name == sfr_names[i])[0][0]
    tremonti_sfrs[i]=sfr[match]

filename = 'SFR_distribution1.pdf'
with PdfPages(filename) as pdf:
    # for i, name, tremonti in zip(len(names), names, tremonti_sfrs):
    for i in range(0,len(names)):
        fig, ax= plt.subplots()
        ax.hist(dfSFR[names[i]], 20, color = 'blue', label='Monte Carlo Values') #remember to switch back to 10
        ax.axvline(tremonti_sfrs[i], color = 'orange',label = "Tremonti SFR: {:.1f}".format(tremonti_sfrs[i]))
        ax.axvline(SFR50[i], color='red', label="50th percentile: {:.1f}".format(SFR50[i]))
        ax.axvline(SFR16[i], color='black', label="16th percentile: {:.1f}".format(SFR16[i]))
        ax.axvline(SFR84[i], color='purple', label="84th percentile: {:.1f}".format(SFR84[i]))
        ax.legend()
        ax.set_xlabel("SFR")
        ax.set_ylabel("Number")
        column = names[i]
        ax.set_title("Distribution of SFRs from Monte Carlo Error Propagation: {}".format(column))
        pdf.savefig()
        plt.close()