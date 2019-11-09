from time import perf_counter as clock
from os import path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import importlib

import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib
##
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_lines
import scipy.integrate as integrate
import scipy.constants as const
from matplotlib.patches import Rectangle
from os import listdir
from matplotlib.backends.backend_pdf import PdfPages
from uncertainties import ufloat
from uncertainties import unumpy as unp
ppcbfits = importlib.import_module('cb_ppxf_tools')


#### Variables that I could generalize
# np.where(yl_fit == np.amax(yl_fit))
# np.median(fluxdata[zmg_mask])


def closest_position(value,array):
    diff_arr = array - value
    smallest = np.amin(np.abs(diff_arr))
    small_pos = np.where(np.abs(diff_arr) == smallest)

    return small_pos


def calc_vel(lam_or,lam_obs): # lambda in units of Angstroms, return vel in km/s
    del_lam = lam_obs-lam_or
    velocity = del_lam*(1/lam_or)*const.c*(1/1000)

    return velocity

#gets the title for the figure based on the file that is being plotted
def get_plot_title(filepath):
    fits_file = fits.open(filepath)

    file_char_list = list(filepath)
    spec_pos = filepath.find('spec-')
    title_list = file_char_list[spec_pos:-5]
    title = ''.join(title_list)
    headers = fits_file['PRIMARY'].header
    hh = str(int(headers['PLUG_RA'] * (24/360)))
    if not (len(hh) == 2):
        hh = '0' + hh
    mm = str(int((headers['PLUG_RA'] * (24/360) * 60) % 60))
    if not (len(mm) == 2):
        mm = '0' + mm
    ss = str(np.round((headers['PLUG_RA'] * (24/360) * 60 * 60) % 60, decimals=2))
    if not (len(str(int(np.round((headers['PLUG_RA'] * (24/360) * 60 * 60) % 60)))) == 2):
        ss = '0' + ss
    title = title + '/J' + hh + mm + ss
    return title




def main(file):
    hdu = fits.open(file)
    t = hdu[1].data
    # z = float(hdu[1].header["Z"]) # SDSS redshift estimate
    # z = hdu['SPECOBJ'].data['Z'][0]
    z= ppcbfits.get_z(filepath)

    mgII = 2799.12
    mgIIl = 2795.5
    mgIIr = 2802.7

    # Only use the wavelength range in common between galaxy and stellar library.
    coadd = hdu['COADD'].data
    fluxdata = coadd['flux']
    wavelength = 10 ** coadd['loglam']
    mask = (wavelength > mgII*(1 + z)-40) & (wavelength < mgII*(1 + z)+40)
    flux = fluxdata[mask]
    galaxy = flux / np.median(flux)  # Normalize spectrum to avoid numerical issues
    wave = wavelength[mask]


    specobj = hdu[2].data #replacing 'SPECOBJ' with index 2 since  spec-6712-56421-0114 doesn't have name SPECOBJ
    spzline = hdu['SPZLINE'].data
    linename = spzline['LINENAME']
    linewave = spzline['LINEWAVE']

    zmg_mask = (wavelength>(mgII*(1 + z)-80)) & (wavelength<(mgII*(1 + z)+80))
    zmg_flux = (fluxdata[zmg_mask] - np.median(fluxdata[zmg_mask]))*(-1)
    zmg_wave = wavelength[zmg_mask]

    mgrndx = np.where(zmg_flux == np.amax(zmg_flux))[0][0]

    spectrum = Spectrum1D(flux=zmg_flux*u.Jy, spectral_axis=zmg_wave*u.AA)

    mgIIl_init = models.Gaussian1D(amplitude=8*u.Jy, mean=zmg_wave[mgrndx]*u.AA, stddev=10*u.AA)
    mgIIr_init = models.Gaussian1D(amplitude=6*u.Jy, mean=zmg_wave[mgrndx-10]*u.AA, stddev=10*u.AA)
    mgl_fit, mgr_fit = fit_lines(spectrum, [mgIIl_init,mgIIr_init], window=10*u.AA)

    std = mgl_fit.stddev/u.AA
    amp = mgl_fit.amplitude/u.Jy
    mean = mgl_fit.mean/u.AA

    yl_fit = mgl_fit(zmg_wave*u.AA)
    yr_fit = mgr_fit(zmg_wave*u.AA)
    gauss_plot = (yr_fit/u.Jy-np.median(fluxdata[zmg_mask]))*(-1)

    square_l = closest_position(zmg_wave[np.where(yl_fit == np.amax(yl_fit))]-(3*std),zmg_wave)[0][0]
    square_r = closest_position(zmg_wave[np.where(yl_fit == np.amax(yl_fit))]+(3*std),zmg_wave)[0][0]


    def get_EW_from_flux(): # recent approach suggested by Alex in blackboard; using 3 sigma to locate absorption
        dlambda = np.array([zmg_wave[x+1]-zmg_wave[x] for x in range(square_l,square_r-1)])
        EqW = np.sum(dlambda *
                     u.AA *
                     (np.median(fluxdata[zmg_mask]) - fluxdata[zmg_mask][square_l:square_r-1]) *
                     (1/np.median(fluxdata[zmg_mask])))

        return EqW


    def get_area_99(): #old method to get equivalent width using area under gaussian
        dlambda = np.array([zmg_wave[x + 1] - zmg_wave[x] for x in range(square_l,square_r-1)])
        gaussian_area = np.sum(dlambda * u.AA * yl_fit[square_l:square_r-1])
        return gaussian_area


    gaussian_area = get_area_99()

    EqW = get_EW_from_flux() # newer equivalent width

    ## Integral method?
    alt_area = integrate.quad(lambda x: amp*(np.e**((-1/2)*(((x-mean)/std)**2))),-np.inf,np.inf)

    # EW = gaussian_area/np.median(fluxdata[zmg_mask]*u.Jy) # old equivalent width calculated with gaussian
    plt.fill_between(zmg_wave*u.AA, 0*u.Jy, np.median(fluxdata[zmg_mask])*u.Jy)


    # get centroid velocity

    centroid_vel = calc_vel(mgII*(1+z), zmg_wave[np.where(yl_fit == np.amax(yl_fit))])
    #print("---Centroid velocity: ", centroid_vel)

    ## Code to find V50 and V90.

    # Determines the area of absorption such that no values above the continuum are included
    # Starting from the calculated center of emission, it returns the spread of the emission line in array index distance.


    def get_abs_spread(data, type):
        cent_indx, cont_median, wave, flux = data[0], data[1], data[2], data[3]
        searching = True
        if type == "blue":
            spread = -1
            while searching:
                if flux[cent_indx + spread] < cont_median:
                    spread -= 1
                if flux[cent_indx + spread] > cont_median:
                    spread += 1
                    searching = False
        if type == "red":
            spread = 1
            while searching:
                if flux[cent_indx + spread] < cont_median:
                    spread += 1
                if flux[cent_indx + spread] > cont_median:
                    spread -= 1
                    searching = False
        return spread


    abs_rspread = get_abs_spread([np.where(yl_fit == np.amax(yl_fit))[0][0],
                                 np.median(fluxdata[zmg_mask]),
                                 zmg_wave,
                                 fluxdata[zmg_mask]],
                                 "red")
    abs_bspread = get_abs_spread([np.where(yl_fit == np.amax(yl_fit))[0][0],
                                 np.median(fluxdata[zmg_mask]),
                                 zmg_wave,
                                 fluxdata[zmg_mask]],
                                 "blue")

    tot_area = EqW * np.median(fluxdata[zmg_mask]) *(1/u.AA) # dimensionless but we know it's angstroms

    #s
    def calc_percentage_pos(data):
        tot_area, percent, cont_median, red_end, wave, flux = data[0], data[1], data[2], data[3], data[4], data[5]

        curr_area = 0
        indx_span = 0
        while curr_area < tot_area*(percent/100):
            indx_span += 1
            dlambda = np.array([wave[x + 1] - wave[x] for x in range(red_end-indx_span, red_end)])
            curr_area = np.sum(dlambda *
                         (cont_median - flux[red_end-indx_span:red_end]))

        return indx_span, curr_area

    v50_indx_span, v50_area = calc_percentage_pos([tot_area,
                                              50,
                                              np.median(fluxdata[zmg_mask]),
                                              np.where(yl_fit == np.amax(yl_fit))[0][0]+abs_rspread,
                                              zmg_wave,
                                              fluxdata[zmg_mask]])

    v98_indx_span, v98_area = calc_percentage_pos([tot_area,
                                              98,
                                              np.median(fluxdata[zmg_mask]),
                                              np.where(yl_fit == np.amax(yl_fit))[0][0]+abs_rspread,
                                              zmg_wave,
                                              fluxdata[zmg_mask]])

    # Calculate velocity
    red_end = np.where(yl_fit == np.amax(yl_fit))[0][0]+abs_rspread

    v50_wav = zmg_wave[red_end - v50_indx_span]
    v98_wav = zmg_wave[red_end - v98_indx_span]

    v50_vel = (v50_wav - (mgII * (1 + z)))*(1/(mgII * (1 + z)))*const.c*(1/1000)
    v98_vel = (v98_wav - (mgIIl * (1 + z)))*(1/(mgIIl * (1 + z)))*const.c*(1/1000)


    # plotting

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    matplotlib.rcParams.update({'font.size': 6})
    plot_title = get_plot_title(file)
    fig.suptitle(plot_title)


    rect = Rectangle((zmg_wave[np.where(yl_fit == np.amax(yl_fit))]-(EqW/2)/u.AA, 0), EqW/u.AA,
                     np.median(fluxdata[zmg_mask]), alpha=0.4, color='red', label='EW Rectangle')

    ax.axvspan(zmg_wave[np.where(yl_fit == np.amax(yl_fit))]-(3*std),
               zmg_wave[np.where(yl_fit == np.amax(yl_fit))]+(3*std), alpha=0.3, color='C5', label='Area of emission')

    ax.add_patch(rect)
    ax.grid(True)

    ax.axvline(x=zmg_wave[np.where(yl_fit == np.amax(yl_fit))], color='blue', label='Centroid', alpha=0.4)
    ax.axhline(y=np.median(fluxdata[zmg_mask]), color='cyan', label='Continuum', alpha=0.4)

    #ax.plot(zmg_wave, zmg_flux, color='black', linewidth=0.3)
    # ax.plot(zmg_wave, yl_fit, color='green', linewidth=0.3)
    #ax.plot(zmg_wave, yr_fit, color='orange', linewidth=0.3)
    ax.plot(zmg_wave, gauss_plot, color='m', linewidth=0.3)
    ax.plot(zmg_wave, fluxdata[zmg_mask], color='black', linewidth=0.3)
    ax.set_title("General Information Plot; Centroid Vel: "+ str(centroid_vel) + r"$\frac{km}{s}$")
    ax.tick_params(labelsize=6)
    ax.set_xlabel("$\AA ngstr \ddot{o} ms$")
    ax.set_ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")


    ax.legend()

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(zmg_wave, fluxdata[zmg_mask], color='black', linewidth=0.3)
    ax2.axvline(x=(mgII * (1 + z)), color='red', label='Expected Centroid', alpha=0.4)
    ax2.axvline(x=zmg_wave[np.where(yl_fit == np.amax(yl_fit))], color='green', label='Actual Centroid', alpha=0.4)
    ax2.axhline(y=np.median(fluxdata[zmg_mask]), color='cyan', label='Continuum', alpha=0.4)
    ax2.axvspan(zmg_wave[red_end - v50_indx_span],
               zmg_wave[red_end], alpha=0.3, color='m', label='Area of v50')
    ax2.set_title("v50 Info Plot; $v=$" + str(v50_vel) + r"$\frac{km}{s}$")
    ax2.set_xlabel("$\AA ngstr \ddot{o} ms$")
    ax2.set_ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")


    ax2.legend()
    ax2.grid(True)


    ax3 = fig.add_subplot(2, 2, 4)
    ax3.plot(zmg_wave, fluxdata[zmg_mask], color='black', linewidth=0.3)
    ax3.axvline(x=(mgIIl * (1 + z)), color='blue', label='Expected MgII Blue', alpha=0.4)
    ax3.axvline(x=(mgII * (1 + z)), color='red', label='Expected Centroid', alpha=0.4)
    ax3.axvline(x=zmg_wave[np.where(yl_fit == np.amax(yl_fit))], color='green', label='Emission Centroid', alpha=0.4)
    ax3.axhline(y=np.median(fluxdata[zmg_mask]), color='cyan', label='Continuum', alpha=0.4)
    ax3.axvspan(zmg_wave[red_end - v98_indx_span],
               zmg_wave[red_end], alpha=0.3, color='m', label='Area of v98')
    ax3.set_title("v98 Info Plot; $v=$" + str(v98_vel) + r"$\frac{km}{s}$")
    ax3.set_xlabel("$\AA ngstr \ddot{o} ms$")
    ax3.set_ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")


    ax3.legend()
    ax3.grid(True)



galfolpath = 'C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hires/galaxies_fits_DR15/'
fitsfiles = [f for f in listdir(galfolpath)]
# spec-1305-52757-0191.fits is interesting. It might atcually have MgII
# blacklist = ['spec-0326-52375-0382.fits',
#              'spec-0566-52238-0319.fits',
#              'spec-0639-52146-0388.fits',
#              'spec-0986-52443-0580.fits',
#              'spec-1305-52757-0191.fits']

blacklist = ['spec-0326-52375-0382.fits',
             'spec-0566-52238-0319.fits',
             'spec-0639-52146-0388.fits',
             'spec-0986-52443-0580.fits']

with PdfPages('V50_V98.pdf') as pdf:
    for filename in fitsfiles:
        if filename not in blacklist:
            filepath = galfolpath + filename
            print('\nPlotting file: ', filename, '..........')
            main(filepath)
            pdf.savefig()
            plt.close('all')
            print('Success')


print('-----------Finished-----------!')

#plt.show(block=False)
