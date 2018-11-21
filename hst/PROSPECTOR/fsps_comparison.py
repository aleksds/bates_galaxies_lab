import os
import numpy as np
import fsps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                            sfh=0, logzsol=0.0, dust_type=2, dust2=0)
print('done with fsps.StellarPopulation')

wave30, spec30 = sp.get_spectrum(tage=0.03)
wave8, spec8 = sp.get_spectrum(tage=0.008)
wave6, spec6 = sp.get_spectrum(tage=0.006)
print('done with sp.get_spectrum')

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
lsun = 3.8e33

bands = np.array(['galex_fuv','galex_nuv','wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])

phottable = 'Photometry/coarse_final/atest.txt'
catalog = ascii.read(phottable)
ad_values = ascii.read('../autogalfit/ad_mag_size_table.dat')

mfuv = np.array(ad_values['mfuv'])
mnuv = np.array(ad_values['mnuv'])
ebv = np.array(ad_values['ebv'])

filename = 'fsps_comparison.pdf'

with PdfPages(filename) as pdf:
    
    for j in range(0,len(catalog)):

        distance = cosmo.luminosity_distance(catalog['z'][j]).to('cm') / u.cm

        mags30 = sp.get_mags(tage=0.03, bands=bands, redshift=catalog['z'][j])
        mags8 = sp.get_mags(tage=0.008, bands=bands, redshift=catalog['z'][j])
        mags6 = sp.get_mags(tage=0.006, bands=bands, redshift=catalog['z'][j])
        
        m814 = -2.5*np.log10(catalog['f_814'][j]*1e-9)
        m475 = -2.5*np.log10(catalog['f_475'][j]*1e-9)
        m160 = -2.5*np.log10(catalog['f_1600'][j]*1e-9)
        m814_cor = m814 - ebv[j]*1.536
        m475_cor = m475 - ebv[j]*3.248
        m160_cor = m160 - ebv[j]*0.512
        mnuv_cor = mnuv[j] - ebv[j]*7.06
        mfuv_cor = mfuv[j] - ebv[j]*4.37
        mass30 = 10**((m814_cor - mags30[3])/(-2.5))
        mass8 = 10**((m814_cor - mags8[3])/(-2.5))
        mass6 = 10**((m814_cor - mags6[3])/(-2.5))

        flux_model30 = spec30 * lsun * mass30 / (4. * np.pi * distance**2) * (1+catalog['z'][j])
        flux_model8 = spec8 * lsun * mass8 / (4. * np.pi * distance**2) * (1+catalog['z'][j])
        flux_model6 = spec6 * lsun * mass6 / (4. * np.pi * distance**2) * (1+catalog['z'][j])

        phot_wave = np.array([1516.,2267.,4750.,8140.,16000.])
        mags30_flux = 10**(mags30/(-2.5))*mass30*3631.*1e-23
        mags8_flux = 10**(mags8/(-2.5))*mass8*3631.*1e-23
        mags6_flux = 10**(mags6/(-2.5))*mass6*3631.*1e-23
        
        phot_wave_rest = phot_wave / (1+catalog['z'][j])
        phot_mag = np.array([mfuv_cor, mnuv_cor, m475_cor, m814_cor, m160_cor])
        phot_flux =10**(phot_mag / (-2.5))*3631*1e-23
        
        fig = plt.figure()
        plt.suptitle(catalog['ID'][j])
        
        ax = fig.add_subplot(2,2,1)
        plt.xlim(0.2500,1.8000)
        #plt.ylim(2e-28,5e-27)
        #plt.xscale("log", nonposx='clip')
        #plt.yscale("log", nonposy='clip')
        plt.ylim(0.,3.)
        plt.xlabel(r'Observed Wavelength [$\mu$m]')
        plt.ylabel(r'Flux [10$^{-27}$ erg/s/cm$^2$/Hz]')
        ax.xaxis.set_label_position('top')
        
        plt.plot(wave30*(1+catalog['z'][j])/1e4,flux_model30/1e-27,color='red', linewidth=0.1)
        plt.plot(wave8*(1+catalog['z'][j])/1e4,flux_model8/1e-27,color='green', linewidth=0.1)
        plt.plot(wave6*(1+catalog['z'][j])/1e4,flux_model6/1e-27,color='blue', linewidth=0.1)
        
        plt.scatter(phot_wave/1e4, mags30_flux/1e-27, marker='o', color='red')
        plt.scatter(phot_wave/1e4, mags8_flux/1e-27, marker='v', color='green')
        plt.scatter(phot_wave/1e4, mags6_flux/1e-27, marker='s', color='blue')

        plt.scatter(phot_wave/1e4, phot_flux/1e-27, color='black', marker='*')

        ax = fig.add_subplot(2,2,3)
        ax.set_xlim([1000,5500])
        #ax.set_ylim([1e-29,5e-27])
        #ax.set_xscale("log", nonposx='clip')
        #ax.set_yscale("log", nonposy='clip')
        plt.ylim(0.,3)
        ax.set_xlabel(r'Observed Wavelength [$\AA$]')
        ax.set_ylabel(r'Flux [10$^{-27}$ erg/s/cm$^2$/Hz]')

        plt.plot(wave30*(1+catalog['z'][j]),flux_model30/1e-27,color='red', linewidth=0.1)
        plt.plot(wave8*(1+catalog['z'][j]),flux_model8/1e-27,color='green', linewidth=0.1)
        plt.plot(wave6*(1+catalog['z'][j]),flux_model6/1e-27,color='blue', linewidth=0.1)
        
        plt.scatter(phot_wave, mags30_flux/1e-27, marker='o', color='red')
        plt.scatter(phot_wave, mags8_flux/1e-27, marker='v', color='green')
        plt.scatter(phot_wave, mags6_flux/1e-27, marker='s', color='blue')

        plt.scatter(phot_wave, phot_flux/1e-27, color='black', marker='*')
        
        # make optical color-color plot
        ax = fig.add_subplot(2,2,2)
        ax.set_xlabel('[F814W] - [F160W]')
        ax.set_ylabel('[F475W] - [F814W]')
        ax.set_xlim([-1.0,1.0])
        ax.set_ylim([-1.0,1.0])
        ax.xaxis.set_label_position('top')
        ax.yaxis.set_label_position('right')
        
        # plot points for the three models
        plt.scatter(phot_mag[3]-phot_mag[4], phot_mag[2]-phot_mag[3], color='black', marker='*', label='data')
        plt.errorbar(phot_mag[3] - phot_mag[4], phot_mag[2] - phot_mag[3], xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='black',linewidth=0.5)
        plt.scatter(mags30[3]-mags30[4], mags30[2]-mags30[3], color='red', marker='o', label='30 Myr')
        plt.scatter(mags8[3]-mags8[4], mags8[2]-mags8[3], color='green', marker='v', label='8 Myr')
        plt.scatter(mags6[3]-mags6[4], mags6[2]-mags6[3], color='blue', marker='s', label='6 Myr')

        # GALEX color-color panel
        ax = fig.add_subplot(2,2,4)
        ax.set_xlabel('[NUV] - [F475W]')
        ax.set_ylabel('[FUV] - [NUV]')
        ax.set_xlim([-0.5,3.0])
        ax.set_ylim([-0.5,3.0])
        ax.yaxis.set_label_position('right')

        
        # plot data point with error bar        
        #ax.errorbar(phot_mag[1] - phot_mag[2], phot_mag[0] - phot_mag[1], xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='magenta')
        plt.scatter(phot_mag[1]-phot_mag[2], phot_mag[0]-phot_mag[1], color='black', marker='*', label='data')
        plt.errorbar(phot_mag[1]-phot_mag[2], phot_mag[0]-phot_mag[1], xerr = 0.1, yerr=0.3, color='black', linewidth=0.5)

        # plot data point after galactic extinction correction
        #mnuv_cor = mnuv[j] - ebv[j]*7.06
        #mfuv_cor = mfuv[j] - ebv[j]*4.37
        #ax.errorbar(mnuv_cor - m475_cor, mfuv_cor- mnuv_cor, xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='purple')

        # plot points for the three models
        ax.scatter(mags30[1]-mags30[2], mags30[0]-mags30[1], color='red', marker='o', label='30 Myr')
        ax.scatter(mags8[1]-mags8[2], mags8[0]-mags8[1], color='green', marker='v', label='8 Myr')
        ax.scatter(mags6[1]-mags6[2], mags6[0]-mags6[1], color='blue', marker='s', label='6 Myr')

        pdf.savefig()
        plt.close()

os.system('open %s &' % filename)
