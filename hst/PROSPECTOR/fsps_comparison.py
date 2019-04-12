import os
import numpy as np
import fsps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from datetime import datetime

start = datetime.now()

dust_arr = [0.2, 0.4, 0.6, 0.8, 1]

for a in range(0, len(dust_arr)):
    dust = dust_arr[a]

    #dust = 0.4 #try different values here. what values are usable?

    sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                    sfh=0, logzsol=0.0, dust_type=2, dust2=dust)
    print('done with fsps.StellarPopulation', datetime.now()-start)

    age = np.array([30,8,6,5]) #try different values here. What values are usable here?  

    wave_one, spec_one = sp.get_spectrum(tage=age[0]/1e3)
    print('done with one sp.get_spectrum', datetime.now()-start)
    wave_two, spec_two = sp.get_spectrum(tage=age[1]/1e3)
    print('done with two sp.get_spectrum', datetime.now()-start)
    wave_thr, spec_thr = sp.get_spectrum(tage=age[2]/1e3)
    print('done with thr sp.get_spectrum', datetime.now()-start)
    wave_fou, spec_fou = sp.get_spectrum(tage=age[3]/1e3)
    print('done with four sp.get_spectrum', datetime.now()-start)


    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    lsun = 3.8e33

    bands = np.array(['galex_fuv','galex_nuv','wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])
    phottable = 'Photometry/coarse_final/atest.txt'
    catalog = ascii.read(phottable)
    ad_values = ascii.read('../autogalfit/ad_mag_size_table.dat')

    mfuv = np.array(ad_values['mfuv'])
    mnuv = np.array(ad_values['mnuv'])
    ebv = np.array(ad_values['ebv'])

    filename = 'fsps_comparison_'+str(dust)+'_'+str(age[0])+'_'+str(age[1])+'_'+str(age[2])+'_'+str(age[3])+'.pdf'

    with PdfPages(filename) as pdf:
    
        for j in range(0,len(catalog)):

            distance = cosmo.luminosity_distance(catalog['z'][j]).to('cm') / u.cm

            mags_one = sp.get_mags(tage=age[0]/1e3, bands=bands, redshift=catalog['z'][j])
            mags_two = sp.get_mags(tage=age[1]/1e3, bands=bands, redshift=catalog['z'][j])
            mags_thr = sp.get_mags(tage=age[2]/1e3, bands=bands, redshift=catalog['z'][j])
            mags_fou = sp.get_mags(tage=age[3]/1e3, bands=bands, redshift=catalog['z'][j])

        
            m814 = -2.5*np.log10(catalog['f_814'][j]*1e-9)
            m475 = -2.5*np.log10(catalog['f_475'][j]*1e-9)
            m160 = -2.5*np.log10(catalog['f_1600'][j]*1e-9)
            m814_cor = m814 - ebv[j]*1.536
            m475_cor = m475 - ebv[j]*3.248
            m160_cor = m160 - ebv[j]*0.512
            #mnuv_cor = mnuv[j] - ebv[j]*7.06
            #mfuv_cor = mfuv[j] - ebv[j]*4.37
            mnuv_cor = mnuv[j] - ebv[j]*6.892
            mfuv_cor = mfuv[j] - ebv[j]*6.738
            mass_one = 10**((m814_cor - mags_one[3])/(-2.5))
            mass_two = 10**((m814_cor - mags_two[3])/(-2.5))
            mass_thr = 10**((m814_cor - mags_thr[3])/(-2.5))
            mass_fou = 10**((m814_cor - mags_fou[3])/(-2.5))
        
            flux_model_one = spec_one * lsun * mass_one / (4. * np.pi * distance**2) * (1+catalog['z'][j])
            flux_model_two = spec_two * lsun * mass_two / (4. * np.pi * distance**2) * (1+catalog['z'][j])
            flux_model_thr = spec_thr * lsun * mass_thr / (4. * np.pi * distance**2) * (1+catalog['z'][j])
            flux_model_fou = spec_fou * lsun * mass_fou / (4. * np.pi * distance**2) * (1+catalog['z'][j])
        
            phot_wave = np.array([1516.,2267.,4750.,8140.,16000.])
            mags_one_flux = 10**(mags_one/(-2.5))*mass_one*3631.*1e-23
            mags_two_flux = 10**(mags_two/(-2.5))*mass_two*3631.*1e-23
            mags_thr_flux = 10**(mags_thr/(-2.5))*mass_thr*3631.*1e-23
            mags_fou_flux = 10**(mags_fou/(-2.5))*mass_fou*3631.*1e-23
        
            phot_wave_rest = phot_wave / (1+catalog['z'][j])
            phot_mag = np.array([mfuv[j], mnuv[j], m475, m814, m160])
            phot_mag_cor = np.array([mfuv_cor, mnuv_cor, m475_cor, m814_cor, m160_cor])
            phot_flux =10**(phot_mag_cor / (-2.5))*3631*1e-23
        
            fig = plt.figure()
            plt.suptitle(catalog['ID'][j])

            # panel showing HST photometry and 0.25 - 1.8 micron model stellar populations
            ax = fig.add_subplot(2,2,1)
            plt.xlim(0.2500,1.8000)
            #plt.ylim(2e-28,5e-27)
            #plt.xscale("log", nonposx='clip')
            #plt.yscale("log", nonposy='clip')
            plt.ylim(0.,3.)
            plt.xlabel(r'Observed Wavelength [$\mu$m]')
            plt.ylabel(r'Flux [10$^{-27}$ erg/s/cm$^2$/Hz]')
            ax.xaxis.set_label_position('top')
        
            plt.plot(wave_one*(1+catalog['z'][j])/1e4,flux_model_one/1e-27,color='red', linewidth=0.1)
            plt.plot(wave_two*(1+catalog['z'][j])/1e4,flux_model_two/1e-27,color='orange', linewidth=0.1)
            plt.plot(wave_thr*(1+catalog['z'][j])/1e4,flux_model_thr/1e-27,color='green', linewidth=0.1)
            plt.plot(wave_fou*(1+catalog['z'][j])/1e4,flux_model_fou/1e-27,color='blue', linewidth=0.1)

        
            plt.scatter(phot_wave/1e4, mags_one_flux/1e-27, marker='o', color='red')
            plt.scatter(phot_wave/1e4, mags_two_flux/1e-27, marker='v', color='orange')
            plt.scatter(phot_wave/1e4, mags_thr_flux/1e-27, marker='s', color='green')
            plt.scatter(phot_wave/1e4, mags_fou_flux/1e-27, marker='<', color='blue')
        
            plt.scatter(phot_wave/1e4, phot_flux/1e-27, color='black', marker='*')

            # panel showing GALEX and F475W photometry and 0.1 - 0.55 micron model stellar populations
            ax = fig.add_subplot(2,2,3)
            ax.set_xlim([1000,5500])
            #ax.set_ylim([1e-29,5e-27])
            #ax.set_xscale("log", nonposx='clip')
            #ax.set_yscale("log", nonposy='clip')
            plt.ylim(0.,3)
            ax.set_xlabel(r'Observed Wavelength [$\AA$]')
            ax.set_ylabel(r'Flux [10$^{-27}$ erg/s/cm$^2$/Hz]')

            plt.plot(wave_one*(1+catalog['z'][j]),flux_model_one/1e-27,color='red', linewidth=0.1)
            plt.plot(wave_two*(1+catalog['z'][j]),flux_model_two/1e-27,color='orange', linewidth=0.1)
            plt.plot(wave_thr*(1+catalog['z'][j]),flux_model_thr/1e-27,color='green', linewidth=0.1)
            plt.plot(wave_fou*(1+catalog['z'][j]),flux_model_fou/1e-27,color='blue', linewidth=0.1)
        
            plt.scatter(phot_wave, mags_one_flux/1e-27, marker='o', color='red')
            plt.scatter(phot_wave, mags_two_flux/1e-27, marker='v', color='orange')
            plt.scatter(phot_wave, mags_thr_flux/1e-27, marker='s', color='green')
            plt.scatter(phot_wave, mags_fou_flux/1e-27, marker='<', color='blue')

            plt.scatter(phot_wave, phot_flux/1e-27, color='black', marker='*')
        
            # make optical color-color plot
            ax = fig.add_subplot(2,2,2)
            ax.set_xlabel('[F814W] - [F160W]')
            ax.set_ylabel('[F475W] - [F814W]')
            ax.set_xlim([-1.0,1.0])
            ax.set_ylim([-1.0,1.0])
            ax.xaxis.set_label_position('top')
            ax.yaxis.set_label_position('right')

            # plot data point with error bar
            plt.scatter(phot_mag_cor[3]-phot_mag_cor[4], phot_mag_cor[2]-phot_mag_cor[3], color='black', marker='*', label='data')
            plt.errorbar(phot_mag_cor[3] - phot_mag_cor[4], phot_mag_cor[2] - phot_mag_cor[3], xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='black',linewidth=0.5)
            plt.scatter(phot_mag[3]-phot_mag[4], phot_mag[2]-phot_mag[3], color='black', marker='.')

            # plot points for the three models
            plt.scatter(mags_one[3]-mags_one[4], mags_one[2]-mags_one[3], color='red', marker='o', label=str(age[0])+' Myr')
            plt.scatter(mags_two[3]-mags_two[4], mags_two[2]-mags_two[3], color='orange', marker='v', label=str(age[1])+' Myr')
            plt.scatter(mags_thr[3]-mags_thr[4], mags_thr[2]-mags_thr[3], color='green', marker='s', label=str(age[2])+' Myr')
            plt.scatter(mags_fou[3]-mags_fou[4], mags_fou[2]-mags_fou[3], color='blue', marker='<', label=str(age[3])+' Myr')

        

        
            # GALEX color-color panel
            ax = fig.add_subplot(2,2,4)
            ax.set_xlabel('[NUV] - [F475W]')
            ax.set_ylabel('[FUV] - [NUV]')
            ax.set_xlim([-0.5,3.0])
            ax.set_ylim([-0.5,3.0])
            ax.yaxis.set_label_position('right')

        
            # plot data point with error bar        
            #ax.errorbar(phot_mag[1] - phot_mag[2], phot_mag[0] - phot_mag[1], xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='magenta')
            plt.scatter(phot_mag_cor[1]-phot_mag_cor[2], phot_mag_cor[0]-phot_mag_cor[1], color='black', marker='*', label='data')
            plt.errorbar(phot_mag_cor[1]-phot_mag_cor[2], phot_mag_cor[0]-phot_mag_cor[1], xerr = 0.1, yerr=0.3, color='black', linewidth=0.5)
            plt.scatter(phot_mag[1]-phot_mag[2], phot_mag[0]-phot_mag[1], color='black', marker='.')
        
            # plot data point after galactic extinction correction
            #mnuv_cor = mnuv[j] - ebv[j]*7.06
            #mfuv_cor = mfuv[j] - ebv[j]*4.37
            #ax.errorbar(mnuv_cor - m475_cor, mfuv_cor- mnuv_cor, xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='purple')

            # plot points for the three models
            plt.scatter(mags_one[1]-mags_one[2], mags_one[0]-mags_one[1], color='red', marker='o', label=str(age[0])+' Myr')
            plt.scatter(mags_two[1]-mags_two[2], mags_two[0]-mags_two[1], color='orange', marker='v', label=str(age[1])+' Myr')
            plt.scatter(mags_thr[1]-mags_thr[2], mags_thr[0]-mags_thr[1], color='green', marker='s', label=str(age[2])+' Myr')
            plt.scatter(mags_fou[1]-mags_fou[2], mags_fou[0]-mags_fou[1], color='blue', marker='<', label=str(age[3])+' Myr')

            plt.legend(loc='lower right', prop={'size': 9})

        
            pdf.savefig()
            plt.close()

os.system('open %s &' % filename)

print('done ', datetime.now()-start)
