import os
import numpy as np
import fsps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii

sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=2, dust2=0.0)
sp_dust = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=2, dust2=0.8)
sp_must = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=1, dust2=0.8)

print('done with fsps.StellarPopulation')

#bands = np.array(['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])
bands = np.array(['galex_fuv','galex_nuv','wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])

phottable = 'Photometry/coarse_final/atest.txt'
catalog = ascii.read(phottable)

ad_values = ascii.read('../autogalfit/ad_mag_size_table.dat')
mfuv = np.array(ad_values['mfuv'])
mnuv = np.array(ad_values['mnuv'])
ebv = np.array(ad_values['ebv'])

filename = 'fsps_explore.pdf'

with PdfPages(filename) as pdf:
    fig = plt.figure()
    # ages from 4 Myr to 30 Myr
    t = np.array([0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.015,0.02,0.023,0.025,0.03])
    colors = ['magenta','violet','indigo','blue','green','cyan','orange','red','magenta','violet', 'indigo', 'blue', 'green', 'cyan', 'orange', 'red']

    # plot spectra associated with each model
    for i in range(0,len(t)):
         wave, spec = sp.get_spectrum(tage=t[i])
         plt.plot(wave,spec,color=colors[i], label=str(int(t[i]*1e3))+' Myr')
         print('done with get_spectrum for tage=', t[i])

    plt.xlim(900,30000)
    plt.ylim(1e-15,1e-12)
    plt.xscale("log", nonposx='clip')
    plt.yscale("log", nonposy='clip')
    plt.xlabel('Wavelength [Angstroms]')
    plt.ylabel('Luminosity')
    plt.legend(loc='lower right', prop={'size': 10})
    pdf.savefig()
    plt.close()

    print('done making spectrum plots')
    
    for j in range(0,len(catalog)):

        mags = np.zeros([len(t),len(bands)])
        mags_dust = np.zeros([len(t),len(bands)])
        mags_must = np.zeros([len(t),len(bands)])

        for i in range(0,len(mags)):
            mags[i] = sp.get_mags(tage=t[i], bands=bands, redshift=catalog['z'][j])
            mags_dust[i] = sp_dust.get_mags(tage=t[i], bands=bands, redshift=catalog['z'][j])
            mags_must[i] = sp_must.get_mags(tage=t[i], bands=bands, redshift=catalog['z'][j])
            print('done with get_mags for ', catalog['ID'][j], t[i])
        fig = plt.figure()

        m814 = -2.5*np.log10(catalog['f_814'][j]*1e-9)
        mass = 10**((m814 - mags[:,3])/(-2.5))
        mass_dust = 10**((m814 - mags_dust[:,1])/(-2.5))
        mass_must = 10**((m814 - mags_must[:,3])/(-2.5))
        print(catalog['ID'][j])
        print(t)
        print(mass)
        print(mass_dust)
        print(mass_must)

        # make optical color-color plot
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('[F814W] - [F160W]')
        ax.set_ylabel('[F475W] - [F814W]')
        ax.set_xlim([-1.0,2.0])
        ax.set_ylim([-1.0,2.0])

        # plot points for the three models
        plt.scatter(mags[:,3]-mags[:,4], mags[:,2]-mags[:,3], color='blue', marker='o', label='no dust')
        plt.scatter(mags_dust[:,3]-mags_dust[:,4], mags_dust[:,2]-mags_dust[:,3], color='green', marker='v', label='CZ: 0.8 mag')
        plt.scatter(mags_must[:,3]-mags_must[:,4], mags_must[:,2]-mags_must[:,3], color='red', marker='s', label='PL: 0.8 mag')

        # add time labels for each model
        for i in range(0,len(mags)):
            plt.text(mags[i,3]-mags[i,4], mags[i,2]-mags[i,3],str(int(t[i]*1e3)))
            plt.text(mags_dust[i,3]-mags_dust[i,4], mags_dust[i,2]-mags_dust[i,3],str(int(t[i]*1e3)))
            plt.text(mags_must[i,3]-mags_must[i,4], mags_must[i,2]-mags_must[i,3],str(int(t[i]*1e3)))

        # plot the data
        plt.scatter(-2.5*np.log10(catalog['f_814'][j]*1e-9) - -2.5*np.log10(catalog['f_1600'][j]*1e-9), -2.5*np.log10(catalog['f_475'][j]*1e-9) - -2.5*np.log10(catalog['f_814'][j]*1e-9), color='magenta', marker='*', s=500)

        # apply galactic extinction correction and re-plot
        m814 = -2.5*np.log10(catalog['f_814'][j]*1e-9)
        m475 = -2.5*np.log10(catalog['f_475'][j]*1e-9)
        m160 = -2.5*np.log10(catalog['f_1600'][j]*1e-9)
        m814_cor = m814 - ebv[j]*1.536
        m475_cor = m475 - ebv[j]*3.248
        m160_cor = m160 - ebv[j]*0.512
        plt.scatter(m814_cor - m160_cor, m475_cor - m814_cor, color='purple', marker='*', s=500)
        
        plt.suptitle(catalog['ID'][j])
        ax.legend(loc='lower right', prop={'size': 10})
        pdf.savefig()
        plt.close()

        # GALEX color-color panel
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('[NUV] - [F475W]')
        ax.set_ylabel('[FUV] - [NUV]')
        ax.set_xlim([-0.5,3.0])
        ax.set_ylim([-0.5,3.0])
        plt.suptitle(catalog['ID'][j])

        # plot data point with error bar        
        ax.errorbar(mnuv[j] - -2.5*np.log10(catalog['f_475'][j]*1e-9), mfuv[j] - mnuv[j], xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='magenta')

        # plot data point after galactic extinction correction
        mnuv_cor = mnuv[j] - ebv[j]*7.06
        mfuv_cor = mfuv[j] - ebv[j]*4.37
        ax.errorbar(mnuv_cor - m475_cor, mfuv_cor- mnuv_cor, xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='purple')

        # plot points for the three models
        ax.scatter(mags[:,1]-mags[:,2], mags[:,0]-mags[:,1], color='blue', marker='o', label='no dust')
        ax.scatter(mags_must[:,1]-mags_must[:,2], mags_must[:,0]-mags_must[:,1], color='red', marker='s', label='PL: 0.8 mag')
        ax.scatter(mags_dust[:,1]-mags_dust[:,2], mags_dust[:,0]-mags_dust[:,1], color='green', marker='s', label='CZ: 0.8 mag')

        # add time labels for each model
        for i in range(0,len(mags)):
            plt.text(mags[i,1]-mags[i,2], mags[i,0]-mags[i,1],str(int(t[i]*1e3)))
            plt.text(mags_dust[i,1]-mags_dust[i,2], mags_dust[i,0]-mags_dust[i,1],str(int(t[i]*1e3)))
            plt.text(mags_must[i,1]-mags_must[i,2], mags_must[i,0]-mags_must[i,1],str(int(t[i]*1e3)))
        
        pdf.savefig()
        plt.close()

os.system('open %s &' % filename)
