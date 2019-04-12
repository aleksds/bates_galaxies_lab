import os
import numpy as np
import fsps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii

sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=2, dust2=0.0)
sp_dust = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=2, dust2=0.4)
sp_must = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=2, dust2=0.8)

bands = np.array(['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])

phottable = 'Photometry/coarse_final/atest.txt'
catalog = ascii.read(phottable)
    
filename = 'fsps_example.pdf'

with PdfPages(filename) as pdf:
    fig = plt.figure()

    t = np.array([0.004,0.006,0.007,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4])
    colors = ['violet', 'indigo', 'blue', 'green', 'cyan', 'orange', 'red', 'darkviolet', 'midnightblue', 'darkgoldenrod', 'chocolate']
    for i in range(0,len(t)):
         wave, spec = sp.get_spectrum(tage=t[i])
         plt.plot(wave,spec,color=colors[i], label=str(int(t[i]*1e3))+' Myr')

    plt.xlim(3000,10000)
    plt.ylim(1e-15,1e-13)
    plt.xscale("log", nonposx='clip')
    plt.yscale("log", nonposy='clip')
    plt.xlabel('Wavelength [Angstroms]')
    plt.ylabel('Luminosity')
    plt.legend(loc='lower right', prop={'size': 10}, framealpha=0.1)
    pdf.savefig()
    plt.close()

    for j in range(0,len(catalog)):

        mags = np.zeros([len(t),len(bands)])
        mags_dust = np.zeros([len(t),len(bands)])
        mags_must = np.zeros([len(t),len(bands)])
        for i in range(0,len(mags)):
            mags[i] = sp.get_mags(tage=t[i], bands=bands, redshift=catalog['z'][j])
            mags_dust[i] = sp_dust.get_mags(tage=t[i], bands=bands, redshift=catalog['z'][j])
            mags_must[i] = sp_must.get_mags(tage=t[i], bands=bands, redshift=catalog['z'][j])
        fig = plt.figure()

        m814 = -2.5*np.log10(catalog['f_814'][j]*1e-9)
        mass = 10**((m814 - mags[:,1])/(-2.5))
        mass_dust = 10**((m814 - mags_dust[:,1])/(-2.5))
        mass_must = 10**((m814 - mags_must[:,1])/(-2.5))
        
        plt.scatter(mags[:,1]-mags[:,2], mags[:,0]-mags[:,1], color='blue', marker='o', label='no dust')
        for i in range(0,len(mags)):
            plt.text(mags[i,1]-mags[i,2], mags[i,0]-mags[i,1],str(int(t[i]*1e3)))
        plt.scatter(mags_dust[:,1]-mags_dust[:,2], mags_dust[:,0]-mags_dust[:,1], color='green', marker='v', label='0.4 mag')
        for i in range(0, len(mags_dust)):
            plt.text(mags_dust[i,1]-mags_dust[i,2], mags_dust[i,0]-mags_dust[i,1],str(int(t[i]*1e3)))
        plt.scatter(mags_must[:,1]-mags_must[:,2], mags_must[:,0]-mags_must[:,1], color='red', marker='s', label='0.8 mag')
        for i in range(0, len(mags_must)):
            plt.text(mags_must[i,1]-mags_must[i,2], mags_must[i,0]-mags_must[i,1],str(int(t[i]*1e3)))
        plt.scatter(-2.5*np.log10(catalog['f_814'][j]*1e-9) - -2.5*np.log10(catalog['f_1600'][j]*1e-9), -2.5*np.log10(catalog['f_475'][j]*1e-9) - -2.5*np.log10(catalog['f_814'][j]*1e-9), color='magenta', marker='*', s=500)
        plt.xlabel('[F814W] - [F160W]')
        plt.ylabel('[F475W] - [F814W]')
        plt.xlim([-0.4,1.4])
        plt.ylim([-0.4,1.4])
        plt.title(catalog['ID'][j])
        plt.legend(loc='lower right', prop={'size': 10})
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        plt.scatter(mass,t*1e3,color='blue', marker='o', label='no dust')
        plt.scatter(mass_dust,t*1e3,color='green', marker='v', label='0.4 mag')
        plt.scatter(mass_must,t*1e3,color='red', marker='s', label='0.8 mag')
        plt.xlabel('SSP Stellar Mass')
        plt.ylabel('SSP Age [Myr]')
        plt.title(catalog['ID'][j])
        plt.ylim([0,120])
        plt.xlim([1e9,1e11])
        plt.xscale("log", nonposx='clip')
        plt.legend(loc='upper left', prop={'size': 10})
        pdf.savefig()
        plt.close()
    
os.system('open %s &' % filename)
