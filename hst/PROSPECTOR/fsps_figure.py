import os
import numpy as np
import fsps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii

#from matplotlib import rc
#rc('text', usetex=True)

sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=0, dust2=0.0)
sp_dust = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                sfh=0, logzsol=0.0, dust_type=0, dust2=1.0)
#sp_must = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
#                                sfh=0, logzsol=0.0, dust_type=0, dust2=0.8)

bands = np.array(['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])

data = ascii.read('../autogalfit/bgl_phot.dat')

filename = 'fsps_figure.pdf'

with PdfPages(filename) as pdf:

    #t = np.array([0.005,0.008,0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.2,0.3])
    #t = np.array([0.010,0.020,0.023,0.025,0.030,0.050,0.100,0.150,0.200, 0.300])
    #colors = ['violet', 'indigo', 'blue', 'green', 'cyan', 'orange', 'red', 'darkviolet', 'midnightblue', 'darkgoldenrod', 'chocolate']
    #colors = ['darkviolet', 'violet', 'purple', 'indigo', 'midnightblue', 'blue', 'cyan', 'green', 'orange', 'red']
    #colors = ['cyan', 'indigo', 'violet', 'darkviolet', 'midnightblue', 'blue', 'purple', 'green', 'orange', 'red']
    t = np.array([0.005,0.010,0.025,0.050,0.100,0.200])
    colors = ['violet', 'magenta', 'blue', 'green', 'orange', 'red']
    #markers = ['.', 'o', 'v', '8', 's', 'p', '*', 'h', 'x', 'D']
    markers = ['o', 'v', 'p', 's', 'h', 'D']

    #fig = plt.figure()

    #for i in range(0,len(t)):
    #     wave, spec = sp.get_spectrum(tage=t[i])
    #     plt.plot(wave,spec,color=colors[i], label=str(int(t[i]*1e3))+' Myr')

    #plt.xlim(3000,10000)
    #plt.ylim(1e-15,1e-13)
    #plt.xscale("log", nonposx='clip')
    #plt.yscale("log", nonposy='clip')
    #plt.xlabel('Wavelength [Angstroms]')
    #plt.ylabel('Luminosity')
    #plt.legend(loc='lower right', prop={'size': 10}, framealpha=0.1)
    #pdf.savefig()
    #plt.close()


    fig = plt.figure()
    
    for j in range(0,len(data)):

        mags = np.zeros([len(t),len(bands)])
        mags_dust = np.zeros([len(t),len(bands)])
        #mags_must = np.zeros([len(t),len(bands)])
        for i in range(0,len(mags)):
            mags[i] = sp.get_mags(tage=t[i], bands=bands, redshift=data['z'][j])
            mags_dust[i] = sp_dust.get_mags(tage=t[i], bands=bands, redshift=data['z'][j])
            #mags_must[i] = sp_must.get_mags(tage=t[i], bands=bands, redshift=data['z'][j])
        #fig = plt.figure()
 
        ax = fig.add_subplot(3,4,j+1)
        plt.rc('xtick',labelsize=8)
        plt.rc('ytick',labelsize=8)
        
        m475 = data['m475'][j] - data['ebv'][j] * 3.248
        m814 = data['m814'][j] - data['ebv'][j] * 1.536
        m160 = data['m160'][j] - data['ebv'][j] * 0.512
        print(data['Galaxy'][j], m475, m814, m160)

        oir_color_unc = np.sqrt(data['u814'][j]**2+data['u160'][j]**2)
        uvo_color_unc = np.sqrt(data['u475'][j]**2+data['u814'][j]**2)
        
        #m814 = -2.5*np.log10(catalog['f_814'][j]*1e-9)
        mass = 10**((m814 - mags[:,1])/(-2.5))
        mass_dust = 10**((m814 - mags_dust[:,1])/(-2.5))
        #mass_must = 10**((m814 - mags_must[:,1])/(-2.5))

        ax.scatter(m814 - m160, m475 - m814, label='data', marker='+', s=20, color='black')
        ax.errorbar(x=m814 - m160, y=m475 - m814, xerr=oir_color_unc, yerr=uvo_color_unc, color='black', elinewidth=2)
        #ax.scatter(mags[:,1]-mags[:,2], mags[:,0]-mags[:,1], color='blue', marker='o', label='model', s=10)
        for i in range(0,len(mags)):
            #ax.text(mags[i,1]-mags[i,2], mags[i,0]-mags[i,1],str(int(t[i]*1e3)), fontsize=5)
            ax.scatter(mags[i,1]-mags[i,2], mags[i,0]-mags[i,1], color=colors[i], marker=markers[i], label=str(int(t[i]*1e3))+' Myr', s=10)
        #ax.scatter(mags_dust[:,1]-mags_dust[:,2], mags_dust[:,0]-mags_dust[:,1], color='green', marker='v', label='0.4 mag')
        #for i in range(0, len(mags_dust)):
        #    plt.text(mags_dust[i,1]-mags_dust[i,2], mags_dust[i,0]-mags_dust[i,1],str(int(t[i]*1e3)))
        oir_color = mags[3,1]-mags[3,2]
        uvo_color = mags[3,0]-mags[3,1]
        oir_color_dust = mags_dust[3,1]-mags_dust[3,2]
        uvo_color_dust = mags_dust[3,0]-mags_dust[3,1]
        #plt.scatter(mags_must[:,1]-mags_must[:,2], mags_must[:,0]-mags_must[:,1], color='red', marker='s', label='0.8 mag')
        #for i in range(0, len(mags_must)):
        #    plt.text(mags_must[i,1]-mags_must[i,2], mags_must[i,0]-mags_must[i,1],str(int(t[i]*1e3)))
        #plt.scatter(-2.5*np.log10(catalog['f_814'][j]*1e-9) - -2.5*np.log10(catalog['f_1600'][j]*1e-9), -2.5*np.log10(catalog['f_475'][j]*1e-9) - -2.5*np.log10(catalog['f_814'][j]*1e-9), color='magenta', marker='*', s=500)
        #ax.scatter(m814 - m160, m475 - m814, color='magenta', marker='*', s=20)

        if j==0:
            ax.set_ylabel('[F475W] - [F814W]', fontsize=8)
            ax.arrow(0.5, 0.5, oir_color_dust - oir_color, uvo_color_dust - uvo_color, head_width=0.05, head_length=0.1, length_includes_head=1, fc='k', ec='k')
            ax.text(0.6,0.7,r"$\hat{\tau}=1$", fontsize=8, rotation=55)
            #ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1.0, 2.5))
        if j==1:
            ax.legend(loc='upper center', fontsize=6, bbox_to_anchor=(0.75, 1.05,0.7,0.2), ncol=7, fancybox=True)
            ax.set_yticklabels(['','','',''])
        #if j==4:
        #    ax.set_ylabel('[F475W] - [F814W]', fontsize=8)
        #if j==8:
        #    ax.set_ylabel('[F475W] - [F814W]', fontsize=8)
        ylabel = [0,4,8]
        if j in ylabel:
            ax.set_ylabel('[F475W] - [F814W]', fontsize=8)
        else:
            ax.set_yticklabels(['','','',''])
        ax.set_xlabel('[F814W] - [F160W]', fontsize=8)
        #plt.xlim([-0.4,1.4])
        ax.set_xlim([-0.8,1.3])
        ax.set_ylim([-0.2,1.6])
        #plt.title(data['Galaxy'][j])
        ax.text(0.5,1.3,data['Galaxy'][j], fontsize=8)
        #ax.legend(loc='lower right', prop={'size': 10})
        #extraticks = [-0.5, 0.5]
        #ax.set_xticks(list(ax.get_xticks()) + extraticks)
    pdf.savefig()
    plt.close()

        #fig = plt.figure()
        #plt.scatter(mass,t*1e3,color='blue', marker='o', label='no dust')
        #plt.scatter(mass_dust,t*1e3,color='green', marker='v', label='0.4 mag')
        #plt.scatter(mass_must,t*1e3,color='red', marker='s', label='0.8 mag')
        #plt.xlabel('SSP Stellar Mass')
        #plt.ylabel('SSP Age [Myr]')
        #plt.title(data['Galaxy'][j])
        #plt.ylim([0,120])
        #plt.xlim([1e9,1e11])
        #plt.xscale("log", nonposx='clip')
        #plt.legend(loc='upper left', prop={'size': 10})
        #pdf.savefig()
        #plt.close()
    
os.system('open %s &' % filename)
