import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from matplotlib.ticker import ScalarFormatter
from sedpy import observate
from astropy.io import fits


projpath = os.getcwd()+'/Brown2014/An_Atlas_of_Galaxy_SEDs/An_Atlas_of_Galaxy_SEDs/'

files = glob.glob(projpath+'*.dat')

table = fits.open('hizea_photo_galex_wise_v1.0.fit')
ngal = len(table[1].data)
w3_mag = table[1].data['AB_MAG'][0:ngal][:,9]
w3_mag_err = table[1].data['AB_MAG_ERR'][0:ngal][:,9]
w4_mag = table[1].data['AB_MAG'][0:ngal][:,10]
w4_mag_err = table[1].data['AB_MAG_ERR'][0:ngal][:,10]
redshift = table[1].data['Z']
names = table[1].data['SHORT_NAME'][0:ngal]

w3minusw4_unc = np.sqrt(w3_mag_err ** 2 + w4_mag_err ** 2)
valid = w3minusw4_unc<0.5
print(valid)


def plot_phot():
    ax.scatter(redshift[valid], w3_mag[valid] - w4_mag[valid], marker='*')
    ax.errorbar(redshift[valid], w3_mag[valid] - w4_mag[valid], yerr=w3minusw4_unc[valid], color='purple', linestyle="None")

    for i in range(0,len(redshift[valid])):
        plt.text(redshift[valid][i], w3_mag[valid][i]-w4_mag[valid][i], names[valid][i], fontsize=6)


name = 'brown_seds.pdf'

with PdfPages(name) as pdf:
    for i in range(0,len(files)):
    #for i in range(0,1):


        table = ascii.read(files[i])

        wave = np.array(table['col1'])
        flam = np.array(table['col2'])
        flux = wave * flam

        fnu = flam * wave**2 / 3e18 #* 1e46
        
        loc = files[0].find('spec')

        gal = files[i][loc-9:loc-1]
        #print(gal)

        filternames = ['wise_w{}'.format(n) for n in ['1', '2', '3', '4']]
        wise_filters = observate.load_filters(filternames)
        zarray = np.arange(41)/100.+0.4
        w3_minus_w4 = np.zeros(len(zarray))
        for j in range(0,len(zarray)):
            model_mags = observate.getSED(wave*(1+zarray[j]), flam, filterlist=wise_filters)
            w3_minus_w4[j] = model_mags[2] - model_mags[3]
            #print('model_mags:', model_mags)
            #print('[W3]-[W4]:', w3_minus_w4[j])

        
        top = np.max(w3_minus_w4)
        bot = np.min(w3_minus_w4)

        if (top > 1.2):

            print(gal)
            
            ## lambda * f_lambda figure
            #fig = plt.figure()
            #ax = fig.add_subplot(111)
            #
            #plt.plot(wave/1e4, flux)
            #
            #plt.title(gal)
            #
            #prange = wave > 3e4
            #fmin = np.min(flux[prange])
            #fmax = np.max(flux[prange])
            #
            #plt.yscale('log')
            #plt.ylim([fmin/2,fmax*2])
            #
            #plt.xscale('log')
            #plt.xlim(3,40)
            #
            #plt.xlabel(r'Wavelength [$\mu$m]')
            #plt.ylabel(r'Flux $\lambda f_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$]')
            #
            #ax.set_xticks([3, 5, 8, 10, 20, 30])
            #ax.xaxis.set_major_formatter(ScalarFormatter())
            #
            #labels = [item.get_text() for item in ax.get_xticklabels()]
            #labels[0] = '3'
            #labels[1] = '5'
            #labels[2] = '8'
            #labels[3] = '10'
            #labels[4] = '20'
            #labels[5] = '30'
            #
            #pdf.savefig()
            #plt.close()
            
            # 2nd plot -- fnu figure
            fig = plt.figure()
            ax = fig.add_subplot(211)
            
            plt.plot(wave/1e4, fnu)
            
            plt.title(gal)
            
            prange = wave > 3e4
            fmin = np.min(fnu[prange])
            fmax = np.max(fnu[prange])
            
            plt.yscale('log')
            plt.ylim([fmin/2,fmax*2])
            
            plt.xscale('log')
            plt.xlim(3,40)
            
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel(r'Flux $f_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
            
            ax.set_xticks([3, 5, 8, 10, 20, 30])
            ax.xaxis.set_major_formatter(ScalarFormatter())
            
            labels = [item.get_text() for item in ax.get_xticklabels()]
            labels[0] = '3'
            labels[1] = '5'
            labels[2] = '8'
            labels[3] = '10'
            labels[4] = '20'
            labels[5] = '30'
            
            ax.axvspan(12/1.4, 12/1.8, facecolor='g', alpha=0.5)
            ax.axvspan(22/1.4, 22/1.8, facecolor='r', alpha=0.5)
            
            #pdf.savefig()
            #plt.close()
            
            ## 3rd plot -- flambda figure
            #fig = plt.figure()
            #
            #ax = fig.add_subplot(111)
            #
            #plt.plot(wave/1e4, flam)
            #
            #plt.title(gal)
            #
            #prange = wave > 3e4
            #fmin = np.min(flam[prange])
            #fmax = np.max(flam[prange])
            #
            #plt.yscale('log')
            #plt.ylim([fmin/2,fmax*2])
            #
            #plt.xscale('log')
            #plt.xlim(3,40)
            #
            #plt.xlabel(r'Wavelength [$\mu$m]')
            #plt.ylabel(r'Flux $f_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')
            #
            #ax.set_xticks([3, 5, 8, 10, 20, 30])
            #ax.xaxis.set_major_formatter(ScalarFormatter())
            #
            #labels = [item.get_text() for item in ax.get_xticklabels()]
            #labels[0] = '3'
            #labels[1] = '5'
            #labels[2] = '8'
            #labels[3] = '10'
            #labels[4] = '20'
            #labels[5] = '30'
            #
            #plt.tight_layout()
            #
            #pdf.savefig()
            #plt.close()
            
            # plot of WISE W3-W4 clor
            #fig = plt.figure()
            ax = fig.add_subplot(212)
            
            plt.plot(zarray, w3_minus_w4)
            
            #plt.title(gal)
            
            #plt.ylim([-1.5,2.5])
            plt.ylim([0,2.5])
            plt.xlim([0.4,0.8])
            
            plt.xlabel(r'Redshift')
            plt.ylabel(r'[W3] - [W4] color')

            plot_phot()
            
            #ax.axhline(y=top)
            #ax.axhline(y=bot)
            
            #plt.text(0.45, top+0.15, 'max='+f'{top:.2f}')
            #plt.text(0.45, bot-0.3, 'min='+f'{bot:.2f}')
            
            plt.tight_layout()
            
            pdf.savefig()
            plt.close()
        
        
os.system('open %s &' % name)
