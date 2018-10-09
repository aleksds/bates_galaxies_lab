import os
import numpy as np
import fsps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.constants import G
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

bands = np.array(['wfc3_uvis_f475w','wfc3_uvis_f814w','wfc3_ir_f160w'])

phottable = 'Photometry/coarse_final/atest.txt'
catalog = ascii.read(phottable)

sdss_files = ['sdss/spec-0761-54524-0409.fits', 'sdss/spec-0566-52238-0319.fits', 'sdss/spec-0483-51924-0628.fits', 'sdss/spec-1305-52757-0191.fits', 'sdss/spec-0581-52356-0196.fits', 'sdss/spec-0518-52282-0605.fits', 'sdss/spec-0913-52433-0300.fits', 'sdss/spec-6712-56421-0114.fits', 'sdss/spec-1054-52516-0153.fits', 'sdss/spec-1577-53495-0564.fits', 'sdss/spec-0639-52146-0388.fits', 'sdss/spec-0732-52221-0445.fits']
sdss_norm = np.zeros(len(catalog))+2e-21
sdss_norm[2] = 3e-21
sdss_norm[5] = 1e-21
sdss_norm[10] = 3e-21

filename = 'fsps_3models.pdf'

vflow = np.array([1228,1206,2470,1778,1828,1830,875,1211,829,2416,1456,606])

jc_values = ascii.read('../autogalfit/jc_mag_size_table.dat')
re_pix = np.array(jc_values['re'])
z = catalog['z']
re_kpc = np.zeros(len(z))

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
for i in range(0,len(z)):
    re_kpc[i] = re_pix[i] / cosmo.arcsec_per_kpc_proper(z[i]) * 0.025 * u.arcsec / u.kpc 
    #print(re_pix[i], re_kpc[i])


#                J0826,  J0901,  J0905,  J0944,  J1107, J1219, J1341, J1506, J1558, J1613, J2116, J2140
#age  = np.array([0.006, 0.0078, 0.0066, 0.0054, 0.0052, 0.006, 0.005, 0.006, 0.006, 0.008, 0.007, 0.006])
#dust = np.array([1.1,   0.38,   0.51,   0.54,   1.5,    1.6,   1.3,   0.68,  1.2,   1.23,  0.7,   1.7])
#logm = np.array([9.8,   9.56,   9.51,   9.59,  10.0,   10.1,   9.97,  9.7,   9.8,  10.11,  9.48,  9.81])

age = np.zeros([len(catalog),3])
dust = np.zeros([len(catalog),3])
logm = np.zeros([len(catalog),3])

# J0826
age[0] = np.array([0.035,  0.0079, 0.006])
dust[0] = np.array([0.1,   0.1,   1.0])
logm[0] = np.array([10.1,  9.4,  9.7])

# J0901
age[1] = np.array([0.04,  0.008, 0.0078])
dust[1] = np.array([0.25,  0.2,   0.3])
logm[1] = np.array([10.0,   9.3,  9.3])

# J0905
age[2] = np.array([0.045,  0.0074, 0.0065])
dust[2] = np.array([0.0,   0.15,   0.48])
logm[2] = np.array([10.3,   9.5,  9.6])

# J0944
age[3] = np.array([0.035,  0.008, 0.008])
dust[3] = np.array([0.4,   0.4,  0.4])
logm[3] = np.array([10.1,  9.4,  9.4])

# J1107
age[4] = np.array([0.006,  0.0075, 0.0052])
dust[4] = np.array([0.65,   0.05,   1.05])
logm[4] = np.array([9.3,    9.1,    9.6])

# J1219
age[5] = np.array([0.050,  0.0075, 0.006])
dust[5] = np.array([0.6,   0.7,   1.3])
logm[5] = np.array([10.3,  9.5,   9.7])

# J1341
age[6] = np.array([0.05,  0.007, 0.005])
dust[6] = np.array([0.0,   0.0,   0.9])
logm[6] = np.array([10.4,  9.4,   10.0])

# J1506
age[7] = np.array([0.04,  0.0077, 0.006])
dust[7] = np.array([0.0,   0.0,  0.8])
logm[7] = np.array([10.2,  9.4,  9.7])

# J1558
age[8] = np.array([0.028,  0.0083, 0.006])
dust[8] = np.array([0.1,   0.1,  1.2])
logm[8] = np.array([9.7,  9.1,  9.6])

# J1613
age[9] = np.array([0.04,  0.008, 0.008])
dust[9] = np.array([0.9,   0.8,   0.8])
logm[9] = np.array([10.5,  9.6,   9.6])

# J2116
age[10] = np.array([0.06,  0.007, 0.007])
dust[10] = np.array([0.14,   0.5,  0.5])
logm[10] = np.array([10.4,   9.6,  9.6])

# J2140
age[11] = np.array([0.0064, 0.0072, 0.0071])
dust[11] = np.array([1.0,   0.7,  0.7])
logm[11] = np.array([9.8,   9.7,  9.7])

colors=['red', 'green', 'blue']


with PdfPages(filename) as pdf:

    for m in range(0, 12):#len(catalog)):

        # generate stellar population models
        sp_zero = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                         sfh=0, logzsol=0.0, dust_type=2, dust2=dust[m,0])
        mags_zero = sp_zero.get_mags(tage=age[m,0], bands=bands, redshift=catalog['z'][m])
        wave_zero, spec_zero = sp_zero.get_spectrum(tage=age[m,0])

        sp_one = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                        sfh=0, logzsol=0.0, dust_type=2, dust2=dust[m,1])
        mags_one = sp_one.get_mags(tage=age[m,1], bands=bands, redshift=catalog['z'][m])
        wave_one, spec_one = sp_one.get_spectrum(tage=age[m,1])

        sp_two = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                        sfh=0, logzsol=0.0, dust_type=2, dust2=dust[m,2])        
        mags_two = sp_two.get_mags(tage=age[m,2], bands=bands, redshift=catalog['z'][m])
        wave_two, spec_two = sp_two.get_spectrum(tage=age[m,2])

        # make a figure
        fig = plt.figure()
        plt.suptitle(catalog['ID'][m])

        # fnu panel
        ax = fig.add_subplot(2,2,1)
        ax.set_xlim(2000,20000)
        ax.set_ylim(1e-15,1e-12)
        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')
        #ax.set_xlabel('Wavelength [Angstroms]')
        ax.set_ylabel('Flux (f_nu)')
        #ax.set_xticklabels(())
        ax.set_yticklabels(())
        
        ax.plot(wave_zero,spec_zero,color='red', label=str(round(age[m,0]*1e3,1))+' Myr, dust '+str(dust[m,0])+' mag')
        ax.plot(wave_one,spec_one,color='green', label=str(round(age[m,1]*1e3,1))+' Myr, dust '+str(dust[m,1])+' mag')
        ax.plot(wave_two,spec_two,color='blue', label=str(round(age[m,2]*1e3,1))+' Myr, dust '+str(dust[m,2])+' mag')    

        # flambda panel
        ax = fig.add_subplot(2,2,2)
        ax.set_xlim(2500,5000)
        ax.set_ylim(1e-22,1e-19)
        #plt.xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')
        ax.set_xlabel('Wavelength [Angstroms]')
        ax.set_ylabel('Flux (f_lambda)')
        ax.set_yticklabels(())

        ax.plot(wave_zero,spec_zero/wave_zero**2,color='red', label=str(round(age[m,0]*1e3,1))+' Myr, '+str(dust[m,0])+' mag')
        ax.plot(wave_one,spec_one/wave_one**2,color='green', label=str(round(age[m,1]*1e3,1))+' Myr, '+str(dust[m,1])+' mag')
        ax.plot(wave_two,spec_two/wave_two**2,color='blue', label=str(round(age[m,2]*1e3,1))+' Myr, '+str(dust[m,2])+' mag')    

        #sdss_file = 'sdss/spec-0761-54524-0409.fits'
        hdu_sdss = fits.open(sdss_files[m])
        sdss_flux = hdu_sdss[1].data['flux']
        sdss_wave = 10.**(hdu_sdss[1].data['loglam'])
        
        ax.plot(sdss_wave/(1.+catalog['z'][m]), sdss_flux * sdss_norm[m], color='black', lw=0.1)
        
        # color-color panel
        ax = fig.add_subplot(2,2,3)
        ax.set_xlabel('[F814W] - [F160W]')
        ax.set_ylabel('[F475W] - [F814W]')
        #ax.set_xlim([0.0,0.8])
        #ax.set_ylim([0.0,0.8])
        ax.set_xlim([-0.5,1.5])
        ax.set_ylim([-0.5,1.5])
        #ax.scatter(-2.5*np.log10(catalog['f_814'][m]*1e-9) - -2.5*np.log10(catalog['f_1600'][m]*1e-9), -2.5*np.log10(catalog['f_475'][m]*1e-9) - -2.5*np.log10(catalog['f_814'][m]*1e-9), color='magenta', marker='*', s=500)
        ax.errorbar(-2.5*np.log10(catalog['f_814'][m]*1e-9) - -2.5*np.log10(catalog['f_1600'][m]*1e-9), -2.5*np.log10(catalog['f_475'][m]*1e-9) - -2.5*np.log10(catalog['f_814'][m]*1e-9), xerr = 0.05*np.sqrt(2), yerr=0.05*np.sqrt(2), color='magenta')

        ax.scatter(mags_zero[1]-mags_zero[2], mags_zero[0]-mags_zero[1], color='red', marker='s',label=str(round(age[m,0]*1e3,1))+' Myr, '+str(dust[m,0])+' mag, '+str(round(logm[m,0],2)))
        ax.scatter(mags_one[1]-mags_one[2], mags_one[0]-mags_one[1], color='green', marker='v',label=str(round(age[m,1]*1e3,1))+' Myr, '+str(dust[m,1])+' mag, '+str(round(logm[m,1],2)))
        ax.scatter(mags_two[1]-mags_two[2], mags_two[0]-mags_two[1], color='blue', marker='o',label=str(round(age[m,2]*1e3,1))+' Myr, '+str(dust[m,2])+' mag, '+str(round(logm[m,2],2)))
        
        ax.legend(loc='lower right', prop={'size': 8})

        # check masses
        m814 = -2.5*np.log10(catalog['f_814'][m]*1e-9)
        mass_zero = 10**((m814 - mags_zero[1])/(-2.5))
        mass_one = 10**((m814 - mags_one[1])/(-2.5))
        mass_two = 10**((m814 - mags_two[1])/(-2.5))
        print(catalog['ID'][m], np.log10(mass_zero), np.log10(mass_one), np.log10(mass_two))
        masses = np.array([mass_zero,mass_one,mass_two])
        print('factor of ', np.max(masses)/np.min(masses))      

        vesc_zero = np.sqrt(G*mass_zero*const.M_sun/(re_kpc[m]*u.kpc)).to('km/s')
        vesc_one = np.sqrt(G*mass_one*const.M_sun/(re_kpc[m]*u.kpc)).to('km/s')
        vesc_two = np.sqrt(G*mass_two*const.M_sun/(re_kpc[m]*u.kpc)).to('km/s')
        print(vflow[m], vesc_zero, vesc_one, vesc_two)

        plt.text(1.8,1.0,'vout='+str(vflow[m]))
        plt.text(1.8,0.7,'vesc_red='+str(int(vesc_zero*u.s/u.km)))
        plt.text(1.8,0.4,'vesc_green='+str(int(vesc_one*u.s/u.km)))
        plt.text(1.8,0.1,'vesc_blue='+str(int(vesc_two*u.s/u.km)))
                
        pdf.savefig()
        plt.close()
    
os.system('open %s &' % filename)
