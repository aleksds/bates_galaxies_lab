# Aleks Diamond-Stanic
# 20190404
# Goal: Starting with a suite of suite of output fits files from GALFIT, perform aperture photometry of the data and model to measure total and residual flux

# import relevant packages
import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.io import ascii
from astropy import units as u
from photutils import CircularAperture
from photutils import aperture_photometry
import img_scale
from astropy.cosmology import FlatLambdaCDM 


# calculate magnitude for a given flux in Jy
def mag(flux):
    return -2.5*np.log10(flux/3631)

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

# the code requires the user to specify an input directory that has output fits files from GALFIT
dir = sys.argv[1]

# find all of the relevant fits files
files = glob.glob(dir+'/*output.fits')

# define relevant variables
#xcen = 50
#ycen = 50
xcen = 100
ycen = 100
fmin = -2
fmax = 12
radii = np.arange(100)+1
#radii = np.array([1,2,3,4,5,6,7,8,9,10,15,20,30,40,50])
#radii = np.array([5,10,15,20,30,40,60])
#radii = np.array([5,15,30,45])
#radii = np.array([5,10,20,50])
#radii = np.array([8,10,12,15,25,50])

# arrays for v band
vdat_flux = np.zeros([len(files),len(radii)])
vunc_flux = np.zeros([len(files),len(radii)])
vmod_flux = np.zeros([len(files),len(radii)])
vres_flux = np.zeros([len(files),len(radii)])

vdat_annulus = np.zeros([len(files),len(radii)])
vunc_annulus = np.zeros([len(files),len(radii)])
vmod_annulus = np.zeros([len(files),len(radii)])
vres_annulus = np.zeros([len(files),len(radii)])

# arrays for u band
udat_flux = np.zeros([len(files),len(radii)])
uunc_flux = np.zeros([len(files),len(radii)])
umod_flux = np.zeros([len(files),len(radii)])
ures_flux = np.zeros([len(files),len(radii)])

udat_annulus = np.zeros([len(files),len(radii)])
uunc_annulus = np.zeros([len(files),len(radii)])
umod_annulus = np.zeros([len(files),len(radii)])
ures_annulus = np.zeros([len(files),len(radii)])

# arrays for j band
jdat_flux = np.zeros([len(files),len(radii)])
junc_flux = np.zeros([len(files),len(radii)])
jmod_flux = np.zeros([len(files),len(radii)])
jres_flux = np.zeros([len(files),len(radii)])

jdat_annulus = np.zeros([len(files),len(radii)])
junc_annulus = np.zeros([len(files),len(radii)])
jmod_annulus = np.zeros([len(files),len(radii)])
jres_annulus = np.zeros([len(files),len(radii)])

# array to store information from headers
vfnu = np.zeros(len(files))
vexp = np.zeros(len(files))

ufnu = np.zeros(len(files))
uexp = np.zeros(len(files))

jfnu = np.zeros(len(files))
jexp = np.zeros(len(files))

data = ascii.read('bgl_phot.dat')

# let's produce plots that will be saved in a pdf file
name = 'galfit_image_'+dir+'_'+str(len(radii))+'_'+str(np.max(radii))+'.pdf'

with PdfPages(name) as pdf:
    
    # loop over each output fits file from GALFIT
    for i in range(0, len(files)):
        
        # read in the image
        hdu = fits.open(files[i])
        
        # these images have 18 different extensions:
        # hdu[0] is INPUT_V
        # hdu[1] is INPUT_U
        # hdu[2] is INPUT_J
        # hdu[3] is MODEL_V
        # hdu[4] is MODEL_U
        # hdu[5] is MODEL_J
        # hdu[6] is RESIDUAL_V
        # hdu[7] is RESIDUAL_U
        # hdu[8] is RESIDUAL_J
        # hdu[9] is COMPONENT_1_sersic_V
        # hdu[10] is COMPONENT_1_sersic_U
        # hdu[11] is COMPONENT_1_sersic_J
        # hdu[12] is SIGMA_V
        # hdu[13] is SIMGA_U
        # hdu[14] is SIGMA_J
        # hdu[15] is PSF_V
        # hdu[16] is PSF_U
        # hdu[17] is PSF_J

        # F814W corresponds roughly to rest-frame V-band
        vdat, vdat_head = hdu[0].data-0.53, hdu[0].header
        vfnu[i] = vdat_head['PHOTFNU']
        vexp[i] = vdat_head['EXPTIME']
        vunc, vunc_head = hdu[12].data, hdu[12].header
        vmod, vmod_head = hdu[3].data-0.53, hdu[3].header
        vres, vres_head = hdu[6].data-0.53, hdu[6].header

        # F475W corresponds roughly to rest-frame U-band
        udat, udat_head = hdu[1].data-1.14, hdu[1].header
        ufnu[i] = udat_head['PHOTFNU']
        uexp[i] = udat_head['EXPTIME']
        uunc, uunc_head = hdu[13].data, hdu[13].header
        umod, umod_head = hdu[4].data-1.14, hdu[4].header
        ures, ures_head = hdu[7].data-1.14, hdu[7].header        

        # F160W corresponds roughly to rest-frame J-band
        jdat, jdat_head = hdu[2].data+0.89, hdu[2].header
        jfnu[i] = jdat_head['PHOTFNU']
        jexp[i] = jdat_head['EXPTIME']
        junc, junc_head = hdu[14].data, hdu[14].header
        jmod, jmod_head = hdu[5].data+0.89, hdu[5].header
        jres, jres_head = hdu[8].data+0.89, hdu[8].header  

        redshift = data['z'][i]
        arcsecperkpc = cosmo.arcsec_per_kpc_proper(redshift)
        radpix_15kpc = arcsecperkpc * 12.5 / 0.05 * u.kpc / u.arcsec
        
        # loop over each aperture
        for j in range(0,len(radii)):

            # define the circular aperture
            aperture = CircularAperture([xcen, ycen], radii[j])

            # perform aperture photmetry in that aperture and save the flux value
            vdat_table = aperture_photometry(vdat, aperture)
            vdat_flux[i,j] = vdat_table['aperture_sum'][0]
            udat_table = aperture_photometry(udat, aperture)
            udat_flux[i,j] = udat_table['aperture_sum'][0]
            jdat_table = aperture_photometry(jdat, aperture)
            jdat_flux[i,j] = jdat_table['aperture_sum'][0]
                
            # repeat for the sigma image
            vunc_table = aperture_photometry(vunc, aperture)
            vunc_flux[i,j] = vunc_table['aperture_sum'][0]
            uunc_table = aperture_photometry(uunc, aperture)
            uunc_flux[i,j] = uunc_table['aperture_sum'][0]
            junc_table = aperture_photometry(junc, aperture)
            junc_flux[i,j] = junc_table['aperture_sum'][0]
                        
            # repeat for the model image
            vmod_table = aperture_photometry(vmod, aperture)
            vmod_flux[i,j] = vmod_table['aperture_sum'][0]
            umod_table = aperture_photometry(umod, aperture)
            umod_flux[i,j] = umod_table['aperture_sum'][0]
            jmod_table = aperture_photometry(jmod, aperture)
            jmod_flux[i,j] = jmod_table['aperture_sum'][0]
                        
            # repeat for the residual image
            vres_table = aperture_photometry(vres, aperture)
            vres_flux[i,j] = vres_table['aperture_sum'][0]
            ures_table = aperture_photometry(ures, aperture)
            ures_flux[i,j] = ures_table['aperture_sum'][0]
            jres_table = aperture_photometry(jres, aperture)
            jres_flux[i,j] = jres_table['aperture_sum'][0]

            if j==0:
                vdat_annulus[i,j] = vdat_flux[i,j]
                udat_annulus[i,j] = udat_flux[i,j]
                jdat_annulus[i,j] = jdat_flux[i,j]

                vunc_annulus[i,j] = vunc_flux[i,j]
                uunc_annulus[i,j] = uunc_flux[i,j]
                junc_annulus[i,j] = junc_flux[i,j]
                
                vmod_annulus[i,j] = vmod_flux[i,j]
                umod_annulus[i,j] = umod_flux[i,j]
                jmod_annulus[i,j] = jmod_flux[i,j]
                
                vres_annulus[i,j] = vres_flux[i,j]
                ures_annulus[i,j] = ures_flux[i,j]
                ures_annulus[i,j] = jres_flux[i,j]
            else:
                vdat_annulus[i,j] = vdat_flux[i,j] - vdat_flux[i,j-1]
                udat_annulus[i,j] = udat_flux[i,j] - udat_flux[i,j-1]
                jdat_annulus[i,j] = jdat_flux[i,j] - jdat_flux[i,j-1]

                vunc_annulus[i,j] = vunc_flux[i,j] - vunc_flux[i,j-1]
                uunc_annulus[i,j] = uunc_flux[i,j] - uunc_flux[i,j-1]
                junc_annulus[i,j] = junc_flux[i,j] - junc_flux[i,j-1]
                
                vmod_annulus[i,j] = vmod_flux[i,j] - vmod_flux[i,j-1]
                umod_annulus[i,j] = umod_flux[i,j] - umod_flux[i,j-1]
                jmod_annulus[i,j] = jmod_flux[i,j] - jmod_flux[i,j-1]
                
                vres_annulus[i,j] = vres_flux[i,j] - vres_flux[i,j-1]
                ures_annulus[i,j] = ures_flux[i,j] - ures_flux[i,j-1]
                jres_annulus[i,j] = jres_flux[i,j] - jres_flux[i,j-1]


        fig = plt.figure()

        ax = fig.add_subplot(3,3,1)
        plt.title('F475W Data')
        rotated = np.flip(np.rot90(udat, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,2)
        plt.title('F475W Model')
        rotated = np.flip(np.rot90(umod, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,3)
        plt.title('F475W Resid')
        rotated = np.flip(np.rot90(ures, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,4)
        plt.title('F814W Data')
        rotated = np.flip(np.rot90(vdat, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,5)
        plt.title('F814W Model')
        rotated = np.flip(np.rot90(vmod, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,6)
        plt.title('F814W Resid')
        rotated = np.flip(np.rot90(vres, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,7)
        plt.title('F814W Data')
        rotated = np.flip(np.rot90(jdat, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,8)
        plt.title('F814W Model')
        rotated = np.flip(np.rot90(jmod, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')

        ax = fig.add_subplot(3,3,9)
        plt.title('F814W Resid')
        rotated = np.flip(np.rot90(jres, 2), 1)
        plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
        
        pdf.savefig()
        plt.close()
                
        # start our figure
        fig = plt.figure()

        # have the title at the top of the page be the name of the galfit output 
        plt.suptitle(vdat_head['LOGFILE'])


        #aper = 39
        default = 39
        aper = round(float(radpix_15kpc)-1.)

        print(vdat_head['LOGFILE'], radpix_15kpc, aper+1)
        
        # plot flux vs radius for the data, uncertainty, model, and residual
        ax = fig.add_subplot(2,2,1)
        ax.scatter(radii, mag(vdat_flux[i]*vfnu[i]/vexp[i]))#/1e5)
        print('vmag: ', mag(vdat_flux[i]*vfnu[i]/vexp[i]))
        print(mag(vdat_flux[i]*vfnu[i]/vexp[i])[aper], np.std(mag(vdat_flux[i]*vfnu[i]/vexp[i])[aper-4:aper+5]))
        print(mag(vdat_flux[i]*vfnu[i]/vexp[i])[default], np.std(mag(vdat_flux[i]*vfnu[i]/vexp[i])[default-4:default+5]))
        #plt.errorbar(radii,vdat_flux[i]/1e5,yerr=vunc_flux[i]/1e5)
        ax.scatter(radii, mag(vunc_flux[i]*vfnu[i]/vexp[i]), marker='s', color='red')
        ax.scatter(radii, mag(vmod_flux[i]*vfnu[i]/vexp[i]), marker='o', color='green')
        ax.scatter(radii, mag(vres_flux[i]*vfnu[i]/vexp[i]), marker='+', color='black')
        ax.set_ylim(27, 18)
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')

        # add labels for axes and a title for the name of the image extension
        #plt.xlabel('Radius [pixels]')#, fontsize=14)
        plt.ylabel('Flux [image units]')#, fontsize=14)
        plt.title(vdat_head['EXTNAME'])

        # repeat for F475W
        ax = fig.add_subplot(2,2,2)
        ax.scatter(radii, mag(udat_flux[i]*ufnu[i]/uexp[i]))
        print('umag: ', mag(udat_flux[i]*ufnu[i]/uexp[i]))
        print(mag(udat_flux[i]*ufnu[i]/uexp[i])[aper], np.std(mag(udat_flux[i]*ufnu[i]/uexp[i])[aper-4:aper+5]))
        print(mag(udat_flux[i]*ufnu[i]/uexp[i])[default], np.std(mag(udat_flux[i]*ufnu[i]/uexp[i])[default-4:default+5]))        
        ax.scatter(radii, mag(uunc_flux[i]*ufnu[i]/uexp[i]), marker='s', color='red')
        ax.scatter(radii, mag(umod_flux[i]*ufnu[i]/uexp[i]), marker='o', color='green')
        ax.scatter(radii, mag(ures_flux[i]*ufnu[i]/uexp[i]), marker='+', color='black')
        #ax.set_ylim(fmin, fmax)
        ax.set_ylim(27, 18)
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')
        
        # add labels for axes and a title for the name of the image extension
        #plt.xlabel('Radius [pixels]')
        #plt.ylabel('Flux [image units]')
        plt.title(udat_head['EXTNAME'])

        # repeat for F160W
        ax = fig.add_subplot(2,2,3)
        ax.scatter(radii, mag(jdat_flux[i]*jfnu[i]/jexp[i]))
        print('jmag: ', mag(jdat_flux[i]*jfnu[i]/jexp[i]))
        print(mag(jdat_flux[i]*jfnu[i]/jexp[i])[aper], np.std(mag(jdat_flux[i]*jfnu[i]/jexp[i])[aper-4:aper+5]))
        print(mag(jdat_flux[i]*jfnu[i]/jexp[i])[default], np.std(mag(jdat_flux[i]*jfnu[i]/jexp[i])[default-4:default+5]))
        ax.scatter(radii, mag(junc_flux[i]*jfnu[i]/jexp[i]), marker='s', color='red')
        ax.scatter(radii, mag(jmod_flux[i]*jfnu[i]/jexp[i]), marker='o', color='green')
        ax.scatter(radii, mag(jres_flux[i]*jfnu[i]/jexp[i]), marker='+', color='black')
        #ax.set_ylim(fmin, fmax*3)
        ax.set_ylim(27, 18)
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')
        
        # add labels for axes and a title for the name of the image extension
        plt.xlabel('Radius [pixels]')
        plt.ylabel('Flux for '+jdat_head['EXTNAME'])
        #plt.title(jdat_head['EXTNAME'])


        # plot F160W / F814W ratio and F814W / F475W ratio for residual flux
        ax = fig.add_subplot(2,2,4)
        #ax.scatter(radii, jres_flux[i] / vres_flux[i], marker='+', color='red')
        ax.scatter(radii, mag(vres_flux[i]*vfnu[i]/vexp[i]) - mag(jres_flux[i]*jfnu[i]/jexp[i]), marker='+', color='red')
        #ax.scatter(radii, vres_flux[i] / ures_flux[i], marker='+', color='green')
        ax.scatter(radii, mag(ures_flux[i]*ufnu[i]/uexp[i]) - mag(vres_flux[i]*vfnu[i]/vexp[i]), marker='+', color='green')
        #ax.set_ylim(-10, 10)
        ax.set_ylim(-1,3)
        plt.xlabel('Radius [pixels]')
        #plt.ylabel('jres / vres')
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')
        # save this page of the pdf file
        pdf.savefig()
        plt.close()

        # start page 2
        fig = plt.figure()

        # have the title at the top of the page be the name of the galfit output 
        plt.suptitle(vdat_head['LOGFILE'])

        # plot flux vs radius for the data, uncertainty, model, and residual
        ax = fig.add_subplot(2,2,1)
        ax.scatter(radii, mag(vdat_annulus[i]*vfnu[i]/vexp[i]))
        ax.scatter(radii, mag(vunc_annulus[i]*vfnu[i]/vexp[i]), marker='s', color='red')
        ax.scatter(radii, mag(vmod_annulus[i]*vfnu[i]/vexp[i]), marker='o', color='green')
        ax.scatter(radii, mag(vres_annulus[i]*vfnu[i]/vexp[i]), marker='+', color='black')
        plt.ylabel('Flux [image units]')#, fontsize=14)
        plt.title(vdat_head['EXTNAME'])
        ax.set_ylim(27,18)
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')

        # repeat for F475W
        ax = fig.add_subplot(2,2,2)
        ax.scatter(radii, mag(udat_annulus[i]*ufnu[i]/uexp[i]))
        ax.scatter(radii, mag(uunc_annulus[i]*ufnu[i]/uexp[i]), marker='s', color='red')
        ax.scatter(radii, mag(umod_annulus[i]*ufnu[i]/uexp[i]), marker='o', color='green')
        ax.scatter(radii, mag(ures_annulus[i]*ufnu[i]/uexp[i]), marker='+', color='black')
        plt.title(udat_head['EXTNAME'])
        ax.set_ylim(27,18)
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')
        
        # repeat for F160W
        ax = fig.add_subplot(2,2,3)
        ax.scatter(radii, mag(jdat_annulus[i]*jfnu[i]/jexp[i]))
        ax.scatter(radii, mag(junc_annulus[i]*jfnu[i]/jexp[i]), marker='s', color='red')
        ax.scatter(radii, mag(jmod_annulus[i]*jfnu[i]/jexp[i]), marker='o', color='green')
        ax.scatter(radii, mag(jres_annulus[i]*jfnu[i]/jexp[i]), marker='+', color='black')
        plt.xlabel('Radius [pixels]')
        plt.ylabel('Flux for '+jres_head['EXTNAME'])
        ax.set_ylim(27,18)
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')
        
        # plot F160W / F814W ratio and F814W / F475W ratio for residual flux
        ax = fig.add_subplot(2,2,4)
        #ax.scatter(radii, jres_annulus[i] / vres_annulus[i], marker='+', color='red')
        #ax.scatter(radii, vres_annulus[i] / ures_annulus[i], marker='+', color='green')
        uv_color = mag(ures_annulus[i]*ufnu[i]/uexp[i]) - mag(vres_annulus[i]*vfnu[i]/vexp[i])
        vj_color = mag(vres_annulus[i]*vfnu[i]/vexp[i]) - mag(jres_annulus[i]*jfnu[i]/jexp[i])
        ax.scatter(radii, vj_color, marker='+', color='red')
        ax.scatter(radii, uv_color, marker='+', color='green')
        #ax.set_ylim(-10, 10)
        ax.set_ylim(-1,3)
        plt.xlabel('Radius [pixels]')
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')
        
        pdf.savefig()
        plt.close()

        # start page 3
        fig = plt.figure()

        ax = fig.add_subplot(2,2,1)
        ax.scatter(radii, vres_annulus[i] / vunc_annulus[i], color='green')
        ax.scatter(radii, ures_annulus[i] / uunc_annulus[i], color='blue')
        ax.scatter(radii, jres_annulus[i] / junc_annulus[i], color='red')
        plt.ylabel('SNR')
        plt.xlabel('Radius [pixels]')
        plt.axvline(x=aper+1, color='black')
        plt.axvline(x=default+1, color='red')

        ax = fig.add_subplot(2,2,4)
        #uv_color = -2.5*np.log10(ures_annulus[i] / vres_annulus[i])
        #vj_color = -2.5*np.log10(vres_annulus[i] / jres_annulus[i])
        ax.scatter(vj_color, uv_color)
        plt.xlabel('[F814W]-[F160W]')
        plt.ylabel('[F475W]-[F814W]')
        ax.set_xlim([-0.5,3.0])
        ax.set_ylim([-0.5,3.0])
        
        pdf.savefig()
        plt.close() 
        
# open the pdf file
os.system('open %s &' % name)

