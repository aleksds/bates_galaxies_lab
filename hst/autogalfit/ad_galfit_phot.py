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
from photutils import CircularAperture
from photutils import aperture_photometry

# the code requires the user to specify an input directory that has output fits files from GALFIT
dir = sys.argv[1]
gal = sys.argv[2]

# find all of the relevant fits files
files = glob.glob(dir+'/'+gal+'*output.fits')

# define relevant variables
xcen = 50
ycen = 50
fmin = -2
fmax = 12
#radii = np.arange(40)+1
#radii = 10.**(np.arange(21)/20*1.9-0.3)
#radii = np.array([1.,2.,3.,4.,5.,10.,15.,20.,25.,30.,35.,40.,45.])
radii = np.array([5.,10.,20.,30.,40.])

# arrays for v band
vdat_flux = np.zeros([len(files),len(radii)])
vunc_flux = np.zeros([len(files),len(radii)])
vmod_flux = np.zeros([len(files),len(radii)])
vres_flux = np.zeros([len(files),len(radii)])
vres_sub = np.zeros([len(files),len(radii)])

# arrays for u band
udat_flux = np.zeros([len(files),len(radii)])
uunc_flux = np.zeros([len(files),len(radii)])
umod_flux = np.zeros([len(files),len(radii)])
ures_flux = np.zeros([len(files),len(radii)])
ures_sub = np.zeros([len(files),len(radii)])

# arrays for j band
jdat_flux = np.zeros([len(files),len(radii)])
junc_flux = np.zeros([len(files),len(radii)])
jmod_flux = np.zeros([len(files),len(radii)])
jres_flux = np.zeros([len(files),len(radii)])
jres_sub = np.zeros([len(files),len(radii)])

# let's produce plots that will be saved in a pdf file
name = 'galfit_phot_'+dir+'_'+gal+'.pdf'

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
        vdat, vdat_head = hdu[0].data, hdu[0].header
        vunc, vunc_head = hdu[12].data, hdu[12].header
        vmod, vmod_head = hdu[3].data, hdu[3].header
        vres, vres_head = hdu[6].data, hdu[6].header

        # F475W corresponds roughly to rest-frame U-band
        udat, udat_head = hdu[1].data, hdu[1].header
        uunc, uunc_head = hdu[13].data, hdu[13].header
        umod, umod_head = hdu[4].data, hdu[4].header
        ures, ures_head = hdu[7].data, hdu[7].header        

        # F160W corresponds roughly to rest-frame J-band
        jdat, jdat_head = hdu[2].data, hdu[2].header
        junc, junc_head = hdu[14].data, hdu[14].header
        jmod, jmod_head = hdu[5].data, hdu[5].header
        jres, jres_head = hdu[8].data, hdu[8].header  
        
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

            # define annular values
            if j == 0:
                vres_sub[i,j] = vres_flux[i,j]
                ures_sub[i,j] = ures_flux[i,j]
                jres_sub[i,j] = jres_flux[i,j]
            else:
                vres_sub[i,j] = vres_flux[i,j]-vres_flux[i,j-1]
                ures_sub[i,j] = ures_flux[i,j]-ures_flux[i,j-1]
                jres_sub[i,j] = jres_flux[i,j]-jres_flux[i,j-1]
                          
        # start our figure
        fig = plt.figure()

        # have the title at the top of the page be the name of the galfit output 
        plt.suptitle(vdat_head['LOGFILE'])

        # plot flux vs radius for the data, uncertainty, model, and residual
        ax = fig.add_subplot(2,2,1)
        ax.scatter(radii, vdat_flux[i]/1e5)
        #plt.errorbar(radii,vdat_flux[i]/1e5,yerr=vunc_flux[i]/1e5)
        ax.scatter(radii, vunc_flux[i]/1e5, marker='s', color='red')
        ax.scatter(radii, vmod_flux[i]/1e5, marker='o', color='green')
        ax.scatter(radii, vres_flux[i]/1e5, marker='+', color='black')
        ax.set_ylim(fmin, fmax)

        # add labels for axes and a title for the name of the image extension
        #plt.xlabel('Radius [pixels]')#, fontsize=14)
        plt.ylabel('Flux [image units]')#, fontsize=14)
        plt.title(vdat_head['EXTNAME'])

        # repeat for F475W
        ax = fig.add_subplot(2,2,2)
        ax.scatter(radii, udat_flux[i]/1e5)
        ax.scatter(radii, uunc_flux[i]/1e5, marker='s', color='red')
        ax.scatter(radii, umod_flux[i]/1e5, marker='o', color='green')
        ax.scatter(radii, ures_flux[i]/1e5, marker='+', color='black')
        ax.set_ylim(fmin, fmax)
        
        # add labels for axes and a title for the name of the image extension
        #plt.xlabel('Radius [pixels]')
        #plt.ylabel('Flux [image units]')
        plt.title(udat_head['EXTNAME'])

        # repeat for F160W
        ax = fig.add_subplot(2,2,3)
        ax.scatter(radii, jdat_flux[i]/1e5)
        ax.scatter(radii, junc_flux[i]/1e5, marker='s', color='red')
        ax.scatter(radii, jmod_flux[i]/1e5, marker='o', color='green')
        ax.scatter(radii, jres_flux[i]/1e5, marker='+', color='black')
        ax.set_ylim(-5, 35)
        
        # add labels for axes and a title for the name of the image extension
        plt.xlabel('Radius [pixels]')
        plt.ylabel('Flux for '+jdat_head['EXTNAME'])
        #plt.title(jdat_head['EXTNAME'])


        # plot F160W / F814W ratio and F814W / F475W ratio for residual flux
        ax = fig.add_subplot(2,2,4)
        ax.scatter(radii, jres_flux[i] / vres_flux[i], marker='D', color='red', label='F160W/F814W integrated')
        ax.scatter(radii, jres_sub[i] / vres_sub[i], marker='*', color='red', label='F160W/F814W annulus')
        ax.scatter(radii, vres_flux[i] / ures_flux[i], marker='D', color='green', label='F814W/F475W integrated')
        ax.scatter(radii, vres_sub[i] / ures_sub[i], marker='*', color='green', label='F814W/F475W annulus')
        ax.set_ylim(-5, 25)
        plt.xlabel('Radius [pixels]')
        #plt.ylabel('jres / vres')
        plt.legend(loc='upper left', fontsize=7)
        
        # save this page of the pdf file
        pdf.savefig()
        plt.close()


    # additional page that explores colors vs radius as a function of model 
    colors=['red','orange','yellow','green','blue','indigo','violet'] 
    fig = plt.figure()

    ax = fig.add_subplot(2,2,1)
    ax.set_ylim(-5, 25)
    plt.ylabel('F160W/F814W integrated')
    for i in range(0, len(files)):
        ax.scatter(radii, jres_flux[i] / vres_flux[i], marker='D', color=colors[i])

    ax = fig.add_subplot(2,2,2)
    plt.title('F160W/F814W annulus')
    for i in range(0, len(files)):
        ax.scatter(radii, jres_sub[i] / vres_sub[i], marker='*', color=colors[i])
    ax.set_ylim(-5, 25)

    ax = fig.add_subplot(2,2,3)
    plt.ylabel('F814W/F475W integrated')
    for i in range(0, len(files)):
        ax.scatter(radii, vres_flux[i] / ures_flux[i], marker='D', color=colors[i])
    ax.set_ylim(-5, 25)

    ax = fig.add_subplot(2,2,4)
    plt.xlabel('F814W/F475W annulus')
    for i in range(0, len(files)):
        ax.scatter(radii, vres_sub[i] / ures_sub[i], marker='*', color=colors[i])
    ax.set_ylim(-5, 25)

    pdf.savefig()
    plt.close()
    
    # additional page that looks at colors in magnitude space
    colors=['red','orange','yellow','green','blue','indigo','violet'] 
    fig = plt.figure()

    ax = fig.add_subplot(2,2,1)
    ax.set_ylim(-0.5,3.5)
    ax.set_xlim(-0.5,3.5)
    plt.title('Integrated colors')
    plt.ylabel('[F475W]-[F814W]')
    plt.xlabel('[F814W]-[F160W]')
    for i in range(0, len(files)):
        ax.scatter(2.5*np.log10(jres_flux[i][2:len(radii)] / vres_flux[i][2:len(radii)]), 2.5*np.log10(vres_flux[i][2:len(radii)] / ures_flux[i][2:len(radii)]), marker='D', color=colors[i])


    ax = fig.add_subplot(2,2,4)
    ax.set_ylim(-0.5,3.5)
    ax.set_xlim(-0.5,3.5)
    plt.title('Annular colors')
    plt.ylabel('[F475W]-[F814W]')
    plt.xlabel('[F814W]-[F160W]')
    for i in range(0, len(files)):
        ax.scatter(2.5*np.log10(jres_sub[i][2:len(radii)] / vres_sub[i][2:len(radii)]), 2.5*np.log10(vres_sub[i][2:len(radii)] / ures_sub[i][2:len(radii)]), marker='*', color=colors[i])
        
    pdf.savefig()
    plt.close

# open the pdf file
os.system('open %s &' % name)

