# Christian Bradna
# Copied from bgl_dap_plates

# Imports


import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import glob
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


def main(plate_rec, min_p, max_p):
    global oii_flux, denomination, denomination, vertical_distance_map
    min_par = float(min_p)
    max_par = float(max_p)

    # array to store interesting galaxies
    interst_main = []

    # function for plotting maps of relevant quantities
    def daplot(quantity, qmin, qmax):
        ax.set_xlim(0, size)
        ax.set_ylim(0, size)
        ax.set_xticklabels(())
        ax.set_yticklabels(())
        plt.imshow(quantity, origin='lower', interpolation='nearest',
                   norm=colors.LogNorm(vmin=qmin, vmax=qmax), cmap=cm.coolwarm)
        plt.colorbar()

    # minimum and maximum emission-line fluxes for plot ranges
    fmin = 1e-19
    fmax = 1e-16

    # define a standard cosmology
    cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

    # directory for relevant Data Reduction Pipeline (drp) and Data Analysis Pipeline (DAP) information
    mpl5_dir = os.environ['MANGADIR_MPL5']  # Be aware that this directory syntax might need revision
    drp = fits.open(mpl5_dir + 'drpall-v2_0_1.fits')
    drpdata = drp[1].data

    # specific the plate of interest
    from os.path import isdir, join
    all_plates = [f for f in os.listdir(mpl5_dir + 'SPX-GAU-MILESHC/') if isdir(join(mpl5_dir + 'SPX-GAU-MILESHC/', f))]
    # print('List of plates: \n',all_plates,'\n')
    plate = plate_rec
    lines = glob.glob(mpl5_dir + 'SPX-GAU-MILESHC/*/*/*' + str(plate) + '*MAPS-SPX-GAU-MILESHC.fits*')
    filename = 'MLP5_dap_multi_' + str(plate) + '_quicklook.pdf'
    current_p = 'temp.pdf'

    # for bookkeeping purposes, here's an array of emission-line names
    eml = ['OIId---3728', 'Hb-----4862', 'OIII---4960', 'OIII---5008', 'OI-----6302', 'OI-----6365', 'NII----6549',
           'Ha-----6564', 'NII----6585', 'SII----6718', 'SII----6732']

    with PdfPages(current_p) as pdf:
        for i in range(0, len(lines)):
            name = lines[i]

            # THIS PART NEEDS ADJUSTING. THERE IS ONE PLATE WITH MORE THAN 4 NUMBERS
            # parse the plate and galaxy information from the filename
            start = name.find('manga-')
            plt1 = start + 6
            plt2 = plt1 + 4
            if str(name[plt2]) != str('-'):
                plt2 = plt2 + 1
            plate = int(name[plt1:plt2])

            end = name.find('-MAPS-SPX-GAU-MILESHC')
            gal2 = end
            gal1 = plt2 + 1
            galaxy = int(name[gal1:gal2])
            print(plate, galaxy)
            denomination = str(plate) + "-" + str(galaxy)

            # find the index of this galaxy in the drpall file
            match = np.where((drpdata.field('plate') == plate) & (drpdata.field('ifudsgn') == str(galaxy)))

            # read in the appropriate DAP file and Halpha flux information
            hdu_dap = fits.open(name)
            # hdu_dap.info()
            dap_ha_sflux = hdu_dap['EMLINE_SFLUX'].data[7, :, :]
            dap_ha_sivar = hdu_dap['EMLINE_SFLUX_IVAR'].data[7, :, :]

            # create arrays that correspond to x and y coordinates for each spaxel
            size = len(hdu_dap['EMLINE_GFLUX'].data[0, :])
            xpos = np.zeros([size, size])
            ypos = np.zeros([size, size])

            for j in range(0, size):
                for k in range(0, size):
                    xpos[j, k] = k
                    ypos[j, k] = j

            # read in the position angle and axis from from the NSA
            theta = drpdata.field('nsa_sersic_phi')[match][0]
            ba = drpdata.field('nsa_sersic_ba')[match][0]
            inc = np.arccos(ba)

            # redefine x=0 and y=0 to be in the center of the field
            xprime = xpos - np.median(xpos)
            yprime = ypos - np.median(ypos)

            # compute the on-sky x and y coordiates defined by the major axis
            trad = np.radians(theta - 90)
            xproj = xprime * np.cos(trad) + yprime * np.sin(trad)
            yproj = xprime * np.sin(trad) * (-1.) + yprime * np.cos(trad)
            zproj = yproj / np.sin(inc) * np.cos(inc)

            # calculate the radius of each pixel in the plane of the disk [units: pixels]
            radpix = np.sqrt(xproj ** 2 + (yproj / ba) ** 2)

            # figure out the conversion between pixel size and kpc
            z = drpdata.field('nsa_z')[match][0]
            radkpc = radpix / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)

            # compute values along the x and y axis of the disk and the z axis above the disk
            xproj_kpc = xproj / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
            yproj_kpc = yproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
            zproj_kpc = zproj / ba / cosmo.arcsec_per_kpc_proper(z) * (0.5 * u.arcsec)
            cen = abs(xproj_kpc) < (1. * u.kpc)

            radkpc_map = np.zeros([size, size])
            xproj_kpc_map = np.zeros([size, size])
            yproj_kpc_map = np.zeros([size, size])
            zproj_kpc_map = np.zeros([size, size])

            for j in range(0, size):
                for k in range(0, size):
                    if (dap_ha_sivar[j, k] > 0.):
                        radkpc_map[j, k] = radkpc[j, k] / u.kpc
                        yproj_kpc_map[j, k] = yproj_kpc[j, k] / u.kpc
                        zproj_kpc_map[j, k] = zproj_kpc[j, k] / u.kpc
                    else:
                        radkpc_map[j, k] = radkpc[j, k] / u.kpc * (-1.)

            # hdu_dap['EMLINE_GFLUX'].header
            # C01     = 'OIId---3728'        / Species in cube channel 1
            # C02     = 'Hb-----4862'        / Species in cube channel 2
            # C03     = 'OIII---4960'        / Species in cube channel 3
            # C04     = 'OIII---5008'        / Species in cube channel 4
            # C05     = 'OI-----6302'        / Species in cube channel 5
            # C06     = 'OI-----6365'        / Species in cube channel 6
            # C07     = 'NII----6549'        / Species in cube channel 7
            # C08     = 'Ha-----6564'        / Species in cube channel 8
            # C09     = 'NII----6585'        / Species in cube channel 9
            # C10     = 'SII----6718'        / Species in cube channel 10
            # C11     = 'SII----6732'        / Species in cube channel 11

            # hdu_dap['EMLINE_GFLUX'].data.shape # (11, 74, 74)

            # for j in range(0, len(hdu_dap['EMLINE_GVEL'].data)):
            # focus on the [O II] emisison line
            for j in range(0, 1):
                print(eml[j])

                # DAP: emission-line quantities
                dap_vel = hdu_dap['EMLINE_GVEL'].data[j, :, :]
                dap_vel_ivar = hdu_dap['EMLINE_GVEL_IVAR'].data[j, :, :]
                dap_vel_err = np.sqrt(1. / dap_vel_ivar)
                dap_sig = hdu_dap['EMLINE_GSIGMA'].data[j, :, :]
                dap_sig_ivar = hdu_dap['EMLINE_GSIGMA'].data[j, :, :]
                dap_sig_err = np.sqrt(1. / dap_sig_ivar)
                dap_flux = hdu_dap['EMLINE_GFLUX'].data[j, :, :]
                dap_mask = hdu_dap['EMLINE_GFLUX_MASK'].data[j, :, :]
                dap_ivar = hdu_dap['EMLINE_GFLUX_IVAR'].data[j, :, :]
                dap_err = np.sqrt(1. / dap_ivar)
                dap_sflux = hdu_dap['EMLINE_SFLUX'].data[j, :, :]
                dap_smask = hdu_dap['EMLINE_SFLUX_MASK'].data[j, :, :]
                dap_sivar = hdu_dap['EMLINE_SFLUX_IVAR'].data[j, :, :]
                dap_serr = np.sqrt(1. / dap_sivar)  # *1.e-4

                # Galaxy rating

                # Masking out ratio plot. Removing x < 1 (and nans).

                # for subarr in oii_to_ha:
                #     for element in subarr:
                #         if element < 0:
                #             element = 0
                oii_to_ha = (dap_sflux / dap_ha_sflux)
                bad = (oii_to_ha < 0.7)
                oii_to_ha[bad] = 0
                nans = np.isnan(oii_to_ha)
                oii_to_ha[nans] = 0
                infs = np.isinf(oii_to_ha)
                oii_to_ha[infs] = 0
                unreas = oii_to_ha > ((np.mean(oii_to_ha))**3)
                oii_to_ha[unreas] = 0

                # Filtering values nearest to the galaxy disk according to a filtering parameter

                filt_parm = ((90 / (inc * 180 / np.pi)) ** 2) * (np.amax(zproj_kpc_map) / ((inc * 180 / np.pi) ** 0.4))
                vertical_distance_map = copy.deepcopy(zproj_kpc_map)

                # for subarr in vertical_distance_map:
                #     vert_max = 0
                #     for element in subarr:
                #         if vert_max < element:
                #             vert_max = element

                # for subarr in vertical_distance_map:
                #     for element in subarr:
                #         if np.abs(element) < filt_parm:
                #             element = 0

                bad = (np.abs(vertical_distance_map) < filt_parm)
                vertical_distance_map[bad] = 0

                # Filtering values of flux < 0.1

                oii_flux = copy.deepcopy(dap_sflux)
                # for subarr in oii_flux:
                #     for element in subarr:
                #         if element < 0:
                #             element = 0

                bad = (oii_flux < 0.1)
                oii_flux[bad] = 0

                # bradna_index = np.abs( (((vertical_distance_map * ((inc*180/np.pi) ** 6)) ** 4) * oii_flux) * ((vertical_distance_map ** 2) * oii_to_ha) *
                #                        (((vertical_distance_map * ((inc*180/np.pi) ** 6)) ** 2) * (oii_flux / dap_serr)) * oii_flux / 1e17)  # removed * ((inc*180/np.pi) ** 5)


                # Actual formula for a rating
                # bradna_index = np.abs( (((vertical_distance_map * ((inc*180/np.pi) ** 6)) ** 4) * ((10*oii_flux)**2)) * ((vertical_distance_map ** 2) * oii_to_ha) *
                #                        (((vertical_distance_map * ((inc*180/np.pi) ** 6)) ** 2) * (oii_flux / dap_serr)) * (oii_flux) / 1e20)  # removed * ((inc*180/np.pi) ** 5)
                normalize_vdist = (inc * 180 / np.pi) / (2 * np.amax(zproj_kpc_map))
                bradna_index = np.abs((10 * oii_flux ** 2) * (oii_to_ha ** 2) * (oii_flux / dap_serr) * (vertical_distance_map ** .5)) # * (normalize_vdist ** 2))

                # Calculating and printing stuff
                param = (bradna_index > 0)
                bradna_index_sum = np.sum(bradna_index[param])
                if bradna_index_sum == 0:
                    bradna_index_sum = 10

                if (min_par < (inc * 180 / np.pi) < max_par) & (np.sum(oii_flux) > 540):
                    state = '>>>>>>>>>>>>>>>>Interesting'
                else:
                    state = 'Not Interesting'
                b_index = math.log(bradna_index_sum, 10)
                print(plate, galaxy, 'B_index: ', b_index, state, 'Rating: ', b_index*(100/57)) # math.log(bradna_index_sum, 10)
                print("OIId flux", np.sum(oii_flux), '\n', 'Flux sum to rating ratio: ', (np.sum(oii_flux)/b_index), 'Sum', bradna_index_sum, '\n' )
                denom_and_rating = str(denomination) + ' Rating: ' + str(b_index*(100/57)) + ' Flux/rating ratio: ' + str((np.sum(oii_flux)/b_index))


                # Plotting interesting galaxies
                if state == '>>>>>>>>>>>>>>>>Interesting':
                    # plot 1: galaxy coordinates in kpc
                    fig = plt.figure()
                    ax = fig.add_subplot(3, 3, 1)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(zproj_kpc_map,
                               origin='lower',
                               interpolation='nearest',
                               cmap=cm.coolwarm)
                    plt.colorbar()
                    plt.title('vertical distance [kpc]', fontsize=10)
                    plt.suptitle(denom_and_rating, fontsize = 10)

                    # plot 2/3: emission-line flux vs vertical distance
                    ax = fig.add_subplot(4, 2, 2)
                    ax.set_yscale("log", nonposy='clip')
                    ax.set_ylim(fmin, fmax)
                    ax.set_xlim(np.max(zproj_kpc_map) * (-1.), np.max(zproj_kpc_map))
                    ax.scatter(zproj_kpc[cen], dap_ha_sflux[cen] * 1.e-17, lw=0, s=1)
                    ax.scatter(zproj_kpc, dap_ha_sflux * 1.e-17, lw=0, s=0.2)
                    ax.scatter(zproj_kpc[cen], dap_sflux[cen] * 1.e-17, lw=0, color='red', s=1)
                    ax.scatter(zproj_kpc, dap_sflux * 1.e-17, lw=0, color='red', s=0.2)
                    plt.title('emission-line flux vs vertical distance', fontsize=9)
                    plt.text(0., fmin * 2., 'i=' + str(round(inc * 180. / np.pi, 2)), fontsize=9)

                    # plot 4: emission-line SFLUX
                    ax = fig.add_subplot(3, 3, 4)
                    daplot(dap_sflux * 1.e-17, fmin, fmax)
                    plt.title(eml[j] + ' SFLUX', fontsize=10)

                    # plot 5: emission-line SFLUX serror
                    ax = fig.add_subplot(3, 3, 5)
                    daplot(dap_serr * 1.e-17, fmin, fmax)
                    plt.title(eml[j] + ' SFLUX Error', fontsize=10)

                    # plot 6: emission-line SFLUX signal-to-noise ratio
                    ax = fig.add_subplot(3, 3, 6)
                    daplot(dap_sflux / dap_serr, 0.1, 10.)
                    plt.title(eml[j] + ' SFLUX S/N', fontsize=10)

                    # plot 7: emission-line flux / Halpha flux
                    ax = fig.add_subplot(3, 3, 7)
                    daplot(dap_sflux / dap_ha_sflux, 0.1, 10.)
                    plt.title(eml[j] + '/HA', fontsize=10)

                    # plot 8: emission-line velocity
                    ax = fig.add_subplot(3, 3, 8)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(dap_vel, origin='lower',
                               interpolation='nearest',
                               cmap=cm.coolwarm, vmin=-250, vmax=250)
                    plt.colorbar()
                    plt.title(eml[j] + ' GVEL', fontsize=10)

                    # plot 9: emission-line dispersion
                    ax = fig.add_subplot(3, 3, 9)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(dap_sig, origin='lower',
                               interpolation='nearest',
                               cmap=cm.YlOrRd, vmin=0, vmax=250)
                    plt.colorbar()
                    plt.title(eml[j] + ' GSIGMA', fontsize=10)

                    pdf.savefig()
                    plt.close()

                    #Plotting the new arrays
                    fig = plt.figure()
                    # plot 1: galaxy coordinates in kpc but without filtered values
                    fig = plt.figure()
                    ax = fig.add_subplot(2, 2, 1)
                    ax.set_xlim(0, size)
                    ax.set_ylim(0, size)
                    ax.set_xticklabels(())
                    ax.set_yticklabels(())
                    plt.imshow(vertical_distance_map,
                               origin='lower',
                               interpolation='nearest',
                               cmap=cm.coolwarm)
                    plt.colorbar()
                    plt.title('filtered vertical distance', fontsize=10)
                    plt.suptitle(denom_and_rating, fontsize = 10)

                    # plot 4: flux
                    ax = fig.add_subplot(2, 2, 2)
                    daplot(oii_flux * np.abs(vertical_distance_map) / np.abs(vertical_distance_map) * 1.e-17, fmin, fmax)
                    plt.title(eml[j] + ' SFLUX', fontsize=10)

                    # plot 7: oii/ha
                    ax = fig.add_subplot(2, 2, 3)
                    daplot(oii_to_ha * np.abs(vertical_distance_map) / np.abs(vertical_distance_map), 0.1, 10.)
                    plt.title(eml[j] + '/HA', fontsize=10)


                    pdf.savefig()
                    plt.close()
                    plt.close(fig)




    return (current_p)


mpl5_dir = os.environ['MANGADIR_MPL5']
from PyPDF2 import PdfFileMerger
import copy
from copy import deepcopy

interst = []
from os.path import isdir, join
import math
all_plates = [f for f in os.listdir(mpl5_dir + 'SPX-GAU-MILESHC/') if isdir(join(mpl5_dir + 'SPX-GAU-MILESHC/', f))]
print('There are ', len(all_plates), ' plates in total.')
start_p = int(input("Start?")) - 1
fin_p = int(input("Finish?"))
min_p = input('Min inclination?')
max_p = input('Max inclination?')
# main_pdf = 'cool_galaxies_main.pdf'
merger = PdfFileMerger()
for each_plate in range(start_p, fin_p):
    temp_pdf = main(all_plates[each_plate], min_p, max_p)
    merger.append(open(temp_pdf, 'rb'))

with open('result.pdf', 'wb') as fout:
    merger.write(fout)

print(
"KEEP IN MIND SOMETIMES YOU GET AN ERROR WHEN THE PROGRAM TRIES TO OPEN THE PDF. JUST GO TO THE DIRECTORY WHERE IT WAS CREATED AND OPEN IT MANUALLY.")
os.system("open %s &" % 'result.pdf')
