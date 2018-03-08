# Fahim Khan
# Copied from cb_oii_gvel_flip.py

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
import copy
from copy import deepcopy


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

# specify the plate of interest
from os.path import isdir, join

all_plates = [f for f in os.listdir(mpl5_dir + 'SPX-GAU-MILESHC/') if isdir(join(mpl5_dir + 'SPX-GAU-MILESHC/', f))]
print('List of plates: \n', all_plates, '\n')
plate = int(input('Please enter plate number: '))
lines = glob.glob(mpl5_dir + 'SPX-GAU-MILESHC/*/*/*' + str(plate) + '*MAPS-SPX-GAU-MILESHC.fits*')
filename = 'MLP5_dap_multi_' + str(plate) + '_quicklook.pdf'

# for bookkeeping purposes, here's an array of emission-line names
eml = ['OIId---3728', 'Hb-----4862', 'OIII---4960', 'OIII---5008', 'OI-----6302', 'OI-----6365', 'NII----6549',
       'Ha-----6564', 'NII----6585', 'SII----6718', 'SII----6732']

with PdfPages(filename) as pdf:
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
                    xproj_kpc_map[j, k] = xproj_kpc[j, k] / u.kpc

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

            # Flipping the gvel plot
            x_axis = copy.deepcopy(xproj_kpc_map)
            y_axis = copy.deepcopy(zproj_kpc_map)
            # neg_yax = copy.deepcopy(y_axis)
            # neg_yax = neg_yax * (-1)
            gvel = copy.deepcopy(dap_vel)
            bad = np.abs(gvel) > 300
            gvel[bad] = 0

            pos_y = y_axis > 0
            neg_y = y_axis < 0
            pos_x = x_axis > 0
            neg_x = x_axis < 0
            y_array = y_axis != 0
            x_array = x_axis != 0
            gvel_invert = np.zeros([size, size])

            for i in dap_vel:
                for v in i:
                    if (v != 0) and (-300 < abs(v) < 300):
                        vy, vx = np.where(dap_vel == v)
                                                                       
                        # zval = y_axis[vy[0]][vx[0]]
                        # zprimy, zprimx = np.where(neg_yax == zval)
                        #
                        # gvel_invert[zprimy[0]][zprimx[0]] = v

                        y_val = y_axis[vy[0]][vx[0]]
                        x_val = x_axis[vy[0]][vx[0]]

                        
                        # if y_val > 0 and x_val > 0:
                        #     dist = np.sqrt(np.square(x_axis[x_array]-x_val)+np.square(y_axis[y_array]-y_val))
                        # elif y_val > 0 and x_val < 0:
                        #     dist = np.sqrt(np.square(x_axis[x_array]+x_val)+np.square(y_axis[y_array]-y_val))
                        # elif y_val < 0 and x_val > 0:
                        #     dist = np.sqrt(np.square(x_axis[x_array]-x_val)+np.square(y_axis[y_array]+y_val))
                        # else:
                        #     dist = np.sqrt(np.square(x_axis[x_array]+x_val)+np.square(y_axis[y_array]+y_val))


                        # if y_val > 0 and x_val > 0:
                        #     dist = np.sqrt(np.square(x_axis-x_val)+np.square(y_axis-y_val))
                        # elif y_val > 0 and x_val < 0:
                        #     dist = np.sqrt(np.square(x_axis+x_val)+np.square(y_axis-y_val))
                        # elif y_val < 0 and x_val > 0:
                        #     dist = np.sqrt(np.square(x_axis-x_val)+np.square(y_axis+y_val))
                        # else:
                        #     dist = np.sqrt(np.square(x_axis+x_val)+np.square(y_axis+y_val))

                        dist = np.sqrt(np.square(x_axis-x_val)+np.square(y_axis-y_val))

                        index = np.where(dist == np.amin(dist)) #getting two index values for some y_val, x_val=0
                        flp_y_val = y_axis[y_array][index[0]]
                        #flp_x_val = x_axis[x_array][index[0]]
                        if np.size(flp_y_val) != 1:
                            first = abs(abs(y_val) - abs(flp_y_val[0]))
                            second = abs(abs(y_val) - abs(flp_y_val[1]))
                            if first < second:
                                flp_vy, flp_vx = np.where(y_axis == flp_y_val[0])  #using x_axis == flp_x_val gives the same values
                                gvel_invert[flp_vy[0]][flp_vx[0]] = v
                            else:
                                flp_vy, flp_vx = np.where(y_axis == flp_y_val[1])
                                gvel_invert[flp_vy[0]][flp_vx[0]] = v
                        else:
                            flp_vy, flp_vx = np.where(y_axis == flp_y_val[0])
                            gvel_invert[flp_vy[0]][flp_vx[0]] = v



                        # flipped_vy, flipped_vx = np.where(y_axis == -1*y_val)
                        # gvel_invert[flipped_vy,flipped_vx] = vx

                        #nearest_y = y_axis - y_val
                        #yminy, yminx = np.where(nearest_y == np.amin(nearest_y))
                        #if y_val > 0:
                            #nearest_y = y_axis[neg_y] - y_val
                        #else:
                            #nearest_y = y_axis[pos_y] - y_val
                        # nearest_y.sort(key=lambda x: x[2], reverse=True)
                        # nearest_y.sort(reverse=True)
                        #nearest_y = sorted(nearest_y, reverse=True)

                        # gvel_invert[new_y_coor][new_x_coor] = v

                        #The line method below

                        # if (float(y_val) and float(x_val)) != 0:
                        #
                        #     y_round = float(str(round(y_val, 1)))
                        #     x_round = float(str(round(x_val, 1)))
                        #
                        #     rounded_y = np.around(y_axis, decimals=1)
                        #     rounded_x = np.around(x_axis, decimals=1)
                        #
                        #     # # Two arrays with points that more or less draw a line
                            # o, k = np.where(rounded_y == (-1 * y_round))
                            # n, m = np.where(rounded_x == x_round)
                            #
                            # if (len(o) >1) and (len(n) > 1):
                            #
                            #     def find_line(y_coors, x_coors):
                            #         slope_values = []
                            #         if len(x_coors) > 2:
                            #             for t in range(0, len(x_coors) - 2):
                            #                 curr_slope = (y_coors[t + 1] - y_coors[t]) / (x_coors[t + 1] - x_coors[t])
                            #                 slope_values.append(curr_slope)
                            #             slope = np.mean(slope_values)
                            #         else:
                            #             slope = (y_coors[1] - y_coors[0]) / (x_coors[1] - x_coors[0])
                            #
                            #
                            #         intercept_values = []
                            #         for t in range(0, len(x_coors) - 1):
                            #             curr_intercept = y_coors[t] - (slope * x_coors[t])
                            #             intercept_values.append(curr_intercept)
                            #         intercept = np.mean(intercept_values)
                            #
                            #         return slope, intercept
                            #
                            #
                            #     x_slope, x_intercept = find_line(n, m)
                            #     y_slope, y_intercept = find_line(o, k)
                            #
                            #     point_x_float = (y_intercept-x_intercept)/(x_slope-y_slope)
                            #     point_y_float = (point_x_float * x_slope) + x_intercept
                            #     new_x_coor = int(round(point_x_float))
                            #     new_y_coor = int(round(point_y_float))
                            #
                            #     if (new_x_coor < size) and (new_y_coor<size):
                            #
                            #         gvel_invert[new_y_coor][new_x_coor] = v

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
            plt.suptitle(name)

            # plot 2/3: emission-line flux vs vertical distance
            ax = fig.add_subplot(3, 3, 2)
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(xproj_kpc_map,
                       origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm)
            plt.colorbar()
            plt.title('horizontal distance [kpc]', fontsize=10)

            # plot inverted gvel
            ax = fig.add_subplot(3, 3, 3)
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(gvel_invert, origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm, vmin=-250, vmax=250)
            plt.colorbar()
            plt.title(eml[j] + 'inv-GVEL', fontsize=10)

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

    print(
    "KEEP IN MIND SOMETIMES YOU GET AN ERROR WHEN THE PROGRAM TRIES TO OPEN THE PDF. JUST GO TO THE DIRECTORY WHERE IT WAS CREATED AND OPEN IT MANUALLY.")
    os.system("open %s &" % filename)

# mpl5_dir = os.environ['MANGADIR_MPL5']
# from os.path import isdir, join
# all_plates = [f for f in os.listdir(mpl5_dir+'SPX-GAU-MILESHC/') if isdir(join(mpl5_dir+'SPX-GAU-MILESHC/', f))]
# print('There are ',len(all_plates),' plates in total.')
# start_p = int(input("Start?"))-1
# fin_p = int(input("Finish?"))
# for each_plate in range(start_p, fin_p):
#     main(all_plates[each_plate])

