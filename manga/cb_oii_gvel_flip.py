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

# specific the plate of interest
from os.path import isdir, join

all_plates = [f for f in os.listdir(mpl5_dir + 'SPX-GAU-MILESHC/') if isdir(join(mpl5_dir + 'SPX-GAU-MILESHC/', f))]
print('List of plates: \n', all_plates, '\n')
plate = int(input('Please enter plate number: '))
lines = glob.glob(mpl5_dir + 'SPX-GAU-MILESHC/*/*/*' + str(plate) + '*MAPS-SPX-GAU-MILESHC.fits*')
filename = 'MLP5_dap_multi_' + str(plate) + '_gvel_flip.pdf'

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




            # CB new approach: making the vel plot circular so we have better symmetry to work with
            x_axis = copy.deepcopy(xproj_kpc_map)
            y_axis = copy.deepcopy(zproj_kpc_map)

            pos_y = y_axis > 0
            neg_y = y_axis < 0
            pos_x = x_axis > 0
            neg_x = x_axis < 0
            gvel_invert = np.zeros([size, size])
            backg_val = abs(dap_vel[0][0])
            gvel_invert = gvel_invert - backg_val


            newx = copy.deepcopy(x_axis)
            newx = newx[np.logical_and(newx != 0, newx !=0)]
            newxt = copy.deepcopy(newx)
            newxt = np.transpose(newxt)
            newxt = newxt[np.logical_and(newxt != 0, newxt != 0)]

            yind, xdum = np.where(x_axis == newx[0])
            ydum, xind = np.where(x_axis == newxt[0])
            yinval = yind[0]
            xinval = xind[0]

            #figuring out the middle of the plot
            midy = int(int(len(dap_vel))/2)
            radius = midy-yinval
            mask_1 = np.zeros([size, size])
            mask_1 = mask_1 + 1

            # creating circle; masking out values that are further away than the radius previously defined
            for i in range(0,len(dap_vel)):
                for valz in range(0, len(dap_vel[i])):
                    dist = ((midy-i)**2 + (midy-valz)**2)**0.5
                    if dist > radius:
                        mask_1[i][valz] = 0



            # Flipping the gvel plot
            x_axis = x_axis * mask_1
            y_axis = y_axis * mask_1
            masked_vel = copy.deepcopy(dap_vel)
            masked_vel = masked_vel * mask_1
            # bad = np.abs(masked_vel) > 1000
            # masked_vel[bad] = 0




            icount = 0
            for i in masked_vel:
                vcount = 0
                for v in i:
                    vy, vx = np.where(masked_vel == v) #Positions in vel plot

                    if v != 0:

                        y_val = y_axis[vy[0]][vx[0]] # Value in z-proj
                        x_val = x_axis[vy[0]][vx[0]] # Value in x-proj

                        # nearest_y = y_axis - y_val
                        # yminy, yminx = np.where(nearest_y == np.amin(np.abs(nearest_y)))
                        if y_val != 0:
                            if y_val > 0: #what if y_val is zero?? cb: solved
                                nearest_y = y_axis[(y_axis < 0)] + y_val #difference in opposite y-axis
                            if y_val < 0:
                                nearest_y = y_axis[(y_axis > 0)] + y_val


                            nearest_y = sorted(nearest_y, key=abs) #sorted abs values of difference from smallest to largest
                            ref_y = nearest_y - y_val #sorted corresponding y_val in the opposite axis
                            fil_ref_x = np.array([])
                            fil_ref_y = np.array([])

                            for vals in ref_y: #getting rid of unwanted ref_y
                                if np.abs((np.abs(vals) - np.abs(y_val))) < 3: #what if y_val is greater than ref_y?? cb:solved maybe
                                    fil_ref_y = np.append(fil_ref_y, vals)

                            fil_ref_y2 = np.array([])

                            count_f = 0
                            count_nf = 0
                            for k in fil_ref_y:

                                ####fix 1
                                kind_neary = np.where(nearest_y == (k+y_val))
                                if y_val > 0 :
                                    kind_yax = np.where((y_axis[(y_axis < 0)]+y_val) == (nearest_y[kind_neary[0][0]]))
                                    refxy, refxx = np.where(y_axis == (y_axis[(y_axis < 0)][kind_yax[0]]))
                                if y_val < 0:
                                    kind_yax = np.where((y_axis[(y_axis > 0)]+y_val) == (nearest_y[kind_neary[0][0]]))
                                    refxy, refxx = np.where(y_axis == (y_axis[(y_axis > 0)][kind_yax[0]]))

                                ####fix 1

                                # version for debugging
                                # if len(refxy) == 0:
                                #     # print("WARNING: Original position of ",k," was not found in z_proj. One pixel lost. Length of refxy ",len(refxy))
                                #     b = 10**100
                                #     fil_ref_x = np.append(fil_ref_x, b)
                                #     count_nf = count_nf + 1
                                #     # print("count of missing pixels ",count_nf)
                                # if len(refxy) != 0:


                                refx_val = x_axis[refxy[0]][refxx[0]]
                                fil_ref_x = np.append(fil_ref_x, refx_val) #corresponding x_val of remaining ref_y
                                fil_ref_y2 = np.append(fil_ref_y2, y_axis[refxy[0]][refxx[0]])
                                count_f = count_f + 1

                                # WORKING version

                                # if len(refxy) == 0:
                                #     # print("WARNING: Original position of ",k," was not found in z_proj. One pixel lost. Length of refxy ",len(refxy))
                                #     b = 10 ** 100
                                #     fil_ref_x = np.append(fil_ref_x, b)
                                #     count_nf = count_nf + 1
                                #     # print("count of missing pixels ",count_nf)
                                # if len(refxy) != 0:
                                #     refx_val = x_axis[refxy[0]][refxx[0]]
                                #     fil_ref_x = np.append(fil_ref_x,
                                #                           refx_val)  # corresponding x_val of remaining ref_y
                                #     count_f = count_f + 1


                                        # print(k, " Appended. lenrefxy ",len(refxy)," Count of appended pixels ",count_f)
                                #Why is the size of fil_ref_x smaller than the size of fill_ref_y??

                            fil_y = copy.deepcopy(fil_ref_y2)
                            fil_x = copy.deepcopy(fil_ref_x)

                            # dist = (np.square(np.abs(fil_y)-abs(y_val)) + np.square(np.abs(fil_x)-abs(x_val)))**0.5 double flipping it if abs
                            dist = (np.square(fil_y+y_val) + np.square(fil_x-x_val))**0.5 #might be worth fiddling around this line. Somehow this line has control over pixel loss

                            if len(dist) == 0:
                                print('Pixel lost. -- ',plate, galaxy, "  Indexes: ",icount,' ',vcount)


                            else:
                                min_dist_ind = np.where(dist == np.amin(dist))
                                new_y_coor, new_x_coor = np.where(y_axis == fil_ref_y2[min_dist_ind[0][0]])
                                if (len(new_y_coor) == 0) or (len(new_x_coor) == 0):
                                    print("Warning. No pixel is being appended for original vel val: ",k,"\n Info: Elements in fil_y, fil_x and dist: ",
                                          len(fil_y)," ; ",len(fil_x)," ; ",len(dist),"\n For plot ",plate, galaxy)



                            gvel_invert[new_y_coor[0]][new_x_coor[0]] = v

                        if y_val == 0:
                            if gvel_invert[vy[0]][vx[0]] != (backg_val*(-1)):
                                print("A pixel is being overwritten and values are potentially being lost.",
                                      "\n Info: Elements in fil_y, fil_x and dist: ",
                                      len(fil_y), " ; ", len(fil_x), " ; ", len(dist), "\n For plot ", plate, galaxy
                                      )
                            gvel_invert[vy[0]][vx[0]] = v
                    else:
                        gvel_invert[icount][vcount] = v

                    vcount = vcount + 1
                icount = icount + 1

            gvinv = copy.deepcopy(gvel_invert)
            gvinv = sorted(np.ndarray.flatten(gvinv))
            gvmask = sorted(np.ndarray.flatten(masked_vel))


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
            assymetry_p = dap_vel - gvel_invert

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
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(assymetry_p, origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm, vmin=-250, vmax=250)
            plt.colorbar()
            plt.title(eml[j] + 'subtracted', fontsize=10)

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
            plt.imshow(masked_vel, origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm, vmin=-250, vmax=250)
            plt.colorbar()
            plt.title(eml[j] + 'masked-GVEL', fontsize=10)

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
