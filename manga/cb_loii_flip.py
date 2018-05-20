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
import scipy
from scipy import ndimage


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
            x_axis_2 = np.zeros([size,size])
            x_axis_2 = x_axis_2 - 65000
            for i in range(0,size):
                for k in range(0,size):
                    if (dap_ha_sivar[i, k] > 0.):
                        x_axis_2[i,k] = xproj_kpc[i, k] / u.kpc
            y_axis_2 = np.zeros([size,size])
            y_axis_2 = y_axis_2 - 35
            for i in range(0,size):
                for k in range(0,size):
                    if (dap_ha_sivar[i, k] > 0.):
                        y_axis_2[i,k] = yproj_kpc[i, k] / u.kpc

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



            # Applying masks
            # x_axis = x_axis * mask_1
            # y_axis = y_axis * mask_1
            masked_vel = copy.deepcopy(dap_vel)
            masked_vel = masked_vel * mask_1
            bad = np.abs(masked_vel) > 1000
            masked_vel[bad] = 0


            def isnumber(s):
                try:
                    float(s)
                except ValueError:
                    return False
                return True


            def find_line(y_coors, x_coors):
                slope_values = []
                if len(x_coors) > 2:
                    for t in range(0, len(x_coors) - 2):
                        curr_slope = (y_coors[t + 1] - y_coors[t]) / (x_coors[t + 1] - x_coors[t])
                        if isnumber(curr_slope):
                            slope_values.append(curr_slope)
                    slope = np.mean(slope_values)
                else:
                    slope = (y_coors[1] - y_coors[0]) / (x_coors[1] - x_coors[0])
                    print('Single value slope for ',galaxy,'-',plate)

                intercept_values = []
                for t in range(0, len(x_coors) - 1):
                    curr_intercept = y_coors[t] - (slope * x_coors[t])
                    intercept_values.append(curr_intercept)
                intercept = np.mean(intercept_values)

                return slope, intercept




            flipped_y = np.flip(y_axis_2,1)


            # Drawing a line at x = 0
            xi, xo = np.where(np.around(y_axis_2, decimals=0) == 0.0) #round x axis, find the zeros (should be a horizontal line through
            #the major axis
            xi_f,xo_f = np.where(np.around(flipped_y, decimals=0) == 0.0)


            mainx_slope, mainx_intr = find_line(xi,xo) #this is the slope of the line passing through the
            # 0 points in the xproj map, and the y-intercept of such line. This describes the line equation.
            slope_f,int_f = find_line(xi_f,xo_f)
            if np.isinf(mainx_slope) or np.isinf(slope_f):
                print('Infinite slope for ',galaxy,'-',plate,' Main slope isinf ',np.isinf(mainx_slope),' flip slope isinf ',np.isinf(slope_f))


            if mainx_slope>0:
                arbit_val = 5
            if mainx_slope<0:
                arbit_val = -5

            x_o = arbit_val/mainx_slope
            x_f = arbit_val/slope_f
            theta_o = np.abs(np.arctan(arbit_val/x_o))
            theta_f = np.abs(np.arctan(arbit_val/x_f))

            phi = np.pi-(theta_o+theta_f)
            phi = np.degrees(phi) *-1


            #calculating a second angle






            def method_2():

                gvel_invert = copy.deepcopy(dap_vel)
                gvel_invert = np.flip(gvel_invert,1)

                center = round(size/2)
                arbit_vel = dap_vel[center+5,center+5]
                vy, vx = np.where(np.flip(dap_vel,1) == dap_vel[center+5,center+5])
                x_fl = abs(vx[0]-center)
                y_fl = abs(vy[0]-center)
                theta_or = np.abs(np.arctan(1))
                theta_fl = np.abs(np.arctan(y_fl/x_fl))
                if ((vx[0]-center)>0) and (abs(vy[0]-center)>0):
                    phi_2 = theta_or - theta_fl
                if ((vx[0]-center)<0) and (abs(vy[0]-center)>0):
                    phi_2 = -1*(np.pi-(theta_or+theta_fl))
                if ((vx[0]-center)<0) and (abs(vy[0]-center)<0):
                    phi_2 = -1*((np.pi/2)-theta_or+(np.pi/2)+theta_fl)
                if ((vx[0]-center)>0) and (abs(vy[0]-center)<0):
                    phi_2 = -1*(theta_or+theta_fl)

                phi_2 = np.degrees(phi_2)



            gvel_invert = scipy.ndimage.rotate(gvel_invert,phi_2,reshape=False)




            # A = np.asarray([[0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            #      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]])








            # icount = 0
            # for i in dap_vel:
            #     vcount = 0
            #     for v in i:
            #         # vy, vx = np.where(dap_vel == v) #Positions in vel plot
            #         vy = [icount]
            #         vx = [vcount]
            #
            #
            #         y_val = y_axis[vy[0]][vx[0]]  # Value in z-proj
            #         x_val = x_axis_2[vy[0]][vx[0]] # Value in x-proj
            #
            #         # Run only if value is not -65000 for x prok
            #         if x_val != (-65000):
            #             print(plate, '-', galaxy, ' Analizing ', v)
            #             potential_zeros_in_line = []
            #             particular_intercept = vy[0] - (mainx_slope*vx[0])
            #             main_ax_yo, main_ax_yi = np.where(np.around(y_axis_2, decimals=0) == 0.0)
            #             for zero1 in range(0,len(main_ax_yo)):
            #                 y_eq = round(main_ax_yo[zero1])
            #                 x_eq = round(mainx_slope*main_ax_yi[zero1] + particular_intercept)
            #                 if (abs(y_eq-x_eq))<3:
            #                     potential_zeros_in_line.append([main_ax_yo[zero1],main_ax_yi[zero1]])
            #                         # up to this point, we have candidates for the spaxels we can consider
            #                         # as zeros for our line.
            #             # for now, let's look at only the first candidate in potential_zeros
            #             yind_diff = potential_zeros_in_line[0][0] - vy[0]
            #             xind_diff = potential_zeros_in_line[0][1] - vx[0]
            #             new_yind_cord = vy[0] + (2*yind_diff)
            #             new_xind_cord = vx[0] + (2*xind_diff)
            #
            #
            #             gvel_invert[new_yind_cord][new_xind_cord] = v
            #             print(plate,'-',galaxy,' Appended ',v,' @ ',new_yind_cord,', ',new_xind_cord)
            #         else:
            #             if x_val == -65000:
            #                 gvel_invert[icount][vcount]=v
            #             vcount = vcount + 1
            #     icount =icount+1

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
            plt.imshow(y_axis_2,
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
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(np.flip(dap_vel,1), origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm, vmin=-250, vmax=250)
            plt.colorbar()
            plt.title(eml[j] + 'inv-GVEL', fontsize=10)

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
