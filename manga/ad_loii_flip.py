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
from scipy.stats import linregress
from statistics import mean
from PyPDF2 import PdfFileMerger
import copy
from copy import deepcopy
from os.path import isdir, join
import math

# function for plotting maps of relevant quantities
def daplot(quantity, qmin, qmax):
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    plt.imshow(quantity, origin='lower', interpolation='nearest',
        norm=colors.LogNorm(vmin=qmin, vmax=qmax), cmap=cm.coolwarm)
    plt.colorbar()

def circular_plot():
    pos_y = y_axis > 0
    neg_y = y_axis < 0
    pos_x = x_axis > 0
    neg_x = x_axis < 0
    gvel_invert = np.zeros([size, size])
    backg_val = abs(dap_vel[0][0])
    gvel_invert = gvel_invert - backg_val

    newx = copy.deepcopy(x_axis)
    newx = newx[np.logical_and(newx != 0, newx != 0)]
    newxt = copy.deepcopy(newx)
    newxt = np.transpose(newxt)
    newxt = newxt[np.logical_and(newxt != 0, newxt != 0)]

    yind, xdum = np.where(x_axis == newx[0])
    ydum, xind = np.where(x_axis == newxt[0])
    yinval = yind[0]
    xinval = xind[0]

    # figuring out the middle of the plot
    midy = int(int(len(dap_vel)) / 2)
    radius = midy - yinval
    mask_1 = np.zeros([size, size])
    mask_1 = mask_1 + 1

    # creating circle; masking out values that are further away than the radius previously defined
    for i in range(0, len(dap_vel)):
        for valz in range(0, len(dap_vel[i])):
            dist = ((midy - i) ** 2 + (midy - valz) ** 2) ** 0.5
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
        int(s)
    except OverflowError:
        return False
    return True


# calculating a second angle
def method_2():

    gvel_invert = copy.deepcopy(dap_vel)
    gvel_invert = np.flip(gvel_invert, 1)

    center = round(size / 2)
    arbit_vel = dap_vel[center + 5, center + 5]
    vy, vx = np.where(np.flip(dap_vel, 1) == dap_vel[center + 5, center + 5])
    x_fl = abs(vx[0] - center)
    y_fl = abs(vy[0] - center)
    theta_or = np.abs(np.arctan(1))
    theta_fl = np.abs(np.arctan(y_fl / x_fl))
    if ((vx[0] - center) > 0) and (abs(vy[0] - center) > 0):
        phi_2 = theta_or - theta_fl
    if ((vx[0] - center) < 0) and (abs(vy[0] - center) > 0):
        phi_2 = -1 * (np.pi - (theta_or + theta_fl))
    if ((vx[0] - center) < 0) and (abs(vy[0] - center) < 0):
        phi_2 = -1 * ((np.pi / 2) - theta_or + (np.pi / 2) + theta_fl)
    if ((vx[0] - center) > 0) and (abs(vy[0] - center) < 0):
        phi_2 = -1 * (theta_or + theta_fl)

    phi_2 = np.degrees(phi_2)

    
good_galaxies = [[8262, 1902, 1.2380620297469431, 6], [8597, 6104, 1.2339423537754584, 14],
                 [8987, 1901, 1.0405619373189734, 5], [7975, 9102, 0.94207652187996993, 16], # 2
                 [8313, 1901, 0.87610030276181694, 5], [8082, 3701, 0.82404316225203966, 7],
                 [8325, 3701, 0.81377880208120568, 7], [7990, 12703, 0.76183285172580029, 2], # 1
                 [8081, 1901, 0.75668569935436836, 5], [8132, 9101, 0.72775674210016394, 15],
                 [8465, 9101, 0.72703904475936965, 15], [8249, 1901, 0.72413077462329156, 5],
                 [8318, 12702, 0.68223025124789782, 1], [8711, 1902, 0.67447785224117618, 6], # 2
                 [8485, 6103, 0.64804992508574988, 13], [7977, 12704, 0.64478325464597885, 3],
                 [8552, 12702, 0.63332802500860841, 1], [8077, 6102, 0.61166555745912632, 12],
                 [8455, 12705, 0.60665902350021139, 4], [8606, 9102, 0.60095230121318499, 16], # 2 !
                 [8335, 3704, 0.5921413295274518, 10], [8082, 3702, 0.55895498851102443, 8],
                 [8601, 9101, 0.55708210155470717, 15], [7495, 3702, 0.55397262756957333, 8],
                 [8456, 12702, 0.54878300427726157, 1], [8462, 9102, 0.54659457408426393, 16],
                 [8333, 12703, 0.54564785156918072, 2], [8332, 9102, 0.53994760038850864, 16],
                 [8263, 12701, 0.52385177250381487, 0], [8712, 6103, 0.5222240329778487, 13],
                 [10001, 12704, 0.51866549790105843, 3], [8325, 12702, 0.51780380431022655, 1],
                 [8552, 12701, 0.51443301465537905, 0], [8145, 12701, 0.51000653250695349, 0],
                 [8341, 3702, 0.50556577847344275, 8], [8241, 12703, 0.50072291672206792, 2], # 1
                 [8250, 12705, 0.49687045905805932, 4], [8551, 3701, 0.49274267322591198, 7],
                 [8466, 9102, 0.49219758380982226, 16], [8341, 6102, 0.49134950015833856, 12], # 1
                 [8713, 12703, 0.46926664022890985, 2], [8548, 12704, 0.46613522145060393, 3],
                 [8145, 1902, 0.45785589001837146, 6], [8946, 3702, 0.45166970013429186, 8],
                 [8458, 3701, 0.44909872031265419, 7], [8623, 12702, 0.44645541008106382, 1],
                 [8320, 12703, 0.44334656699529668, 2], [8440, 6102, 0.4406229238036044, 12],
                 [8448, 12704, 0.41875562524863841, 3], [7992, 12703, 0.41691579663955919, 2],
                 [8943, 1901, 0.41668537719354126, 5], [8604, 12702, 0.41536295049237021, 1],
                 [8132, 3702, 0.40979483311293313, 8], [8719, 9102, 0.4058148608155791, 16],
                 [8149, 12703, 0.40352727567485019, 2], [8320, 3702, 0.40350785310729842, 8],
                 [8616, 12705, 0.40216556404162851, 4], [8612, 12704, 0.40119415451056462, 3],
                 [8082, 12703, 0.39749885964028359, 2], [8549, 12704, 0.39563819051111127, 3],
                 [8250, 6104, 0.39376075834508173, 14], [8448, 6102, 0.3886348900788561, 12],
                 [8712, 6104, 0.38858696255156172, 14], [8603, 6103, 0.38047044025724541, 13],
                 [8718, 1902, 0.37933473945768748, 6], [8078, 9101, 0.37612827026536461, 15],
                 [8313, 6104, 0.36440368407641727, 14], [8261, 12703, 0.36221469388166277, 2],
                 [8259, 12703, 0.36028586082488018, 2], [8612, 12703, 0.35937413277483754, 2],
                 [8611, 12703, 0.35627483233325941, 2], [8611, 9102, 0.35523592468213228, 16],
                 [8715, 6102, 0.35217193444625589, 12], [8081, 12704, 0.35079074602751215, 3],
                 [8256, 12701, 0.35024347473059553, 0], [8947, 3701, 0.34942036752239586, 7],
                 [8252, 9101, 0.34532877154005981, 15], [8252, 12703, 0.34454780518737216, 2],
                 [8258, 12703, 0.33444533897881396, 2], [8449, 12703, 0.33444228033386136, 2],
                 [8140, 6104, 0.33101584433969111, 14], [8138, 6101, 0.32832089131312514, 10],
                 [8595, 1902, 0.32573236786465798, 6], [7815, 12705, 0.32466227232060801, 4],
                 [8143, 12703, 0.32288236586503721, 2], [8255, 12704, 0.32215391411148675, 3],
                 [8138, 1901, 0.31822262838105059, 5], [7495, 12702, 0.31440585696262957, 1],
                 [7975, 1902, 0.31333692320892431, 6], [8552, 6102, 0.31286699241033489, 12],
                 [8156, 12705, 0.31268269817353039, 4], [8612, 1902, 0.31040898848055454, 6],
                 [8987, 9102, 0.30370523955520434, 16], [8082, 9102, 0.30248772208171665, 16],
                 [8274, 12701, 0.30240600075983959, 0], [8335, 12703, 0.30135736187721651, 2],
                 [8461, 12705, 0.30118073359538072, 4], [8588, 12705, 0.30116987331950157, 4],
                 [7815, 6101, 0.29971926932590376, 11], [8249, 12701, 0.29895664066967109, 0],
                 [8259, 12702, 0.29762218230251192, 1], [8948, 3703, 0.2966884328247012, 9],
                 [8084, 12705, 0.2955649285238951, 4], [8149, 6102, 0.29463448700424372, 12],
                 [8262, 3704, 0.29215428934451043, 10], [8987, 6101, 0.28849693471158805, 11],
                 [8719, 12705, 0.28825474040551458, 4], [8453, 9101, 0.28775962926045012, 15],
                 [8946, 6103, 0.2877208426018914, 13], [8455, 9101, 0.28650817876732992, 15],
                 [8464, 6103, 0.28383866792175866, 13], [7957, 6104, 0.28328708127569058, 14],
                 [10001, 3703, 0.27757769906511587, 9], [8439, 6103, 0.27730649502509203, 13],
                 [8138, 12702, 0.27631268432885536, 1], [8553, 12705, 0.27535414410516479, 4],
                 [8597, 12705, 0.26929770166779371, 4], [10001, 12701, 0.26876420003406543, 0],
                 [8720, 12704, 0.26805706392909173, 3], [8341, 6103, 0.26751988374891889, 13],
                 [8320, 9101, 0.26695811307338391, 15], [8263, 12702, 0.26587071251098537, 1],
                 [9042, 12701, 0.26561342054825515, 0], [8081, 9102, 0.26407112735283156, 16],
                 [7443, 6103, 0.26360463667661793, 13], [8980, 9102, 0.26352073659300518, 16],
                 [8261, 12705, 0.26240454340594949, 4], [8724, 12705, 0.25875275593591573, 4],
                 [8950, 12705, 0.2583784990196451, 4], [8942, 6103, 0.25695624416265622, 13],
                 [8141, 12701, 0.25629379652586043, 0], [8084, 3701, 0.25525551019225728, 7],
                 [8454, 9102, 0.25481182567349442, 16], [8147, 12705, 0.25415052014129241, 4],
                 [7991, 6103, 0.25015819178772625, 13], [8555, 12703, 0.24929274436012691, 2],
                 [8934, 12703, 0.24818905251821313, 2], [8258, 12702, 0.24633634351117289, 1],
                 [8262, 3701, 0.24559540344491718, 7], [8485, 6102, 0.24398491242784581, 12],
                 [8455, 1902, 0.24259605582790048, 6], [8257, 12702, 0.2423959267636443, 1],
                 [8611, 6104, 0.24039425327370323, 14], [8257, 12705, 0.23798579083944751, 4],
                 [7975, 12703, 0.23750348280683087, 2], [8341, 12701, 0.23429061295890391, 0],
                 [9049, 9101, 0.23403341841580871, 14], [8551, 12704, 0.23295406001910338, 3],
                 [8948, 12701, 0.23287630982921015, 0], [8567, 6102, 0.23202744014786006, 12],
                 [8135, 12705, 0.23075749678707899, 4], [8255, 12701, 0.2293615230425449, 0],
                 [8728, 9101, 0.22888000181898127, 15], [8547, 12702, 0.22867679893701109, 1],
                 [8081, 12703, 0.22746119341691032, 2], [8466, 12701, 0.22426125076937581, 0],
                 [8717, 12705, 0.22163309208751486, 4], [8078, 12704, 0.22126958546293668, 3],
                 [8262, 12702, 0.22121393196967182, 1], [8086, 3703, 0.21755542382551799, 9],
                 [8254, 12702, 0.21486779168359407, 1], [8456, 12705, 0.21064241069995532, 4],
                 [8726, 12703, 0.20796209161341328, 2], [8439, 9101, 0.20770405917525486, 15],
                 [8987, 6104, 0.20468031961645647, 14], [7962, 12705, 0.20415471120735104, 4],
                 [8464, 12703, 0.20227199120915743, 2], [7960, 12702, 0.20193348145973733, 1],
                 [8553, 9101, 0.20098880010776266, 15]]
    
# minimum and maximum emission-line fluxes for plot ranges
fmin = 1e-19
fmax = 1e-16

# define a standard cosmology
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

# directory for relevant Data Reduction Pipeline (drp) and Data Analysis Pipeline (DAP) information
# mpl5_dir = os.environ['MANGADIR_MPL5']  # Be aware that this directory syntax might need revision
# drp = fits.open(mpl5_dir + 'drpall-v2_0_1.fits')
mpl8_dir = os.environ['MANGADIR_MPL8']  # Be aware that this directory syntax might need revision

drp = fits.open(mpl8_dir + 'drpall-v2_5_3.fits')    # read in information from drpall file

drpdata = drp[1].data

# specific the plate of interest
from os.path import isdir, join



# for bookkeeping purposes, here's an array of emission-line names
eml = ['OIId---3728', 'Hb-----4862', 'OIII---4960', 'OIII---5008', 'OI-----6302', 'OI-----6365', 'NII----6549',
   'Ha-----6564', 'NII----6585', 'SII----6718', 'SII----6732']

filename ='scrap.pdf' #'ad_loii_flip.pdf'

with PdfPages(filename) as pdf:
    for i in range(0,50):#len(good_galaxies)):#15):
        #all_plates = [f for f in os.listdir(mpl5_dir + 'SPX-GAU-MILESHC/') if
        #    isdir(join(mpl5_dir + 'SPX-GAU-MILESHC/', f))]
        galaxy_c = good_galaxies[i]
        gal_indx = galaxy_c[3]
        plate = galaxy_c[0]
        ifu = galaxy_c[1]
        denom_and_rating = str(galaxy_c[0]) + '-' + str(galaxy_c[1]) + " Rating: " + str(galaxy_c[2])
        # lines = glob.glob(mpl5_dir + 'SPX-GAU-MILESHC/*/*/*' + str(plate) + '*MAPS-SPX-GAU-MILESHC.fits*')
        # lines = glob.glob(mpl8_dir + 'HYB10-MILESHC-MILESHC/' + str(plate) + '/' + \
        #        str(ifu) + '/manga-' + str(plate) + '-' + \
        #        str(ifu) + '-MAPS-HYB10-MILESHC-MILESHC.fits.gz')
        # name = lines[gal_indx]
        name = mpl8_dir + 'HYB10-MILESHC-MILESHC/' + str(plate) + '/' + \
               str(ifu) + '/manga-' + str(plate) + '-' + \
               str(ifu) + '-MAPS-HYB10-MILESHC-MILESHC.fits.gz'

        # THIS PART NEEDS ADJUSTING. THERE IS ONE PLATE WITH MORE THAN 4 NUMBERS; done
        # parse the plate and galaxy information from the filename
        galaxy = galaxy_c[1]
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
            x_axis_2 = np.zeros([size, size])
            x_axis_2 = x_axis_2 - 2
            for i in range(0, size):
                for k in range(0, size):
                    if (dap_ha_sivar[i, k] > 0.):
                        x_axis_2[i, k] = xproj_kpc[i, k] / u.kpc
            x_axis_2[0] = x_axis_2[0] - 2
            y_axis_2 = np.zeros([size, size])
            y_axis_2 = y_axis_2 - 2
            for i in range(0, size):
                for k in range(0, size):
                    if (dap_ha_sivar[i, k] > 0.):
                        y_axis_2[i, k] = yproj_kpc[i, k] / u.kpc


            y_axis = copy.deepcopy(zproj_kpc_map)


            # Drawing a line at x = 0
            xi, xo = np.where(np.around(x_axis_2,
                                        decimals=0) == 0.0)  # round x axis, find the zeros (should be a horizontal line through
            # the major axis
            xi_f, xo_f = np.where(np.around(np.flip(x_axis_2, 1), decimals=0) == 0.0)

            # mainx_slope = find_slope(xi,xo) #this is the slope of the line passing through the
            # 0 points in the xproj map, and the y-intercept of such line. This describes the line equation.
            # slope_f = find_slope(xi_f,xo_f)
            # alternative
            secnd_m_slope_o_y = -1/(linregress(np.where(np.around(y_axis_2, decimals=0) == 0.0)))[0]
            secnd_m_slope_f_y = -1 / (linregress(np.where(np.around(np.flip(y_axis_2, 1), decimals=0) == 0.0)))[0]




            mainx_slope = linregress(xo, xi)
            mainx_slope = mainx_slope[0]
            slope_f = linregress(xo_f, xi_f)
            slope_f = slope_f[0]
            if np.isinf(mainx_slope) or np.isinf(slope_f):
                print('Infinite slope for ', galaxy, '-', plate, ' Main slope isinf ', np.isinf(mainx_slope),
                      ' flip slope isinf ', np.isinf(slope_f))


            mainx_slope = (mainx_slope+secnd_m_slope_o_y)/2
            slope_f = (slope_f+secnd_m_slope_f_y)/2

            if mainx_slope > 0:
                arbit_val = 5
            if mainx_slope < 0:
                arbit_val = -5

            x_o = arbit_val / mainx_slope
            x_f = arbit_val / slope_f
            theta_o = np.abs(np.arctan(arbit_val / x_o))
            theta_f = np.abs(np.arctan(arbit_val / x_f))
            theta_o = np.abs(np.arctan(mainx_slope))
            theta_f = np.abs(np.arctan(slope_f))

            # if mainx_slope>0:
            phi = np.pi - (theta_o + theta_f)
            # if mainx_slope<0:
            #     phi = theta_o + theta_f
            phi = np.degrees(phi)

            test_val = np.amax(x_axis_2)
            test_val = int(test_val * 0.75)
            tst_o1, tst_o2 = np.where(np.around(x_axis_2, decimals=0) == test_val)
            indx_1 = int(len(tst_o1) / 3)
            indx_2 = int(len(tst_o2) * (2 / 3))
            tst_o = x_axis_2[tst_o1[indx_1], tst_o2[indx_1]]
            tst_oo = x_axis_2[tst_o1[indx_2], tst_o2[indx_2]]


            gvel_invert = copy.deepcopy(dap_vel)
            gvel_invert = np.flip(gvel_invert, 1)
            gvel_invert = scipy.ndimage.rotate(gvel_invert, phi, reshape=False)

            if mainx_slope < 0:
                phi = theta_f + theta_o
            tst_f = scipy.ndimage.rotate(np.flip(x_axis_2, 1), phi, reshape=False)[tst_o1[indx_1], tst_o2[indx_1]]
            tst_ff = scipy.ndimage.rotate(np.flip(x_axis_2, 1), phi, reshape=False)[tst_o1[indx_2], tst_o2[indx_2]]

            gvel_inv_2 = copy.deepcopy(gvel_invert)
            if ((tst_o * tst_f) < 0) or ((tst_oo * tst_ff) < 0):

                gvel_invert = scipy.ndimage.rotate(gvel_invert, 180, reshape=False)
                print("Phi for", plate, " - ", galaxy, ": ", phi, " inverted.")



            print("Phi for", plate, " - ", galaxy, ": ", phi, ". Slopes. M: ", mainx_slope, ". F: ", slope_f, ".")

            # Plotting interesting galaxies
            # plot 1: galaxy coordinates in kpc
            # plot 1: galaxy coordinates in kpc
            assymetry_p = dap_vel - gvel_invert

            fig = plt.figure()
            ax = fig.add_subplot(3, 3, 1)
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(y_axis_2,
                       origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm)
            plt.colorbar()
            plt.title('y_axis_2', fontsize=10)
            plt.suptitle(name)

            # plot 2/3: emission-line flux vs vertical distance
            ax = fig.add_subplot(3, 3, 2)
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(x_axis_2,
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
                       cmap=cm.coolwarm, vmin=-125, vmax=125)
            plt.colorbar()
            plt.title(eml[j] + 'gvel_invert', fontsize=10)

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
                       cmap=cm.coolwarm, vmin=-125, vmax=125)
            plt.colorbar()
            plt.title(eml[j] + 'subtracted', fontsize=10)

            # plot 6: emission-line SFLUX signal-to-noise ratio
            ax = fig.add_subplot(3, 3, 6)
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(np.flip(dap_vel, 1), origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm, vmin=-125, vmax=125)
            plt.colorbar()
            plt.title(eml[j] + 'np.flip(dap_vel)', fontsize=10)

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
                       cmap=cm.coolwarm, vmin=-125, vmax=125)
            plt.colorbar()
            plt.title(eml[j] + ' dap_vel', fontsize=10)

            # plot 9: emission-line dispersion
            ax = fig.add_subplot(3, 3, 9)
            ax.set_xlim(0, size)
            ax.set_ylim(0, size)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
            plt.imshow(gvel_inv_2, origin='lower',
                       interpolation='nearest',
                       cmap=cm.coolwarm, vmin=-125, vmax=125)
            plt.colorbar()
            plt.title(eml[j] + 'gvel_inv_2', fontsize=10)

            pdf.savefig()
            plt.close()
            plt.close(fig)

os.system("open %s &" % filename)
