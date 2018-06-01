# Modified from code written by Christian Bradna
# Aleks Diamond-Stanic
# May 30, 2018

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

def channel_dictionary(hdu, ext):
    """
    Construct a dictionary of the channels in a MAPS file.
    """
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict
    
# define a standard cosmology
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

# directory for relevant Data Reduction Pipeline (drp) and Data Analysis Pipeline (DAP) information
mpl7_dir = os.environ['MANGADIR_MPL7']  # Be aware that this directory syntax might need revision
drp = fits.open(mpl7_dir + 'drpall-v2_4_3.fits')
drpdata = drp[1].data

# minimum and maximum emission-line fluxes for plot ranges
fmin = 1e-19
fmax = 1e-16

# for bookkeeping purposes, here's an array of emission-line names
eml = ['OIId---3728', 'Hb-----4862', 'OIII---4960', 'OIII---5008', 'OI-----6302', 'OI-----6365', 'NII----6549',
       'Ha-----6564', 'NII----6585', 'SII----6718', 'SII----6732']

# specify initial information about plates and galaxies to consider
all_plates = [f for f in os.listdir(mpl7_dir + 'HYB10-GAU-MILESHC/') if isdir(join(mpl7_dir + 'HYB10-GAU-MILESHC/', f))]
print('There are ', len(all_plates), ' plates in total.')
start_p = int(input("Start?")) - 1
fin_p = int(input("Finish?"))
min_p = input('Min inclination?')
max_p = input('Max inclination?')

# Fetching interesting galaxies and rating them
good_galaxies = []
gg_properties = []
prop_names_i = ["Plate", "Galaxy", "Global Quality", "Oii_ha Q 1", "Oii_ha Q 2", "Oii Q 1", "Oii Q 2",
                "Len of flux elements", "Mean of Flux",
                "Oii_ha Q1 SUM", "Oii_ha Q1 Mean", "Vertical Max", "S/N Sum", "S/N Mean", "Ha flux SUM", "Ha Mean",
                "B_Index SUM", "B_Index Mean",
                "Flux Weight Mean", "Flux Ratio Weight Mean", "Vertical Weight Value"]

filename = 'MPL7_dap_multi_automated_quicklook.pdf'
with PdfPages(filename) as pdf:
    # loop over all plates
    for each_plate in range(start_p, fin_p):
        plate_rec = all_plates[each_plate]
        
        min_par = float(min_p)
        max_par = float(max_p)
    
        # array to store interesting galaxies
        interesting_main = []
        properties = []
    
        # specific the plate of interest
        plate = plate_rec
        lines = glob.glob(mpl7_dir + 'HYB10-GAU-MILESHC/*/*/*' + str(plate) + '*MAPS-HYB10-GAU-MILESHC.fits*')

        # loop over each galaxy
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
    
            end = name.find('-MAPS-HYB10-GAU-MILESHC')
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
            emlc = channel_dictionary(hdu_dap, 'EMLINE_SFLUX')
            dap_ha_sflux = hdu_dap['EMLINE_SFLUX'].data[emlc['Ha-6564'],:,:]
            dap_ha_sivar = hdu_dap['EMLINE_SFLUX_IVAR'].data[emlc['Ha-6564'],:,:]
            
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
    
            #insert IF statement for inclination
            if (min_par < (inc * 180 / np.pi) < max_par):
    
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
                    # We are working with oii_flux, oii to ha, vertical distance, inclination, masks and the error for the S/N
    
                    #trying to take the masks, creating a masking array that takes awful values and turns them into 0,
                    spax_masks = copy.deepcopy(dap_mask)
                    # for g in spax_masks:
                    #     for f in g:
                    #         if f >= (2**30):
                    #             f = 0
                    # spax_masks = spax_masks*(1/(2**32))
                    # spax_masks = (1/spax_masks)
                    good_flag = (spax_masks<(2**30))
                    spax_masks[good_flag] = 1
                    bad_spaxels = (spax_masks >= (2**30))
                    spax_masks[bad_spaxels] = 0
    
                    # Filtering values nearest to the galaxy disk according to a filtering parameter
    
                    # filt_parm = ((90 / (inc * 180 / np.pi)) ** 2) * (np.amax(zproj_kpc_map) / ((inc * 180 / np.pi) ** 0.4)) LEVEL 1
                    r_e = (drpdata.field('nsa_sersic_th50')[match][0] * u.arcsec)/cosmo.arcsec_per_kpc_proper(z)
                    Vr_e = r_e * np.cos(inc)
    
                    filt_parm = ((90 / (inc * 180 / np.pi))**2 * (Vr_e))/u.kpc
    
                    disk_filter = copy.deepcopy(zproj_kpc_map)
    
                    bad = (np.abs(disk_filter) < filt_parm)
                    disk_filter[bad] = 0
                    disk_filter = np.abs(disk_filter)
                    # rest = (vertical_distance_map != 0)
                    # vertical_distance_map[rest] = 1
                    disk_filter = disk_filter/disk_filter
                    nans = np.isnan(disk_filter)
                    disk_filter[nans] = 0
                    disk_filter = disk_filter*spax_masks
    
                    # LEVEL 2 Disk Filtering
                    filt_parm = ((90 / (inc * 180 / np.pi)) ** 3 * (Vr_e)) / u.kpc
    
                    disk_filter_2 = copy.deepcopy(zproj_kpc_map)
    
                    bad = (np.abs(disk_filter_2) < filt_parm)
                    disk_filter_2[bad] = 0
                    disk_filter_2 = np.abs(disk_filter_2)
                    # rest = (vertical_distance_map != 0)
                    # vertical_distance_map[rest] = 1
                    disk_filter_2 = disk_filter_2 / disk_filter_2
                    nans = np.isnan(disk_filter_2)
                    disk_filter_2[nans] = 0
                    disk_filter_2 = disk_filter_2 * spax_masks
    
                    # Creating quality number oii_flux level 2
    
                    oii_flux_2 = (copy.deepcopy(dap_sflux)) * disk_filter_2
                    useful_spx_num = len(oii_flux_2[(oii_flux_2 != 0)])
                    bad = (0.3 > oii_flux_2)
                    oii_flux_2[bad] = 0
                    bad = (oii_flux_2 > 3)
                    oii_flux_2[bad] = 0
                    inter_spx_num = len(oii_flux_2[(oii_flux_2 != 0)])  # C parameter in my notes
                    if (useful_spx_num != 0):
                        oii_quality_2 = inter_spx_num / useful_spx_num  # Quality number
                    if (useful_spx_num == 0):
                        oii_quality_2 = 0
    
                    # Creating quality number oii_flux level 1
    
                    oii_flux = (copy.deepcopy(dap_sflux)) * disk_filter
                    useful_spx_num = len(oii_flux[(oii_flux != 0)])
                    bad = (0.25 > oii_flux)
                    oii_flux[bad] = 0
                    bad = (oii_flux > 5)
                    oii_flux[bad] = 0
                    inter_spx_num = len(oii_flux[(oii_flux != 0)])  # C parameter in my notes
                    if (useful_spx_num != 0):
                        oii_quality = inter_spx_num / useful_spx_num  # Quality number
                    if (useful_spx_num == 0):
                        oii_quality = 0
    
                    # Masking out ratio plot. Removing x < 1 (and nans). Creating quality number oii_ha level 2
                    oii_to_ha_2 = (dap_sflux / dap_ha_sflux)
    
                    nans = np.isnan(oii_to_ha_2)
                    oii_to_ha_2[nans] = 0
                    infs = np.isinf(oii_to_ha_2)
                    oii_to_ha_2[infs] = 0
    
                    oii_to_ha_2 = oii_to_ha_2 * disk_filter_2
                    useful_spx_num = len(oii_to_ha_2[(oii_to_ha_2 != 0)])
    
                    bad = (oii_to_ha_2 < 0.7)
                    oii_to_ha_2[bad] = 0
                    unreas = (oii_to_ha_2 > 4)
                    oii_to_ha_2[unreas] = 0
                    unreas = (oii_to_ha_2 > ((1 + np.mean(oii_to_ha_2)) ** 2))
                    oii_to_ha_2[unreas] = 0
                    inter_spx_num = len(oii_to_ha_2[(oii_to_ha_2 != 0)])
                    if (useful_spx_num != 0):
                        oii_ha_quality_2 = inter_spx_num / useful_spx_num
                    if (useful_spx_num == 0):
                        oii_ha_quality_2 = 0
    
                    # Creating quality number oii_ha level 1
    
                    oii_to_ha = (dap_sflux / dap_ha_sflux)
    
                    nans = np.isnan(oii_to_ha)
                    oii_to_ha[nans] = 0
                    infs = np.isinf(oii_to_ha)
                    oii_to_ha[infs] = 0
    
                    oii_to_ha = oii_to_ha * disk_filter
                    useful_spx_num = len(oii_to_ha[(oii_to_ha != 0)])
    
                    bad = (oii_to_ha < 0.5)
                    oii_to_ha[bad] = 0
                    unreas = (oii_to_ha > ((1 + np.mean(oii_to_ha)) ** 3))
                    oii_to_ha[unreas] = 0
                    unreas = (oii_to_ha > 9)
                    oii_to_ha[unreas] = 0
                    inter_spx_num = len(oii_to_ha[(oii_to_ha != 0)])
                    if (useful_spx_num != 0):
                        oii_ha_quality_1 = inter_spx_num / useful_spx_num
                    if (useful_spx_num == 0):
                        oii_ha_quality_1 = 0
    
                    # Vertical distance weight calculator
                    ver_dist = copy.deepcopy(zproj_kpc_map)
                    ver_dist = ver_dist * disk_filter
                    ver_dist = np.abs(ver_dist)
                    if (useful_spx_num != 0):
                        ver_dist = (ver_dist - np.amin(ver_dist[np.nonzero(ver_dist)]))
                        ver_dist[ver_dist < 0] = 0
                        ver_dist = ver_dist * (10 / np.amax(ver_dist))
                    if (useful_spx_num == 0):
                        ver_dist = 0
    
                    # Vertical weighing
                    flux_vertical_w = oii_flux * (ver_dist ** 2)
                    ratio_vertical_w = oii_to_ha * (ver_dist ** 2)
                    flux_w_mean = np.mean(flux_vertical_w)
                    ratio_w_mean = np.mean(ratio_vertical_w)
                    v_weight = (flux_w_mean * ratio_w_mean) / 9
    
                    # Generating some numbers to rank
    
                    normalize_vdist = (inc * 180 / np.pi) / (2 * np.amax(zproj_kpc_map))
                    bradna_index = np.abs(((10 * oii_flux) ** 2) * ((10 * oii_to_ha) ** 3.5) * (
                        oii_flux / dap_serr) * disk_filter)  # * (normalize_vdist ** 2))
                    global_quality = ((oii_ha_quality_1 + oii_ha_quality_2 + oii_quality + oii_quality_2 + v_weight) / 5) * 2 * (np.mean(oii_to_ha))
                    # Calculating and printing stuff
                    param = (bradna_index > 0)
                    bradna_index_sum = np.sum(bradna_index[param])
                    if bradna_index_sum == 0:
                        bradna_index_sum = 10
    
                    if (np.sum(dap_sflux[dap_sflux > 0.2]) > 540) and (global_quality > 0.32):
                        state = '>>>>>>>>>>>>>>>>Has enough information'
                    else:
                        state = 'Not enough information'
                    b_index = math.log(bradna_index_sum, 10)
                    print(plate, galaxy, 'B_index:', b_index, state, 'Rating: ',b_index * (100 / 6))  # math.log(bradna_index_sum, 10)
                    print("Global Quality:", global_quality," Inclination: ",(inc * 180 / np.pi),'\n')
                    # denom_and_rating = str(denomination, ) + ' Rating: ' + str(
                    #     b_index * (100 / 62)) + ' Flux/rating ratio: ' + str((np.sum(oii_flux) / b_index))
    
                    # Appending interesting galaxies to an array
                    if state == '>>>>>>>>>>>>>>>>Has enough information':
                        gal_indx = lines.index(name)
                        # int_gal = np.array([[plate, galaxy, (b_index * (100 / 62)), gal_indx]])
                        # print(int_gal)
                        # np.concatenate((interesting_main, int_gal), axis=0)
                        int_gal = [plate, galaxy, global_quality, gal_indx]  #Changed number that gets the galaxies sorted from (b_index * (100 / 6) to global quality
                        # print(int_gal, ' Inclination: ',inc * 180 / np.pi)
                        ## TEMPORARY TESTING INFORMATIONAL PURPOSE
                        non_zero = oii_flux > 0
                        # print("Length of flux elements",len(oii_flux[non_zero])," Mean value of flux elements ",np.mean(oii_flux))
                        # print("Sum of oii-to-ha level 1 strain ",np.sum(oii_to_ha) ,' Mean: ',np.mean(oii_to_ha))
                        # print("Vertical max", np.amax(np.abs(zproj_kpc_map)))
                        # print("Sum of s/n",np.sum(oii_flux / dap_serr) ,' Mean: ',np.mean(oii_flux / dap_serr))
                        # print("Sum of ha",np.sum(dap_ha_sflux) ,' Mean: ',np.mean(dap_ha_sflux))
                        # print("Index sum: ",bradna_index_sum,' bradna_index mean: ',np.mean(bradna_index))
                        # print("Flux weight mean: ",flux_w_mean,' Ratio weight mean :',ratio_w_mean)
                        # print("Oii Quality: ",oii_quality," Oii_ha Q level 1: ",oii_ha_quality_1," Oii_ha Q level 2",oii_ha_quality_2," Global quality: ",global_quality,'\n')
                        properties_i = [plate, galaxy, global_quality, oii_ha_quality_1, oii_ha_quality_2, oii_quality, oii_quality_2,
                                      len(oii_flux[non_zero]),np.mean(oii_flux),np.sum(oii_to_ha),np.mean(oii_to_ha),np.amax(np.abs(zproj_kpc_map)),
                                      np.sum(oii_flux / dap_serr),np.mean(oii_flux / dap_serr),np.sum(dap_ha_sflux),np.mean(dap_ha_sflux),
                                      bradna_index_sum,np.mean(bradna_index),flux_w_mean,ratio_w_mean,v_weight]
    
                        interesting_main.append(int_gal)
                        properties.append(properties_i)

                        # Plotting interesting galaxies
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
                        plt.title('vertical distance [kpc]', fontsize=8)
                        #plt.suptitle(denom_and_rating, fontsize=8)
                        plt.suptitle(str(plate)+'-'+str(galaxy)+' B_index:'+'{0:.2f}'.format(b_index))

                        # plot 2/3: emission-line flux vs vertical distance
                        ax = fig.add_subplot(3, 3, 2)
                        daplot(dap_sflux * spax_masks * 1.e-17, fmin, fmax)
                        plt.title('Flux w/ mask', fontsize=8)
                     
                     
                        # plot 2/3: emission-line flux vs vertical distance
                        ax = fig.add_subplot(3, 3, 3)
                        daplot(dap_sflux * disk_filter * 1.e-17, fmin, fmax)
                        plt.title('Flux*vertical', fontsize=8)
                     
                        # plot 4: emission-line SFLUX
                        ax = fig.add_subplot(3, 3, 4)
                        daplot(dap_sflux * 1.e-17, fmin, fmax)
                        plt.title(eml[j] + ' SFLUX', fontsize=8)
                     
                        # plot 5: emission-line SFLUX serror
                        ax = fig.add_subplot(3, 3, 5)
                        daplot(oii_to_ha, 0.1, 10.)
                        plt.title('Oii_to_ha Filtered', fontsize=8)
                     
                        # plot 6: emission-line SFLUX signal-to-noise ratio
                        ax = fig.add_subplot(3, 3, 6)
                        daplot(oii_flux_2* 1.e-17, fmin, fmax)
                        plt.title('Oii_2', fontsize=8)
                     
                        # plot 7: emission-line flux / Halpha flux
                        ax = fig.add_subplot(3, 3, 7)
                        daplot(dap_sflux / dap_ha_sflux, 0.1, 10.)
                        plt.title(eml[j] + '/HA', fontsize=8)
                     
                        # plot 8: emission-line velocity
                        ax = fig.add_subplot(3, 3, 8)
                        daplot(oii_flux * 1.e-17,fmin, fmax)
                        plt.title('Oii flux', fontsize=8)
                     
                        # (((10 * oii_flux) ** 2) * (oii_to_ha ** 2) * (oii_flux / dap_serr) * (vertical_distance_map) * spax_masks)
                     
                        # plot 9: emission-line dispersion
                        ax = fig.add_subplot(3, 3, 9)
                        daplot(oii_to_ha_2, 0.1, 10.)
                        plt.title('Oii_to_ha Filtered', fontsize=8)
                     
                     
                        pdf.savefig()
                        plt.close()
                        plt.close(fig)
                        
            else:
                print(plate,galaxy, " - Outside inclination threshold. \n")
    
        #I have to check if this np.delete is really necessary. Coming back this fall, I think I forgot what it does. But it breaks down when there are no single interesting galaxies in
        #a given plate. So I need to make it go through plates even if it doesn't find anything good.
        if len(interesting_main)>0:
            np.delete(interesting_main,0,axis=0)
    
        interesting = interesting_main
        ggproperty = properties 
        #good_galaxies = np.concatenate((good_galaxies, interesting), axis=0)
        good_galaxies = good_galaxies + interesting
        gg_properties = gg_properties + ggproperty
        #return (interesting_main,properties)
    


# Sorting interesting galaxies
# np.delete(good_galaxies,0,axis=0)
# good_galaxies = good_galaxies[np.argsort(good_galaxies[:, 2])]
good_galaxies.sort(key=lambda x: x[2], reverse=True)

# Plotting interesting galaxies in order
for i in good_galaxies:
    print('\n')
    #Name for ith galaxy
    ithname = str(i[0])+str(i[1])
    # Finding galaxy properties by comparing ith name
    for b in range(0, len(gg_properties)):
        bthname = str(gg_properties[b][0])+str(gg_properties[b][1])
        if ithname == bthname:
            element_index = b
            for c in range(0, len(gg_properties[b])):
                print(prop_names_i[c],' ',gg_properties[b][c])

print("KEEP IN MIND SOMETIMES YOU GET AN ERROR WHEN THE PROGRAM TRIES TO OPEN THE PDF. JUST GO TO THE DIRECTORY WHERE IT WAS CREATED AND OPEN IT MANUALLY.")
os.system("open %s &" % filename)
