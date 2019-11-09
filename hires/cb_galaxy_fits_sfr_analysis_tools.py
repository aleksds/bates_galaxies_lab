from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import importlib

from uncertainties.umath import *
from uncertainties import ufloat_fromstr
from uncertainties import ufloat
from matplotlib.backends.backend_pdf import PdfPages
from scipy.constants import c as lightspeed
from scipy.constants import parsec as prsc
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import matplotlib
import matplotlib.gridspec as gridspec

#Vars
ppcbfits = importlib.import_module('cb_ppxf_tools')


#Functions

#looks for a specific lambda in array, returns the position of the closest value
# def cb_match_lambda(lam,element,dec=1):
# def cb_match_lambda(lam,element,dec=1):
#     where = ([],)
#     while len(where[0]) == 0:
#         where = np.where(np.around(lam, decimals=dec) == np.around(element, decimals=dec))
#         dec-=1
#
#     #print('WHERE: ', where, ' TYPE: ', type(where))
#     diff_arr = [np.abs(lam[x] - element) for x in where[0]]
#     min_indx = np.where(diff_arr == np.amin(diff_arr))
#     #print('MIN_INDX ', min_indx[0])
#     closest_position = where[0][min_indx[0][0]]
#
#     return closest_position

def cb_match_lambda(lam, element):
    diff_arr = lam - element
    closest_pos = np.where(np.abs(diff_arr) == np.amin(np.abs(diff_arr)))[0][0]

    return closest_pos

#gets some relevant data from the fits file (to avoid passing the data as argument as much as possible)
def get_fits_data(filepath):
    fits_file = fits.open(filepath)
    coadd = fits_file['COADD'].data
    specobj = fits_file[2].data #replacing 'SPECOBJ' with index 2 since galaxy J1506 doesn't have name SPECOBJ
    spzline = fits_file['SPZLINE'].data
    linewave = spzline['LINEWAVE']
    linename = spzline['LINENAME']
    linez = spzline['LINEZ']
    return coadd,specobj,spzline,linewave,linename,linez

#read relevant quantities
def get_quantities(filepath,emission_name):
    #open file given filepath
    coadd,specobj,spzline,linewave,linename,linez = get_fits_data(filepath)

    hb_ind = np.where(linename == emission_name)

    #HBWAV = spzline['LINEWAVE'][hb_ind]*(1+spzline['LINEZ'][hb_ind])
    HBWAV = spzline['LINEWAVE'][hb_ind]*(1+ppcbfits.get_z(filepath))

    minw = specobj['WAVEMIN']
    maxw = specobj['WAVEMAX']

    wave = 10**coadd['loglam'] #wavelegnths corresponding to each flux value

    return [wave,coadd['flux'],coadd['model'],minw,maxw,HBWAV]

#plot flux and default fit/model already embedded in fits
def pdfplot_flux_deffit(data,outputfile):
    with PdfPages(outputfile) as pdf:
        wave,flux_values,model_val,minw,maxw,HBWAV = data[0],data[1],data[2],data[3],data[4],data[5]

        fig, ax = plt.subplots(dpi=1200)

        plt.axvline(x=HBWAV, color='blue', label='Hb', alpha=0.4)

        ax.plot(wave,flux_values,color='black',linewidth=0.3)
        ax.plot(wave,model_val,color='red',linewidth=0.3, alpha=0.7)

        plt.grid(True)
        plt.xlim(minw,maxw)
        ax.xaxis.set_major_locator(plt.MultipleLocator(500))
        plt.xlabel("$\AA ngstr \ddot{o} ms$")
        plt.ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")
        plt.legend()

        pdf.savefig()
        plt.close('all')

#gets a wavelength spread value based on velocity (km/s); also returns the value of the wavelength of the emission
#that we expect to see in the actual spectrum
def wspread(filepath,emission,vel):
    coadd,specobj,spzline,linewave,linename,linez = get_fits_data(filepath)
    indx = np.where(linename == emission)
    #not sure this one is correct
    # em_wav = linewave[indx]*(1+linez[indx])
    em_wav = linewave[indx]*(1+ppcbfits.get_z(filepath))
    spread = (vel*1000*em_wav)/lightspeed

    return spread, em_wav

#gets the title for the figure based on the file that is being plotted
def get_plot_title(filepath):
    fits_file = fits.open(filepath)

    file_char_list = list(filepath)
    spec_pos = filepath.find('spec-')
    title_list = file_char_list[spec_pos:-5]
    title = ''.join(title_list)
    headers = fits_file['PRIMARY'].header
    hh = str(int(headers['PLUG_RA'] * (24/360)))
    if not (len(hh) == 2):
        hh = '0' + hh
    mm = str(int((headers['PLUG_RA'] * (24/360) * 60) % 60))
    if not (len(mm) == 2):
        mm = '0' + mm
    ss = str(np.round((headers['PLUG_RA'] * (24/360) * 60 * 60) % 60, decimals=2))
    if not (len(str(int(np.round((headers['PLUG_RA'] * (24/360) * 60 * 60) % 60)))) == 2):
        ss = '0' + ss
    title = title + '/J' + hh + mm + ss
    return title

#plots a wavelength range, shows continuum constant and the range over which the emission is calculated to extend
def plot_area_of_interest(plot_data,txt_data,filepath,ppxf_data=False):
    matplotlib.rcParams.update({'font.size': 12})

    #Get relevant stuff for generating a figure
    em_lowlim,em_maxlim,cont_low,cont_high,flux,lam,cont_const = plot_data[0],plot_data[1],plot_data[2],plot_data[3],\
                                                                 plot_data[4],plot_data[5],plot_data[6]
    hbflux,hblum,sfr_ufloat = txt_data[0],txt_data[1],txt_data[2]
    sfr = sfr_ufloat.n
    sfr_u = sfr_ufloat.s
    halum = hblum*2.468
    lam_c = lam[cont_low:cont_high+1]
    flux_c = flux[cont_low:cont_high+1]
    plot_title = get_plot_title(filepath)

    fluxdata = cont_const[2]

    #Generating that figure
    #with PdfPages('Emissionline.pdf') as pdf:
    fig = plt.figure()
    fig.suptitle(plot_title)

    ax = fig.add_subplot(1, 2, 1)
    ax.axvspan(lam[em_lowlim], lam[em_maxlim], alpha=0.5, color='red', label='Area of Emission')
    ax.plot(lam_c, flux_c, color='black', linewidth=0.3)

    if ppxf_data:
        galfit_c = cont_const[0][cont_low:cont_high+1]
        gasfit_c = cont_const[1][cont_low:cont_high+1]
        ax.plot(lam_c, galfit_c, color='orange', linewidth=0.3)
        ax.plot(lam_c, gasfit_c, color='red', linewidth=0.3)
        ax.plot(lam_c, galfit_c-gasfit_c, color='blue', linewidth=0.3)
        ax.plot(lam_c, (galfit_c-gasfit_c)*np.median(fluxdata), color='blue', linewidth=0.3)
        ax.plot(lam_c, fluxdata[cont_low:cont_high+1], color='grey', linewidth=0.3)
    elif not ppxf_data:
        plt.axhline(y=cont_const, color='blue', label='Continuum', alpha=0.4)


    ax.set_xlabel("$\AA ngstr \ddot{o} ms$")
    ax.set_ylabel("Flux [$10^{-17}$ erg/$cm^{2}$/s/$\AA$]")
    plt.legend(loc=2)
    ax.grid(True)

    ax2 = fig.add_subplot(6, 2, 4)
    ax2.text(0,1,'H beta flux: '+str(hbflux), fontsize=11,wrap=True)
    ax2.axes.axis('off')

    ax3 = fig.add_subplot(6, 2, 6)
    ax3.text(0,1,'H beta luminosity: '+str(hblum), fontsize=11,wrap=True)
    ax3.axes.axis('off')

    ax4 = fig.add_subplot(6, 2, 8)
    ax4.text(0,1,'H alpha luminosity: '+str(halum), fontsize=11,wrap=True)
    ax4.axes.axis('off')

    ax5 = fig.add_subplot(6, 2, 10)
    ax5.text(0,1,'Star formation rate: '+str(sfr_ufloat), fontsize=11,wrap=True)
    ax5.axes.axis('off')


        #pdf.savefig(bbox_inches="tight")
        #plt.close('all')

#gets flux value (does the sums and stuff)
def get_flux(filepath,wavel,wspread):
    coadd, specobj, spzline, linewave, linename, linez = get_fits_data(filepath)
    flux = coadd['flux']
    lam = 10**coadd['loglam']
    #next we get indeces to mark regions over which we will sum up
    em_center = cb_match_lambda(lam,wavel)
    em_lowlim = cb_match_lambda(lam,wavel - wspread)
    em_maxlim = cb_match_lambda(lam,wavel + wspread)
    cont_high = cb_match_lambda(lam,wavel + 4*wspread)
    cont_low = cb_match_lambda(lam,wavel - 4*wspread)

    #continuum = np.mean([np.mean(flux[cont_low:em_lowlim]),np.mean(flux[em_maxlim:cont_high])])
    continuum_const = np.median(np.concatenate((flux[cont_low:em_lowlim],flux[em_maxlim:cont_high])))
    dlambda = np.array([lam[x+1]-lam[x] for x in range(em_lowlim,em_maxlim+1)])
    continuum_area = np.sum(float(continuum_const)*dlambda)
    prelim_flux_integral = np.sum(dlambda*flux[em_lowlim:em_maxlim+1])
    flux_val = prelim_flux_integral-continuum_area

    ##print tests
    #print(lam[cont_low:cont_high])

    plot_data = [em_lowlim,em_maxlim,cont_low,cont_high,flux,lam,continuum_const]


    return flux_val, plot_data

#function to calculte flux
def get_lum(filepath,flux_value):
    current_z = ppcbfits.get_z(filepath)
    lum_dis = (cosmo.luminosity_distance(current_z))*(1/u.Mpc)*(100*(prsc*10**6))
    #print('LUMDIST',lum_dis)
    luminosity = (flux_value*10**-17)*4*np.pi*(lum_dis**2)

    return luminosity

#function got get sfr
def get_sfr_hbeta(hbet_lum_u):
    #print('HBETU',hbet_lum_u)
    hbet_lum = ufloat_fromstr(str(hbet_lum_u))
    halpha_lum = hbet_lum*2.468
    #print('HALPHALUM TYPE ',type(halpha_lum))
    if halpha_lum > 0:
        sfrlog = log10(halpha_lum)-41.27
        sfr = 10 ** sfrlog

    elif halpha_lum <= 0:
        sfrlog = -1000
        print("Error: Negative luminosity; corresponding sfr has been set to 1/1000")
        sfr = ufloat(0,0)
    else:
        print("--------------->Halpha lum isn't a number. Type: ",type(halpha_lum))
        sfr = ufloat(0,0)



    return sfr



###Below is code written for testing these tools while writing/debugging


#vars
# filepath = "/Users/cbradna/Documents/spec-0761-54524-0409.fits"
# wfilepath = "C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hires/galaxies_fits/spec-0761-54524-0409.fits"
# outputfile = "plotpdftest.pdf"

# Names of emission lines as they appear inside the fits file
# ['Ly_alpha', 'N_V 1240', 'C_IV 1549', 'He_II 1640',
#            'C_III] 1908', 'Mg_II 2799', '[O_II] 3725', '[O_II] 3727',
#            '[Ne_III] 3868', 'H_epsilon', '[Ne_III] 3970', 'H_delta',
#            'H_gamma', '[O_III] 4363', 'He_II 4685', 'H_beta',
#            '[O_III] 4959', '[O_III] 5007', 'He_II 5411', '[O_I] 5577',
#            '[O_I] 6300', '[S_III] 6312', '[O_I] 6363', '[N_II] 6548',
#            'H_alpha', '[N_II] 6583', '[S_II] 6716', '[S_II] 6730',
#            '[Ar_III] 7135'

#Run routine

# data = get_quantities(wfilepath)
# # pdfplot_flux_deffit(data,outputfile) don't need to call this for figuring out flux
# hb_spread, hb_wav = wspread(wfilepath,'H_beta',400)
# hb_flux, plot_data = get_flux(wfilepath,hb_wav,hb_spread)
# hb_lum = get_lum(wfilepath,hb_flux)
# sfr = get_sfr_hbeta(hb_lum)
#
# with PdfPages('Emissionline.pdf') as pdf:
#     plot_area_of_interest(plot_data,[hb_flux,hb_lum,sfr],wfilepath)
#     pdf.savefig(bbox_inches="tight")
#     plt.close('all')
