from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.constants import c as lightspeed

#Functions

def cb_match_lambda(lam,element,dec=1):
    where = ([],)
    while len(where[0]) == 0:
        where = np.where(np.around(lam, decimals=dec) == np.around(element, decimals=dec))
        dec-=1

    print('WHERE: ', where, ' TYPE: ', type(where))
    diff_arr = [np.abs(lam[x] - element) for x in where[0]]
    min_indx = np.where(diff_arr == np.amin(diff_arr))
    print('MIN_INDX ', min_indx[0])
    closest_position = where[0][min_indx[0][0]]

    return closest_position


def get_fits_data(filepath):
    fits_file = fits.open(filepath)
    coadd = fits_file['COADD'].data
    specobj = fits_file['SPECOBJ'].data
    spzline = fits_file['SPZLINE'].data
    linewave = spzline['LINEWAVE']
    linename = spzline['LINENAME']
    linez = spzline['LINEZ']
    return coadd,specobj,spzline,linewave,linename,linez

#read relevant quantities (pending - implement get_fits_data)
def get_quantities(filepath):
    #open file given filepath
    fits_file = fits.open(filepath)

    flux_values = fits_file['COADD'].data['flux']
    model_val = fits_file['COADD'].data['model']
    HBWAV = fits_file['SPZLINE'].data['LINEWAVE'][15]*(1+fits_file['SPZLINE'].data['LINEZ'][15])
    #note the magic number 15 (happens to be index of Hbeta)

    minw = fits_file['SPECOBJ'].data['WAVEMIN']
    maxw = fits_file['SPECOBJ'].data['WAVEMAX']

    wave = 10**fits_file['COADD'].data['loglam'] #wavelegnths corresponding to each flux value

    return [wave,flux_values,model_val,minw,maxw,HBWAV]

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
    em_wav = linewave[indx]*(1+linez[indx])
    spread = (vel*1000*em_wav)/lightspeed

    return spread, em_wav

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
    continuum = np.median(np.concatenate((flux[cont_low:em_lowlim],flux[em_maxlim:cont_high])))
    dlambda = [lam[x+1]-lam[x] for x in range(em_lowlim,em_maxlim+1)]
    print('DLAMBDA: ',dlambda)
    print('FLUX EMISSION',flux[em_lowlim:em_maxlim+1])
    prelim_flux_integral = np.sum(dlambda*flux[em_lowlim:em_maxlim+1])
    flux_val = prelim_flux_integral-continuum

    return flux_val

#vars
filepath = "/Users/cbradna/Documents/spec-0761-54524-0409.fits"
wfilepath = "C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hires/spec-0761-54524-0409.fits"
outputfile = "plotpdftest.pdf"

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

data = get_quantities(wfilepath)
#pdfplot_flux_deffit(data,outputfile) don't need to call this for figuring out flux
hb_spread, hb_wav = wspread(wfilepath,'H_beta',300)
hb_flux = get_flux(wfilepath,hb_wav,hb_spread)