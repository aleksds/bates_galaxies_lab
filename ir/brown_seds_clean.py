import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from matplotlib.ticker import ScalarFormatter
from sedpy import observate
from astropy.io import fits
from astropy.table import Table


def to_ab(vega_mag, band):
    '''Converts from vega magnitudes to ab mags '''
    zero_mag_fnu = {"W1": 309.540, "W2": 171.787, "W3": 31.674, "W4": 8.363}
    print("vega mag: ", vega_mag)
    f_nu = float(zero_mag_fnu[band] * (10 ** (-vega_mag / 2.5)))
    mag_ab = float(-2.5 * np.log10(f_nu / 3631))

    return mag_ab


def plot_phot(redshift, given_name):
    "plots redshift against w3-w4 for a given galaxy"

    print("Valid:", valid)
    ax.scatter(redshift[valid], w3_mag[valid] - w4_mag[valid], marker='*')
    ax.errorbar(redshift[valid], w3_mag[valid] - w4_mag[valid], yerr=w3minusw4_unc[valid], color='purple',
                linestyle="None")

    for i in range(0, len(redshift[valid])):
        plt.text(redshift[valid][i], w3_mag[valid][i] - w4_mag[valid][i], given_name[valid][i], fontsize=6)


def calc_w3_w4(redshift, given_name,projpath):
    """
    Gets w3,w4 values and uncertainties from unWISE and places them into lists. Then it converts all lists into
    numpy arrays so that they can be further manipulated
    """
    w3_mag = np.zeros(len(redshift))
    w3_mag_err = np.zeros(len(redshift))
    w4_mag = np.zeros(len(redshift))
    w4_mag_err = np.zeros(len(redshift))
    table = []
    count = 0

    # Converting to ab magnitude for each bandpass  data of each galaxy
    for g in given_name:
        print("gal name", g)
        WISEdir = projpath + 'unWISE/%s.fits' % g
        t = Table.read(WISEdir)
        table.append(t)
        w3_mag[count] = float(to_ab(t['w3_mag'], "W3"))
        w3_mag_err[count] = float(t['w3_mag_err'])
        w4_mag[count] = float(to_ab(t['w4_mag'], "W4"))
        w4_mag_err[count] = float(t['w4_mag_err'])

        count += 1

    redshift = np.asarray(redshift)
    given_name = np.asarray(given_name)
    w3_mag = np.asarray(w3_mag)
    w3_mag_err = np.asarray(w3_mag_err)
    w4_mag = np.asarray(w4_mag)
    w4_mag_err = np.asarray(w4_mag_err)
    w3minusw4_unc = np.sqrt(w3_mag_err ** 2 + w4_mag_err ** 2)
    print(w3minusw4_unc)
    valid = w3minusw4_unc < 0.5

    print(valid)

    return given_name, redshift, w3_mag, w3_mag_err, w4_mag, w4_mag_err, w3minusw4_unc, valid


### For Brown
projpath = os.getcwd() + '/Brown2014/An_Atlas_of_Galaxy_SEDs/An_Atlas_of_Galaxy_SEDs/'

###For our galaxies from Petter
# projpath = os.getcwd() + '/'
# print(projpath)

files = glob.glob(projpath + '*.dat')
print("len files", len(files))

def choose_dataset(dataset_name):
    """
    Choose which dataset you want to use and it'll populate the right information.
    """
    if dataset_name == "special 8":
        projpath = os.getcwd() + '/'
        redshift = [0.603, 0.711, 0.514, 0.467, 0.451, 0.661, 0.449, 0.459]
        special_names = ['J0826', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1613', 'J2118']
        special_names, redshift, w3_mag, w3_mag_err, w4_mag, w4_mag_err, w3minusw4_unc, valid = calc_w3_w4(redshift,
                                                                                                        special_names,projpath)
        return projpath, redshift, w3_mag, w3_mag_err, w4_mag, w4_mag_err, w3minusw4_unc, valid

    elif dataset_name == "jwst":
        print("yes it does")
        projpath = os.getcwd() + '/'
        redshift = [0.467, 0.451, 0.608, 0.449, 0.459]
        jwst_names = ['J1107', 'J1219', 'J1506', 'J1613', 'J2118']
        jwst_names, redshift, w3_mag, w3_mag_err, w4_mag, w4_mag_err, w3minusw4_unc, valid = calc_w3_w4(redshift,
                                                                                                        jwst_names,projpath)
        return projpath,jwst_names, redshift, w3_mag, w3_mag_err, w4_mag, w4_mag_err, w3minusw4_unc, valid

    elif dataset_name == "original":
        proj_path = os.getcwd() + '/'
        "This is maybe probably right?"
        table = fits.open('hizea_photo_galex_wise_v1.0.fit')
        print("table type", type(table))
        ngal = len(table[1].data)
        w3_mag = table[1].data['AB_MAG'][0:ngal][:, 9]
        print(type(w3_mag))
        w3_mag_err = table[1].data['AB_MAG_ERR'][0:ngal][:, 9]
        w4_mag = table[1].data['AB_MAG'][0:ngal][:, 10]
        w4_mag_err = table[1].data['AB_MAG_ERR'][0:ngal][:, 10]
        redshift = table[1].data['Z']
        names = table[1].data['SHORT_NAME'][0:ngal]

        w3minusw4_unc = np.sqrt(w3_mag_err ** 2 + w4_mag_err ** 2)
        valid = w3minusw4_unc < 0.5

        return proj_path,names,redshift, w3_mag,w3_mag_err,w4_mag,w4_mag_err,w3minusw4_unc,valid
        # (valid)print

# For 50 Galaxies
#I don't know what used to go here...

def main():
    pass
if __name__ == "main":
    main()

# This is an attempt at modularity. This sets the variables for the whole program, just pick the dataset
#I think it works!
projpath,ds_names, redshift, w3_mag, w3_mag_err, w4_mag, w4_mag_err, w3minusw4_unc, valid = choose_dataset("jwst")

name = 'Brown_SEDs_JWST.pdf'

with PdfPages(name) as pdf:
    print(len(files))
    letter = ord('c')
    best_fit = ['NGC_1275', "UGC_0869", "UGCA_219", "NGC_3690"]
    for i in range(0, len(files)):
        '''
        # want to read header in file, see if the units are wrong, and convert accordingly
        with open(files[i]) as f:
            header1 = f.readline()
            if 'Angstrom' in header1:
                pass
        '''
        table = ascii.read(files[i])
        # table = files[i]
        # table.pprint_all()

        wave = np.array(table['col1'])
        print("wave:", wave)
        flam = np.array(table['col2'])
        flux = wave * flam
        fnu = flam * wave ** 2 / 3e18  # * 1e46
        loc = files[0].find('spec')
        gal = files[i][loc - 9:loc - 1]
        # print(gal)

        filternames = ['wise_w{}'.format(n) for n in ['1', '2', '3', '4']]
        wise_filters = observate.load_filters(filternames)
        zarray = np.arange(41) / 100. + 0.4
        w3_minus_w4 = np.zeros(len(zarray))
        for j in range(0, len(zarray)):
            model_mags = observate.getSED(wave * (1 + zarray[j]), flam, filterlist=wise_filters)
            w3_minus_w4[j] = model_mags[2] - model_mags[3]

        top = np.max(w3_minus_w4)
        bot = np.min(w3_minus_w4)

        if (top > 1.2):
            print(gal)

            # 2nd plot -- fnu figure
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(211)

            plt.plot(wave / 1e4, fnu)
            # The commented out portion here was for quick coordination between plots and Google Sheets

            if ".txt" in files[0]:
                if gal in best_fit:
                    plt.title(f"({chr(letter)}) {gal}, good fit")
                else:
                    plt.title(f"({chr(letter)}) {gal}")


            letter += 1
            prange = wave > 3e4
            fmin = np.min(fnu[prange])
            fmax = np.max(fnu[prange])

            plt.yscale('log')
            plt.ylim([fmin / 2, fmax * 2])

            plt.xscale('log')
            plt.xlim(3, 40)

            #Why do we set the labels here too? Is this part of some other plot?
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel(r'Flux $f_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')

            ax.set_xticks([3, 5, 8, 10, 20, 30])
            ax.xaxis.set_major_formatter(ScalarFormatter())

            labels = [item.get_text() for item in ax.get_xticklabels()]
            labels[0] = '3'
            labels[1] = '5'
            labels[2] = '8'
            labels[3] = '10'
            labels[4] = '20'
            labels[5] = '30'

            ax.axvspan(12 / 1.4, 12 / 1.8, facecolor='g', alpha=0.5)
            ax.axvspan(22 / 1.4, 22 / 1.8, facecolor='r', alpha=0.5)

            ax = fig.add_subplot(212)

            plt.plot(zarray, w3_minus_w4)
            title = "Template Name: " + gal
            plt.title(title)

            plt.ylim([0, 3.5])
            plt.xlim([0.4, 0.8])

            plt.xlabel(r'Redshift')
            plt.ylabel(r'[W3] - [W4] color')

            plot_phot(redshift, ds_names)

            plt.tight_layout()

            pdf.savefig()
            plt.close()

os.system('open %s &' % name)
