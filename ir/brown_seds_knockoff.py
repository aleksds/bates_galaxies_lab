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

### For Brown
projpath = os.getcwd()+'/Brown2014/An_Atlas_of_Galaxy_SEDs/An_Atlas_of_Galaxy_SEDs/'


###For our galaxies from Petter
#projpath = os.getcwd() + '/'
# print(projpath)
files = glob.glob(projpath + '*.dat')
#print("len file", len(files))
#print(type(files))

# For 50 Galaxies

#table = fits.open('hizea_photo_galex_wise_v1.0.fit')
#print("table type", type(table))
'''
ngal = len(table[1].data)
w3_mag = table[1].data['AB_MAG'][0:ngal][:,9]
print(type(w3_mag))
w3_mag_err = table[1].data['AB_MAG_ERR'][0:ngal][:,9]
w4_mag = table[1].data['AB_MAG'][0:ngal][:,10]
w4_mag_err = table[1].data['AB_MAG_ERR'][0:ngal][:,10]
redshift = table[1].data['Z']
names = table[1].data['SHORT_NAME'][0:ngal]

w3minusw4_unc = np.sqrt(w3_mag_err ** 2 + w4_mag_err ** 2)
valid = w3minusw4_unc < 0.5
# (valid)print
'''
# For 8 Galaxies of Interest
# Need to rewrite code for use in block format, not loop format, right now it does WISEdir for every galaxy,
# want a table w/ all of them
'''
projpath = os.getcwd() + '/'
redshift = [0.603, 0.711, 0.514, 0.467, 0.451, 0.661, 0.449, 0.459]
special_names = ['J0826', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1613', 'J2118']
w3_mag = np.zeros(len(redshift))
w3_mag_err = np.zeros(len(redshift))
w4_mag = np.zeros(len(redshift))
w4_mag_err = np.zeros(len(redshift))
table = []
count = 0
'''

####################################
'''For 6 Galaxies accepted to JWST proposal. Names are listed and stuff?
The solution is inelegant but we use the same w3_mag etc and just comment out the top as necessary
'''

projpath = os.getcwd() + '/'
redshift = [0.467, 0.451, 0.437, 0.449, 0.459]
jwst_names = ['J1107', 'J1219', 'J1506', 'J1613', 'J2118']
w3_mag = np.zeros(len(redshift))
w3_mag_err = np.zeros(len(redshift))
w4_mag = np.zeros(len(redshift))
w4_mag_err = np.zeros(len(redshift))
table = []
count = 0

def to_ab(vega_mag, band):

    '''Converts from vega magnitudes to ab mags '''

    # bandpasses = {"W1":2.699,"W2":3.339,"W3":5.147,"W4":6.620}
    zero_mag_fnu = {"W1": 309.540, "W2": 171.787, "W3": 31.674, "W4": 8.363}
    f_nu = float(zero_mag_fnu[band]*(10**(-vega_mag/2.5)))
    # print("type f_nu", type(f_nu))
    # print('f_nu',f_nu)
    mag_ab = float(-2.5 *np.log10(f_nu/3631))
    # print("type mag_ab", type(mag_ab))
    # ('mag_ab',mag_ab)

    return mag_ab

#Converting to ab magnitude for each bandpass  data of each galaxy
#print(float(to_ab(4.852,'W3')))

#for s in special_names:
for j in jwst_names:
    WISEdir = projpath + 'unWISE/%s.fits' % j
    t = Table.read(WISEdir)
    #t.pprint_all()
    table.append(t)
    #print('t type', type(t))
    w3_mag[count] = float(to_ab(t['w3_mag'],"W3"))
    #("t['w3_mag_err']")
    #print(t['w3_mag_err'])
    w3_mag_err[count] = float(t['w3_mag_err'])
    #print("w3_mag_err_ab")
    #print(w3_mag_err[count])
    w4_mag[count] = float(to_ab(t['w4_mag'],"W4"))
    #print("w4_mag_err")
    #(t['w4_mag_err'])
    w4_mag_err[count] = float(t['w4_mag_err'])
    # print("w4_mag_err_ab")
    # print(w4_mag_err[count])

    count += 1

# Can we add all the information into one file?
# Go through each file, store the appropriate information in separate lists?
# Then put it back together, similar to RGB ppm\

"""
redshift = np.asarray(redshift)
special_names = np.asarray(special_names)
w3_mag = np.asarray(w3_mag)
w3_mag_err = np.asarray(w3_mag_err)
w4_mag= np.asarray(w4_mag)
w4_mag_err = np.asarray(w4_mag_err)
w3minusw4_unc = np.sqrt(w3_mag_err ** 2 + w4_mag_err ** 2)
print('w3minusw4_unc')
print(w3minusw4_unc)
valid = w3minusw4_unc < 0.5

print(valid)
"""



#This is the same as above but for jwst names instead. Again this is terrible and messy  but it works
redshift = np.asarray(redshift)
jwst_names = np.asarray(jwst_names)
w3_mag = np.asarray(w3_mag)
w3_mag_err = np.asarray(w3_mag_err)
w4_mag= np.asarray(w4_mag)
w4_mag_err = np.asarray(w4_mag_err)
w3minusw4_unc = np.sqrt(w3_mag_err ** 2 + w4_mag_err ** 2)
#print('w3minusw4_unc')
#print(w3minusw4_unc)
valid = w3minusw4_unc < 0.5

print(valid)


def plot_phot(redshift):
    print("Valid:",valid)
    print('w3_mag[valid]',w3_mag[valid])
    print ('w4_mag[valid]', w4_mag[valid])
    #print('redshift',redshift)
    # print("redshift valid",type(redshift[valid]))
    ax.scatter(redshift[valid], w3_mag[valid] - w4_mag[valid], marker='*')
    ax.errorbar(redshift[valid], w3_mag[valid] - w4_mag[valid], yerr=w3minusw4_unc[valid], color='purple',
                linestyle="None")

    ## For Brown the difference is names v special_names
    #This is for all valid names
    # for i in range(0,len(redshift[valid])):
    # plt.text(redshift[valid][i], w3_mag[valid][i]-w4_mag[valid][i], names[valid][i], fontsize=6)

    ##This is for the 8 Galaxies of interest
    '''
    for i in range(0, len(redshift[valid])):
        #print('w3_mag', w3_mag)
       # print('w3_mag[valid][0]', w3_mag[valid][0])
        plt.text(redshift[valid][i], w3_mag[valid][i] - w4_mag[valid][i], special_names[valid][i], fontsize=6)
    '''
    #This is for the 5 JWST Galaxies


    for i in range(0, len(redshift[valid])):
        # print('w3_mag', w3_mag)
        # print('w3_mag[valid][0]', w3_mag[valid][0])
        plt.text(redshift[valid][i], w3_mag[valid][i] - w4_mag[valid][i], jwst_names[valid][i], fontsize=6)


name = 'brown_seds_jswt_gals.pdf'

with PdfPages(name) as pdf:
    #print("and we got here folks")
    print(len(files))
    letter = ord('c')
    best_fit = ['NGC_1275', "UGC_0869","UGCA_219","NGC_3690"]
    for i in range(0, len(files)):
        # for i in range(0,1):
       # print('we\'re in')
        #print('files[i] type',type(files[i]))
        table = ascii.read(files[i])
        #table = files[i]
       # print("It's here")
        #print('table type',type(table))
        #print(table)
        #table.pprint_all()

        wave = np.array(table['col1'])

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
           # print('we\'re in deeper')
            model_mags = observate.getSED(wave * (1 + zarray[j]), flam, filterlist=wise_filters)
            w3_minus_w4[j] = model_mags[2] - model_mags[3]
            # print('model_mags:', model_mags)
            # print('[W3]-[W4]:', w3_minus_w4[j])

        top = np.max(w3_minus_w4)
        bot = np.min(w3_minus_w4)

        if (top > 1.2):
            print(gal)

            ## lambda * f_lambda figure
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            #
            # plt.plot(wave/1e4, flux)
            #
            # plt.title(gal)
            #
            # prange = wave > 3e4
            # fmin = np.min(flux[prange])
            # fmax = np.max(flux[prange])
            #
            # plt.yscale('log')
            # plt.ylim([fmin/2,fmax*2])
            #
            # plt.xscale('log')
            # plt.xlim(3,40)
            #
            # plt.xlabel(r'Wavelength [$\mu$m]')
            # plt.ylabel(r'Flux $\lambda f_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$]')
            #
            # ax.set_xticks([3, 5, 8, 10, 20, 30])
            # ax.xaxis.set_major_formatter(ScalarFormatter())
            #
            # labels = [item.get_text() for item in ax.get_xticklabels()]
            # labels[0] = '3'
            # labels[1] = '5'
            # labels[2] = '8'
            # labels[3] = '10'
            # labels[4] = '20'
            # labels[5] = '30'
            #
            # pdf.savefig()
            # plt.close()

            # 2nd plot -- fnu figure
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot(211)

            plt.plot(wave / 1e4, fnu)
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


            ## 3rd plot -- flambda figure
            # fig = plt.figure()
            #
            # ax = fig.add_subplot(111)
            #
            # plt.plot(wave/1e4, flam)
            #
            # plt.title(gal)
            #
            # prange = wave > 3e4
            # fmin = np.min(flam[prange])
            # fmax = np.max(flam[prange])
            #
            # plt.yscale('log')
            # plt.ylim([fmin/2,fmax*2])
            #
            # plt.xscale('log')
            # plt.xlim(3,40)
            #
            # plt.xlabel(r'Wavelength [$\mu$m]')
            # plt.ylabel(r'Flux $f_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')
            #
            # ax.set_xticks([3, 5, 8, 10, 20, 30])
            # ax.xaxis.set_major_formatter(ScalarFormatter())
            #
            # labels = [item.get_text() for item in ax.get_xticklabels()]
            # labels[0] = '3'
            # labels[1] = '5'
            # labels[2] = '8'
            # labels[3] = '10'
            # labels[4] = '20'
            # labels[5] = '30'
            #
            # plt.tight_layout()
            #
            # pdf.savefig()
            # plt.close()

            # plot of WISE W3-W4 clor
            # fig = plt.figure()
            ax = fig.add_subplot(212)

            plt.plot(zarray, w3_minus_w4)

            # plt.title(gal)

            plt.ylim([0, 3.5])
            plt.xlim([0.4, 0.8])

            plt.xlabel(r'Redshift')
            plt.ylabel(r'[W3] - [W4] color')

            plot_phot(redshift)

            plt.tight_layout()

            pdf.savefig()
            plt.close()

os.system('open %s &' % name)
