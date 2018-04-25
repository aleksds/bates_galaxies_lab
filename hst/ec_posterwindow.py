import numpy as np
from astropy.io import fits
import glob
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM 
import img_scale
from xlrd import open_workbook

conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

class Galaxy:
    def __init__(self, name, z, x, y):
        self.name = name
        self.z = z
        self.x = x
        self.y = y
        self.lDcm = cosmo.luminosity_distance(self.z)*u.Mpc.to(u.cm) / u.Mpc
        self.radToKpc = conv.arcsec_per_kpc_proper(self.z)*0.05/u.arcsec*u.kpc

# Grabbing and filling galaxy data
cf = ['coarse','fine']

# label describing image processing
rc = ['final','convolved_image']
# matching file name ending whoooooooot.
# not confident the h in needed in whooooot. but lets roll with it for now. its a nice h
fex = ['*sci.fits','.fits']

galaxies = []

wbo = open_workbook('galaxydata.xlsx')
for sheet in wbo.sheets():
    numr = sheet.nrows
    
    for row in range(1,numr):
        galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))



# define the directory that contains the images
dir = os.environ['HSTDIR']
# define a function to plot "postage stamp" images
def plot_image():
    rotated = np.flip(np.rot90(stamp, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
    #plt.tick_params(axis='both', which='major', labelsize=8)
    
filters = ['F475W','F814W','F160W']
res = ['fine','coarse']
type = ['convolved_image', 'final']

with PdfPages('ec_poster.pdf') as pdf:
    for g, gal in enumerate(galaxies):
        fig = plt.figure()
        plt.suptitle(gal.name)
        alldata = []
        for f, fil in enumerate(filters):
            file = glob.glob(dir+gal.name+'*/'+res[1]+'/'+fil+'/'+type[1]+'_'+fil+'*sci.fits')

            hdu = fits.open(file[0])
            data, header = hdu[0].data, hdu[0].header
            width = 100
            dy = width
            dx = width

            stamp = data[round(gal.y-dy):round(gal.y+dy), round(gal.x-dx):round(gal.x+dx)]
            ax = fig.add_subplot(1,len(filters)+1, 1+f)
            #ax = fig.add_subplot(1,4, 1+f)
            plt.axis('off')
            plot_image()
            plt.title(fil)

            if f==0:
                image_b = np.flip(np.rot90(stamp, 2), 1)
            if f==1:
                image_g = np.flip(np.rot90(stamp, 2), 1)
            if f==2:
                image_r = np.flip(np.rot90(stamp, 2), 1)
            
            alldata.append(stamp)

        bx = fig.add_subplot(1, len(filters)+1, 4)
        #bx = fig.add_subplot(1, 4, 4)
        img = np.zeros([width*2,width*2,3])
        img[:,:,0] = img_scale.log(image_r, scale_min=5, scale_max=10000)
        img[:,:,1] = img_scale.log(image_g, scale_min=5, scale_max=10000)
        img[:,:,2] = img_scale.log(image_b, scale_min=5, scale_max=10000)
        plt.imshow(img)
        plt.axis('off')
        plt.title('color')
        pdf.savefig(dpi=1000)
        plt.close()

os.system('open %s &' % 'ec_poster.pdf')
                             
