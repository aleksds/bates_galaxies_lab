#Following the code in ad_posterwindow.py or sg_posterwindow.py, create a new piece of code that reads in the GALFIT output FITS files and shows the data, model, and residual for F475W and F814W for each galaxy (e.g, showing 6 images on on a single page for one galaxy, repeating 12 times)



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

#class Galaxy:
 #   def __init__(self, name, z, x, y):
  #      self.name = name
   #     self.z = z
    #    self.x = x
     #   self.y = y
      #  self.lDcm = cosmo.luminosity_distance(self.z)*u.Mpc.to(u.cm) / u.Mpc
       # self.radToKpc = conv.arcsec_per_kpc_proper(self.z)*0.05/u.arcsec*u.kpc

# Grabbing and filling galaxy data
cf = ['coarse','fine']

# matching file name ending whoooooooot.
# not confident the h in needed in whooooot. but lets roll with it for now. its a nice h
fex = ['*sci.fits','.fits']

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']

#wbo = open_workbook('galaxydata.xlsx')
#for sheet in wbo.sheets():
 #   numr = sheet.nrows
    
   # for row in range(1,numr):
    #    galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))



# define the directory that contains the images
dir = os.environ['HSTDIR']
# define a function to plot "postage stamp" images
def plot_image():
    rotated = np.flip(np.rot90(stamp, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
    #plt.tick_params(axis='both', which='major', labelsize=8)
    
filters = ['_F475W','_F814W']
type = ['data','model','residual']
xcen = [3629, 3934, 3387, 3477, 3573, 3803, 3884, 4147, 3787, 4175, 3567, 4067]
ycen = [4154, 4137, 3503, 3405, 3339, 4170, 4165, 3922, 4186, 3827, 3436, 4054]
folders = ['475coarsefinepsf/', '814coarsefinepsf/']

with PdfPages('cl_poster.pdf') as pdf:
    for g in range(0,len(galaxies)):
        fig = plt.figure()
        plt.suptitle("YOOO")
        alldata = []
        for f in range(0,len(filters)):
            for t in range(0,len(type)):
                file = glob.glob('/Volumes/physics/linux-lab/data/galfit/'+folders[f]+galaxies[g]+filters[f]+'_coarse.fits') 
                                 
                multi = fits.open(file[0])
                
                data, data_header = multi[1].data, multi[1].header
                res, res_header = multi[2].data, multi[2].header
                res, res_header = multi[3].data, multi[3].header
                #define positions for photometry
                positions = [(xcen[g], ycen[g])]

                stamp = data[round(ycen[g]-dy):round(ycen[g]+dy), round(xcen[g]-dx):round(xcen[g]+dx)]
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
        img = np.zeros([120,120,3])

        img[:,:,0] = img_scale.log(image_r, scale_min=5, scale_max=10000)
        img[:,:,1] = img_scale.log(image_g, scale_min=5, scale_max=10000)
        img[:,:,2] = img_scale.log(image_b, scale_min=5, scale_max=10000)
        plt.imshow(img)
        plt.axis('off')
        plt.title('color')
        pdf.savefig(dpi=1000)
        plt.close()

os.system('open %s &' % 'ad_poster.pdf')
                             
