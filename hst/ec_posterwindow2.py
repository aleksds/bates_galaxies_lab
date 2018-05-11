# This code takes cl's version and alters it so that it reads in galfit data that I produced for the J0905 galaxy
# takes data from both the individual and simultaneous runs

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


#wbo = open_workbook('galaxydata.xlsx')
#for sheet in wbo.sheets():
 #   numr = sheet.nrows
    
   # for row in range(1,numr):
    #    galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))



# define the directory that contains the images
dir = os.environ['HSTDIR']
# define a function to plot "postage stamp" images
def plot_image_1():
    rotated = np.flip(np.rot90(stampdata, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
def plot_image_2():
    rotated = np.flip(np.rot90(stampmodel, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
def plot_image_3():
    rotated = np.flip(np.rot90(stampres, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
    

filename = ['F475W','F814W','F814W_F475W']
type = ['data','model','residual']
dx = dy = 100
with PdfPages('ec_poster.pdf') as pdf:
    fig = plt.figure()
    plt.suptitle('J0905 Individual Model')
    alldata = []
    for f in range(0,2):
        for t in range(0,len(type)):
            file = glob.glob('/Volumes/physics/linux-lab/data/galfit/eves_files/ec_J0905_'+filename[f]+'_fine_sersic_output.fits') 
                             
            multi = fits.open(file[0])
            print(file)
            
            data, data_header = multi[1].data, multi[1].header
            model, res_header = multi[2].data, multi[2].header
            res, res_header = multi[3].data, multi[3].header


            stampdata = data[round(401-dy):round(401+dy), round(401-dx):round(401+dx)] 
            stampmodel = model[round(401-dy):round(401+dy), round(401-dx):round(401+dx)]
            stampres = res[round(401-dy):round(401+dy), round(401-dx):round(401+dx)]
            
            if f ==0:
                
                ax = fig.add_subplot(2,3,1)
                plt.axis('off')
                plt.title('F475W Data')
                plot_image_1()
                ax = fig.add_subplot(2,3,2)
                plt.axis('off')
                plt.title('F475W Model')
                plot_image_2()
                ax = fig.add_subplot(2,3,3)
                plt.axis('off')
                plt.title('F475W Residual')
                plot_image_3()

            if f ==1: 
                ax = fig.add_subplot(2,3,4)
                plt.axis('off')
                plt.title('F814W Data')
                plot_image_1()
                ax = fig.add_subplot(2,3,5)
                plt.axis('off')
                plt.title('F814W Model')
                plot_image_2()
                ax = fig.add_subplot(2,3,6)
                plt.axis('off')
                plt.title('F814W Residual')
                plot_image_3()
          
                
    pdf.savefig(dpi=1000)
    plt.close()

    dx = dy = 200
    fig = plt.figure()
    plt.suptitle('J0905 Simultaneous Model')
    alldata = []
    for f in range(0,2):
        for t in range(0,len(type)):
            file = glob.glob('/Volumes/physics/linux-lab/data/galfit/eves_files/sersic/ec_J0905_'+filename[2]+'_sersic_output.fits') 
                             
            multi = fits.open(file[0])
            
            data, data_header = multi[f].data, multi[f].header
            model, res_header = multi[2+f].data, multi[2+f].header
            res, res_header = multi[4+f].data, multi[4+f].header


            stampdata = data[round(201-dy):round(201+dy), round(201-dx):round(201+dx)] 
            stampmodel = model[round(201-dy):round(201+dy), round(201-dx):round(201+dx)]
            stampres = res[round(201-dy):round(201+dy), round(201-dx):round(201+dx)]
            
            if f ==0:
                
                ax = fig.add_subplot(2,3,1)
                plt.axis('off')
                plt.title('F475W Data')
                plot_image_1()
                ax = fig.add_subplot(2,3,2)
                plt.axis('off')
                plt.title('F475W Model')
                plot_image_2()
                ax = fig.add_subplot(2,3,3)
                plt.axis('off')
                plt.title('F475W Residual')
                plot_image_3()

            if f ==1: 
                ax = fig.add_subplot(2,3,4)
                plt.axis('off')
                plt.title('F814W Data')
                plot_image_1()
                ax = fig.add_subplot(2,3,5)
                plt.axis('off')
                plt.title('F814W Model')
                plot_image_2()
                ax = fig.add_subplot(2,3,6)
                plt.axis('off')
                plt.title('F814W Residual')
                plot_image_3()


    pdf.savefig(dpi=1000)
    plt.close()

os.system('open %s &' % 'ec_poster.pdf')

#notes on the images produced: F475W: individual simulation produces a large amount of negative pixels in the centroid, while the simulataneous fit  has a much higher percentage of positve pixels
