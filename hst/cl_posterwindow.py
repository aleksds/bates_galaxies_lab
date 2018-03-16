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
def plot_image_1():
    rotated = np.flip(np.rot90(stampdata, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
def plot_image_2():
    rotated = np.flip(np.rot90(stampmodel, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
def plot_image_3():
    rotated = np.flip(np.rot90(stampres, 2), 1)
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='none')
    
filters = ['_F475W','_F814W']
type = ['data','model','residual']
folders = ['475coarsefinepsf/', '814coarsefinepsf/']
dx = dy = 60
with PdfPages('cl_poster.pdf') as pdf:
    for g in range(0,len(galaxies)):
        fig = plt.figure()
        plt.suptitle(galaxies[g])
        alldata = []
        for f in range(0,len(filters)):
            for t in range(0,len(type)):
                file = glob.glob('/Volumes/physics/linux-lab/data/galfit/'+folders[f]+galaxies[g]+filters[f]+'_coarse.fits') 
                                 
                multi = fits.open(file[0])
                
                data, data_header = multi[1].data, multi[1].header
                model, res_header = multi[2].data, multi[2].header
                res, res_header = multi[3].data, multi[3].header
 

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
              
                
#                if f==0:
 #                   image_data_f = np.flip(np.rot90(stampdata, 2), 1)
  #                  image_model_f = np.flip(np.rot90(stampmodel, 2), 1)
   #                 image_res_f = np.flip(np.rot90(stampres, 2), 1)
    #            if f==1:
     #               image_data_e = np.flip(np.rot90(stampdata, 2), 1)
      #              image_model_e = np.flip(np.rot90(stampmodel, 2), 1)
       #             image_res_e = np.flip(np.rot90(stampres, 2), 1)
            
        #    alldata.append(stampdata)
            
        pdf.savefig(dpi=1000)
        plt.close()

os.system('open %s &' % 'cl_poster.pdf')



        
        #img = np.zeros([6,120,120])
        #img[0,:,:] = img_scale.log(image_data_f, scale_min=5, scale_max=10000) #HEYCHARLESHERESTHEPROBLEMDOOOOODINATOR
       # img[1,:,:] = img_scale.log(image_model_f, scale_min=5, scale_max=10000)
      #  img[2,:,:] = img_scale.log(image_res_f, scale_min=5, scale_max=10000)
     #   img[3,:,:] = img_scale.log(image_data_e, scale_min=5, scale_max=10000)
    #    img[4,:,:] = img_scale.log(image_model_e, scale_min=5, scale_max=10000)
   #     img[5,:,:] = img_scale.log(image_res_e, scale_min=5, scale_max=10000)
  #      plt.imshow(img[0,:,:])
 #       plt.axis('off')
#        plt.title('color')
        


                             
