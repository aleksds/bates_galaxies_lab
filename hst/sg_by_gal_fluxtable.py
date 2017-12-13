# Sophia C W Gottlieb I
# 20170612 produce table of flux values in nanomaggies in a txt file
# and graphs of the log of photons per pixel
#
# import relevant Python modules
import os, sys
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import math
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from scipy import integrate
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
from xlrd import open_workbook

rawconv = int(sys.argv[1])
coarsefine = int(sys.argv[2])
'''

rawconv = 1
coarsefine = 1
'''
# define the directory that contains the images
dir = os.environ['HSTDIR']
conv = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
class Galaxy:
    def __init__(self, name, z, x, y):
        self.name = name
        self.z = z
        self.x = x
        self.y = y
        self.lDcm = cosmo.luminosity_distance(self.z)*u.Mpc.to(u.cm) / u.Mpc
        self.radToKpc = conv.arcsec_per_kpc_proper(self.z)*0.05/u.arcsec*u.kpc
# We grab our galaxy data and make a list of galaxy objects
wbo = open_workbook('galaxydata.xlsx')
#sheet = wbo.sheets()[0]
for sheet in wbo.sheets():
    numr = sheet.nrows
    galaxies = []
    
    for row in range(1,numr):
        galaxies.append(Galaxy(sheet.cell(row,0).value,sheet.cell(row,1).value,sheet.cell(row,2).value,sheet.cell(row,3).value))

        
#setting up arrays with three elements, all zeros - placeholders
filters = ['F475W','F814W','F160W']
wavelengths = np.array([475, 814, 1600])
data = [0 for x in range(len(wavelengths))]
header = [0 for x in range(len(wavelengths))]
fnu = [0 for x in range(len(wavelengths))]
exp = [0 for x in range(len(wavelengths))]
gain = [0 for x in range(len(wavelengths))]

#width MUST be odd.
width = 15
pixr = int((width-1)/2)

# specify the position of the science target and the size of the region around the science target to consider

totalphotons = 0
fluxvalues = [[0 for x in range(len(wavelengths))] for y in range(len(galaxies))]
# define the radii to be used for aperture photometry
radii = np.arange(40)+1
area = [0 for x in range(len(radii))]

# percunc specifies the percent of variation we expect from systematic error... 
# For now, we have chosen 0.05, or 5%.
percUnc = 0.05

#calculate area of each bagel
for i in range(0, len(area)):
    if i == 0:
        area[i] = math.pi*math.pow(radii[0],2)
    else:
        area[i] = math.pi*(math.pow(radii[i],2)-math.pow(radii[i-1],2))

# Now, we loop through all galaxies
res = ['fine','coarse']
proc = ['convolved_image', 'final']
with PdfPages('sg_FLUX.pdf') as pdf: 
    
    for w in galaxies:
        
        print(w)
        fig = plt.figure()
        collection = ['F475W','F814W','F160W']
    
        flux = np.zeros([len(collection),len(radii)]) #*u.Jy
        subflux = np.zeros([len(collection),len(radii)])
        fluxpix = np.zeros([len(collection), width, width])
        for i in range (0, len(collection)):
            
            # read in the images
            #file = glob.glob(dir+w.name+'/'+res[coarsefine]+'/'+str(filters[i])+'/'+proc[rawconv]+'_'+str(filters[i])+'*sci.fits')

            file = glob.glob(dir+w.name+'/'+res[1]+'/'+str(filters[i])+'/'+proc[1]+'_'+str(filters[i])+'*sci.fits')

            #file = glob.glob(dir+galaxies[w]+'_final_'+collection[i]+'*sci.fits')
            hdu = fits.open(file[0])
            data[i], header[i] = hdu[0].data, hdu[0].header
            fnu[i] = header[i]['PHOTFNU']
            exp[i] = header[i]['EXPTIME']
            gain[i] = header[i]['CCDGAIN']
    
            #define positions for photometry
            positions = [(int(w.x), int(w.y))]
            totalphotons = 0
            '''
            # do pixel analysis
            for j in range(0,width):
                
                for k in range(0,width):
                    #print(i,j,k)
                    #fluxpix[i,width-j-1,k] = ((data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])*gain[i]*exp[i])
                    fluxpix[i,width-j-1,k] = math.log10((data[i][j+w.y-pixr+1][k+w.x-pixr+1])*gain[i]*exp[i])
                    totalphotons = totalphotons +  (data[i][j+w.y-pixr+1][k+w.x-pixr+1])*gain[i]*exp[i]
                    #converts units to Jy and then to nanomaggys: Jy is data * fnu / exp and 1 nm = 3.631e-6 Jy
                    #fluxnmaggys[i,width-j-1,k] = data[i][j+ycen[w]-pixr+1][k+xcen[w]-pixr+1])*fnu[i]/exp[i]/(3.631*10**-6)

            plt.imshow(fluxpix[i],cmap='gray')
            plt.colorbar()
            plt.title(w.name+' photon count: '+str(np.log10(totalphotons)) +' ('+ collection[i] + ')')
            plt.xlabel('Pixels')
            plt.ylabel('Pixels')
            pdf.savefig()
            plt.close()
            print(str(np.log10(totalphotons)))
            '''
# NANOMAGGYS
            for j in range(0,len(radii)):
                aperture = CircularAperture(positions, radii[j])
                phot_table = aperture_photometry(data[i], aperture)
                # CONVERT UNITS TO JY TO nMy
                flux[i,j] = phot_table['aperture_sum'][0]*(fnu[i]/exp[i])/(3.631*10**-6)
                #flux[i,j] = phot_table['aperture_sum'][0]*gain[i]*exp[i]
                if j == 0:
                    subflux[i,j] = flux[i,j]
                else:
                    subflux[i,j] = flux[i,j]-flux[i,j-1]

            fluxvalues[galaxies.index(w)]=subflux
print('WHY DONT YOU LOVE ME')

photfolder = 'PROSPECTOR/'+res[rawconv]+'_'+proc[coarsefine]+'/'
if not os.path.exists(photfolder):
    os.makedirs(photfolder)
    
# we shift through the data by galaxy and then by aperture.
for w in range(0,len(galaxies)):
    gal = galaxies[w].name[:5]
    # we create or open our txt file
    f = open(photfolder+ gal+'.txt',"w+")
    print(photfolder+ gal+'.txt')

    # we write our column titles - unsure if these need to stay, just thought it would be nice.
    f.write('ID\t\tf_475\t\tivar_475\t\tf_814\t\tivar_814\t\tf_1600\t\tivar_1600\t\tz\n')
    
    for i in range(0,len(radii)):
        # Building the ID name using if/else for those with single digits.
        
        ID = gal+'_'
        if i < 10:
            ID = ID + '0' + str(i)
        else:
            ID = ID + str(i)
        # writing data divided by tabs (\t character)
        f.write(ID+'\t')
        # we loop through the filters because i am lazy and it is technically good form.
        # at the same time, we grab the inverse variance (squared) which is 1/(flux*res)^2.
        for j in range(0,len(collection)):
            f.write(str(fluxvalues[w][j][i])+'\t')
            ivar = (fluxvalues[w][j][i]*percUnc)**(-2)
            f.write(str(ivar)+'\t')
        # to conclude the line, we include the z value for the galaxy before calling for a new line (\n).
        f.write(str(galaxies[w].z)+'\n')

    f.close()
print('YOURE NOT MY REAL MOM')
