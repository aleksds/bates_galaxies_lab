#code to automatically make .txt files to be fed to prospector
#the plan is to set it up so this code takes an argument that specifies the folder with the galfit output files that it will draw info from
##sys.argv[1] is psf or sersic (model), sys.arg[2] is coarse or fine (plate scale), sys.arg[3] is a number (image region size), sys.arg[4] is coarse or fine (psf), sys.arg[5] is simultaneous or independent or semi
from astropy.io import ascii
import numpy as np
import sys
import os

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
redshifts = [0.603,0.459,0.712,0.514,0.467,0.451,0.451,0.658,0.608,0.402,0.449,0.728,0.752]
filters = ['F475W','F814W','F160W']
time = sys.argv[1]
model = sys.argv[2]
plate = sys.argv[3]
imgsize = sys.argv[4]
psf = sys.argv[5]
togetherness = sys.argv[6]

fluxes = np.zeros([12,3])
ivar = np.zeros([12,3])

if togetherness == 'independent':
    for w in range(0,12):
        for i in range(0,3):
            file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+galaxies[w]+'_'+filters[i]+'_'+model+'_'+plate+'_output.galfit.01.band'
            with open(file) as f:
                content = f.readlines()
            fluxes[w][i] = (10**((np.float(content[47][4:10]))/(-2.5))*10**9)
        
else: # togetherness == 'simultaneous' or 'semi':
    for w in range(0,12):
        file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+galaxies[w]+''+filters[1]+'_'+filters[0]+'_'+filters[2]+'_'+model+'_output.galfit.01.band'
        with open(file) as f:
            content = f.readlines()
        fluxes[w][0] = (10**((np.float(content[47][11:16]))/(-2.5))*10**9)
        fluxes[w][1] = (10**((np.float(content[47][4:9]))/(-2.5))*10**9)
        fluxes[w][2] = (10**((np.float(content[47][4:9]))/(-2.5))*10**9)
        
for w in range(0,12):
    for i in range(0,3):
        ivar[w][i] = 1/(fluxes[w][i]*0.05)**2

file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/pros_input.txt'
text = open(file,'w')
text.write('ID                f_475        ivar_475 f_814        ivar_814 f_160        ivar_160 z\n')
for w in range(0,12):
    text.write(galaxies[w]+'    '+str(round(fluxes[w][0],4))+'    '+str(round(ivar[w][0],4))+'    '+str(round(fluxes[w][1],4))+'    '+str(round(ivar[w][1],4))+'    '+str(round(fluxes[w][2],4))+'    '+str(round(ivar[w][2],4))+'    '+str(redshifts[w])+'\n')
text.close()
