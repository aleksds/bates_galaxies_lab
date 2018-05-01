#automating the creation and running of galfit input files
#trying to make it loop through all of 'withfinepsf' files

import os
import numpy as np
from astropy.io import ascii

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']
longgal = ['J0826+4305','J0901+0314','J0905+5759','J0944+0930','J1107+0417','J1219+0336','J1341+0321','J1506+5402','J1558+3957','J1613+2834','J2116-0634','J2140+1209']
photzeromag = ['25.613','25.027']
for w in range(0,12):
    for i in range(0,2):
        file = galaxies[w]+'_'+filters[i]+'_input.txt'
        text = open(file,'w')
        text.write('#  Input menu file: '+galaxies[w]+'_'+filters[i]+'_coarse_withfinepsf\n') #propably not essential
        text.write('#  Chi^2/nu = 1.845,  Chi^2 = 303988.812,  Ndof = 164804\n') #probably not essential
        text.write('# IMAGE and GALFIT CONTROL PARAMETERS\n') #will not change , probably not essential
        text.write('A) /Volumes/physics/linux-lab/data/hst/'+longgal[w]+'/coarse/'+filters[i]+'/final_'+filters[i]+'_drc_sci.fits\n')
        text.write('B) '+galaxies[w]+'_'+filters[i]+'_coarse.fits\n')
        text.write('C) none\n') #will not change
        text.write('D) /Volumes/physics/linux-lab/data/hst/'+longgal[w]+'/fine/'+filters[i]+'/final_psf.fits\n')
        text.write('E) 2\n') #will not change
        text.write('F) none\n') #will not change
        text.write('G) /Volumes/physics/linux-lab/data/galfit/sersic_constraint.txt\n')

        galcoords = 'galcoords.dat'
        catalog = ascii.read(galcoords)
        
        xcoor = str(catalog[w][1])
        ycoor = str(catalog[w][2])
        xcoorlow = str(catalog[w][1]-200)
        ycoorlow = str(catalog[w][2]-200)
        xcoorhigh = str(catalog[w][1]+200)
        ycoorhigh = str(catalog[w][2]+200)
        
        text.write('H) '+xcoorlow+' '+xcoorhigh+' '+ycoorlow+' '+ycoorhigh+'\n')
        text.write('I) 150    150\n') 
        text.write('J) '+photzeromag[i]+'\n')
        text.write('K) 0.050  0.050\n')
        text.write('O) regular\n') #will not change
        text.write('P) 0\n') #will not change

        text.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
        text.write('# Component number: 1\n')
        text.write(' 0) psf\n')
        text.write(' 1) '+xcoor+' '+ycoor+' 1 1  # position x, y        [pixel]\n')
        text.write(' 3) 20.4776     1\n')
        text.write(' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n')

        text.close()






#os.system('galfitm %s' % file)
