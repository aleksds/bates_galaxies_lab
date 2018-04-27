#automating the creation and running of galfit input files
#in this stage, the code recreates J0826_F475W_coarse_withfinepsf except using a psf fit instead of a sersic fit

import os
import numpy as np

text = open('auto.txt','w')
text.write('#  Input menu file: J0826_F475W_coarse_withfinepsf\n') #propably not essential
text.write('#  Chi^2/nu = 1.845,  Chi^2 = 303988.812,  Ndof = 164804\n') #probably not essential
text.write('# IMAGE and GALFIT CONTROL PARAMETERS\n') #will not change , probably not essential
text.write('A) /usr/netapp/physics/linux-lab/data/hst/J0826+4305/coarse/F475W/final_F475W_drc_sci.fits\n')
text.write('B) J0826_F475W_coarse.fits\n')
text.write('C) none\n') #will not change
text.write('D) /usr/netapp/physics/linux-lab/data/hst/J0826+4305/fine/F475W/final_psf.fits\n')
text.write('E) 2\n') #will not change
text.write('F) none\n') #will not change
text.write('G) /usr/netapp/physics/linux-lab/data/galfit/sersic_constraint.txt\n')
text.write('H) 3428 3828 3954 4354\n')
text.write('I) 150    150\n') 
text.write('J) 26.563\n')
text.write('K) 0.050  0.050\n')
text.write('O) regular\n') #will not change
text.write('P) 0\n') #will not change

text.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
text.write('# Component number: 1\n')
text.write(' 0) psf\n')
text.write(' 1) 3629.7708 4154.9507 1 1  # position x, y        [pixel]\n')
text.write(' 3) 20.4776     1\n')
text.write(' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n')

text.close()


os.system('galfit %s' % text)
