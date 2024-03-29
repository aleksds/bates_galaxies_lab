#Cristopher Thompson
#Last Edit: 05/07/2019


import os
import numpy as np
from astropy.io import ascii
import datetime
import sys
import glob

#Go back and make sure the table is being read in the same way
jc_values = ascii.read('ad_mag_size_table.dat')
#re_vals = np.zeros((7, len(jc_values)))
#re_vals[0] = jc_values['re00']
#re_vals[1] = jc_values['re01']
#re_vals[2] = jc_values['re_small']
#re_vals[3] = jc_values['re']
#re_vals[4] = jc_values['re_large']
#re_vals[5] = jc_values['re05']
#re_vals[6] = jc_values['re06']
#print(re_vals)
ngal = len(jc_values)
#nre = len(re_vals)
# new and deletable: Everything in white can be deleted after we figure out
#the values and everything in comment can be brought back out in color.
re_array = np.array([1/10, 1/5, 1/3, 1/2, 1/1.5, 1/1.2, 1, 1.2, 1.5, 2, 3, 5, 10])
nre = len(re_array)
#tmp = np.zeros((ngal, len(re_vals)))
#for i in range(0,ngal):
#    for j in range(0,len(re_vals)):
#        tmp[i,j] = re_vals[j,i]
#print(tmp)

tmp = np.zeros((ngal, nre))
for i in range(0,ngal):
    for j in range(0,nre):
        tmp[i,j] = jc_values['re'][i] * re_array[j]
print(tmp)

mag_values = ascii.read('ad_mag_size_table.dat')
mag_array = np.array([-0.5, -0.3, -0.2, -0.15, - 0.1, - 0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5])
#mag_array = np.array([-0.5, -0.3, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5])
#mag_array = np.array([-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2,\
#0.3, 0.4, 0.5, 0.6, 0.7])
nmag = len(mag_array)
#'md1f' stands for "magnitude difference of 1st filter"(i.e. 475W). 
md1f = np.zeros((ngal, nmag))
#'md2f' stands for "magnitude difference of 2nd filter"(i,e, 814W).
md2f = np.zeros((ngal, nmag))
for i in range(0, ngal):
    for j in range(0, nmag):
        md1f[i,j] = mag_values['m475'][i] + mag_array[j]
        md2f[i,j] = mag_values['m814'][i] + mag_array[j]
print(md1f)
print(md2f)

galaxies = [
    'J0826',
    'J0901',
    'J0905',
    'J0944',
    'J1107',
    'J1219',
    'J1341',
    'J1506',
    'J1558',
    'J1613',
    'J2116',
    'J2140']

filters2 = ['F475W','F814W']
filters = ['F475W','F814W','F160W']

longgal = [
    'J0826+4305',
    'J0901+0314',
    'J0905+5759',
    'J0944+0930',
    'J1107+0417',
    'J1219+0336',
    'J1341-0321',
    'J1506+5402',
    'J1558+3957',
    'J1613+2834',
    'J2116-0634',
    'J2140+1209']

photzeromag = ['25.613',
               '25.027']

kind = ['independent',
        'simultaneous']

""" sys.argv[1] is psf or sersic (model), sys.arg[2] is coarse or fine (plate
scale), sys.arg[3] is a number (image region size), sys.arg[4] is coarse or
fine (psf), sys.arg[5] is simultaneous or independent or semi, sys.arg[6] is
either none or s_index or re or s_index, re """

model = sys.argv[1]
plate = sys.argv[2]
imgsize = sys.argv[3]
psf = sys.argv[4]
togetherness = sys.argv[5]
constraint = sys.argv[6]

now = datetime.datetime.now()
time = now.strftime("%Y%m%d-%H%M")

hstdir = os.environ['HSTDIR']

#INDEPENDENT FITS
if togetherness == 'independent':
    for w in range(0,ngal):
        galcoords = 'galcoords_'+plate+'.dat'
        catalog = ascii.read(galcoords)

        xcoor = str(catalog[w][1])
        ycoor = str(catalog[w][2])
        xcoorlow = str(int(catalog[w][1]-(np.float(imgsize)/2)))
        ycoorlow = str(int(catalog[w][2]-(np.float(imgsize)/2)))
        xcoorhigh = str(int(catalog[w][1]+(np.float(imgsize)/2)))
        ycoorhigh = str(int(catalog[w][2]+(np.float(imgsize)/2)))
        for i in range(0,len(filters2)):
            sigma_file = glob.glob(hstdir+'sigma/'+longgal[w]+'_'+filters2[i]+'_'+plate+'_sigma.fits')
            print(sigma_file)
            for j in range(0,nmag):
                #if j == 0:
                #    md1f[w,0] = mag_values['m475'][w]
                #    md2f[w,0] = mag_values['m814'][w]
                #else:
                #    md1f[w,j] = mag_values['m475'][w] + re_array[j-1]
                #    md2f[w,j] = mag_values['m814'][w] + re_array[j-1]
                for q in range(0,nre):
                    #if q == 0:
                    #    tmp[w,0] = jc_values['re'][w]
                    #else:
                    #    tmp[w,q] = jc_values['re'][w] + re_array[q-1]
                    if not os.path.exists(time+'_'+model+'_'+plate+'_' \
                    +togetherness+ '_'+imgsize+'_psf'+psf):
                        os.makedirs(time+'_'+model+'_'+plate+'_'\
                        +togetherness+'_' +imgsize+'_psf'+psf)

                    file = time+'_'+model+'_'+plate+'_'+togetherness+'_' \
                        +imgsize+'_psf'+psf+'/'+galaxies[w]+'_'+filters2[i]+ \
                        '_'+model+'_effective_re'+str(q)+'_magnitude'+str(j)+\
                        '_input.txt'
                    text = open(file,'w')
                    text.write('#  Input menu file: '+galaxies[w]+'_' \
                    +filters2[i]+'_'
                    +plate+'\n') #probably not essential
                    text.write('#  Chi^2/nu = ,  Chi^2 = ,  Ndof = \n')
                    #probably not essential
                    text.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
                    #will not change , probably not essential
                    text.write('A) '+hstdir+'/' \
                    +longgal[w]+'/'+plate+'/'+filters2[i]+'/final_' \
                    +filters2[i]+'_drc_sci.fits\n')
                    text.write('B) '+galaxies[w]+'_'+filters2[i]+\
                    '_'+model+'_effective_re'+str(q)+'_magnitude'+str(j)+\
                    '_output.fits\n')
                    text.write('C) '+sigma_file[0]+' \n') #will not change
                    text.write('D) '+hstdir+'/'  \
                    +longgal[w]+'/'+plate+'/'+filters2[i]+'/final_psf.fits\n')
                    text.write('E) 1\n') #will not change
                    text.write('F) none\n') #will not change
                    #if constraint == 'none':
                    #    text.write('G) \n')
                    #if constraint == 's_index':
                    #    text.write('G) /Users/adiamond/physics/linux-lab/data/sersic\
                    #    _index_equaltofour.txt\n')
                    #if constraint == 're':
                    #     text.write('G) /Users/adiamond/physics/linux-lab/data/size_\
                    #     constraint_1to1.txt\n')
                    #if constraint == 's_index, re':
                    #    text.write('G) /Users/adiamond/physics/linux-lab/data/s_\
                    #    index_and_size.txt\n')

                    if model == 'sersic':
                        text.write('G) \n')
                        constraintfile= time+'_'+model+'_'+plate+'_'\
                        +togetherness+'_'+imgsize+'_psf'+psf+'/'+'parameters_constraints_'\
                        +filters2[i]+'_'+str([j])+'_'+str([q])+\
                        '_'+galaxies[w]+'.txt'
                    
                        contraint_text = open(constraintfile, 'w')
                        contraint_text.write('# Component/    parameter\
                        constraint	Comment\n')
                        contraint_text.write('# operation	(see below)\
                        range\n\n')
                        contraint_text.write('1              n             \
                        4 to 4   # Soft constraint: Constrains the\n')
                        contraint_text.write('					# sersic\
                        index n to within values\n')
                        contraint_text.write('				        # from\
                        0.7 to 5.\n')
                        contraint_text.write('1              re         '\
                        +str(tmp[w][q])+' to '+str(tmp[w][q])+' \n')
                        if i==0:
                            contraint_text.write('1              mag     \
                            '+str(md1f[w][j])+' to '+str(md1f[w][j])+' \n')
                        if i==1:
                            contraint_text.write('1              mag     \
                            '+str(md2f[w][j])+' to '+str(md2f[w][j])+' \n')
                        
                        contraint_text.write('1              q            \
                        1 to 1  \n')
                        contraint_text.write('1              pa           \
                        0 to 0  \n')
                        contraint_text.close()
                        text.write('G) parameters_constraints_'+filters2[i]+\
                        '_'+str([j])+'_'+str([q])+'_'+galaxies[w]+'.txt\n')
                        


                        text.write('H) '+xcoorlow+' '+xcoorhigh+' '+ycoorlow+\
                        ' \
                        '+ycoorhigh+'\n')
                        text.write('I) 150    150\n')
                        text.write('J) '+photzeromag[i]+'\n')
                        
                        if plate == 'coarse':
                            text.write('K) 0.050  0.050\n')
                        if plate == 'fine':
                            text.write('K) 0.025  0.025\n')
                        
                        text.write('O) regular\n') #will not change
                        text.write('P) 0\n') #will not change
                        text.write('W) input,sigma,psf,component,model,residual\n')
                        text.write('#   par)    par value(s)    fit\
                        toggle(s)   # parameter description \n')
                        text.write('# Component number: 1\n') 
                        
                        if model == 'psf':
                            text.write(' 0) psf\n')
                            text.write(' 1) '+xcoor+' '+ycoor+'\
                             1 1  # position x, y [pixel]\n')
                            if i==0:
                                mag = mag_values [w][1] #the m475 value
                                text.write(' 3) '+str(md1f[w][j])+'     1/n')
                            if i==1:
                                mag = mag_values[w][2] #the m814 values
                                text.write(' 3) '+str(md2f[w][j])+'     1\n')
                                #change it to where it is pulling from Aleks'
                                #data table.Just one by itself
                            text.write(' Z) 0 #Skip this model in output\
                            image? \
                            (yes=1, no=0)\n')
                        if model == 'sersic':
                            text.write(' 0) sersic\n')
                            text.write(' 1) '+xcoor+' '+ycoor+\
                            ' 1 1  # position x, y [pixel]\n')
                            if i==0:
                                mag = mag_values[w][1] #the m475 value
                                text.write(' 3) '+str(md1f[w][j])+'     1\n')
                            if i==1:
                                mag = mag_values[w][2] #the m814 values
                                text.write(' 3) '+str(md2f[w][j])+'     1\n')
                            text.write(' 4) '+str(tmp[w][q])+'     1\n')
                            text.write(' 5) 4.0      1\n')
                            text.write(' 6) 0      0\n')
                            text.write(' 7) 0      0\n')
                            text.write(' 8) 0      0\n')
                            text.write(' 9) 1.0      1\n')
                            text.write(' 10) 0     1\n')
                            text.write(' Z) 0      \n')
                        
                        text.close()

#SIMULTANEOUS FITS OF VARYING DEGREE
else:
    togetherness == 'simultaneous'
    for w in range(0,ngal):
        galcoords = 'galcoords_'+plate+'.dat'
        catalog = ascii.read(galcoords)

        xcoor = str(catalog[w][1])
        ycoor = str(catalog[w][2])
        xcoorlow = str(int(catalog[w][1]-(np.float(imgsize)/2)))
        ycoorlow = str(int(catalog[w][2]-(np.float(imgsize)/2)))
        xcoorhigh = str(int(catalog[w][1]+(np.float(imgsize)/2)))
        ycoorhigh = str(int(catalog[w][2]+(np.float(imgsize)/2)))
        for q in range(0,nre):
            if not os.path.exists(time+'_'+model+'_'+plate+'_' \
            +togetherness+'_'+imgsize+'_psf'+psf):
                os.makedirs(time+'_'+model+'_'+plate+'_'+togetherness+'_'\
                +imgsize+'_psf'+psf)
            file = time+'_'+model+'_'+plate+'_'+togetherness+'_'\
            +imgsize+'_psf'+psf+'/'+galaxies[w]+\
            '_F814W,F475W,F160W_'+model+'_effective_re'+str(q)+'_input.txt'

            text = open(file,'w')

            text.write('#  Input menu file: '+galaxies[w]+\
            '_filters[475,814,160]_'+plate+'\n') #propably not essential
            text.write('#  Chi^2/nu = ,  Chi^2 = ,  Ndof = \n')
            #probably not essential
            text.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
            #will not change , probably not essential
            text.write('A) /Users/cthomps2/physics/linux-lab/data/hst/'\
            +longgal[w]+'/'+plate+'/'+filters[1]+'/final_'+filters[1]+\
            '_drc_sci.fits,/Users/cthomps2/physics/linux-lab/data/hst/'+longgal[w]+\
            '/'+plate+'/'+filters[0]+'/final_'+filters[0]+\
            '_drc_sci.fits,/Users/cthomps2/physics/linux-lab/data/hst/'+longgal[w]+\
            '/'+plate+'/'+filters[2]+'/final_'+filters[2]+'_drz_sci.fits\n')
            text.write('A1) V,U,J\n')
            text.write('A2) 814.000,475.000,160.000\n')
            text.write('B) '+galaxies[w]+'_F814W_F475W_F160W_'+model+\
            '_effective_re'+str(q)+'output.fits\n')
            text.write('C) none,none,none      0.000\n')
            text.write('D) /Users/cthomps2/physics/linux-lab/data/hst/'+longgal[w]+\
            '/'+psf+'/'+filters[1]+'/final_psf.fits,/Users/cthomps2/physics/linux-\
            lab/data/hst/'+longgal[w]+'/'+psf+'/'+filters[0]+'/final_psf.fits\
            ,/Users/cthomps2/physics/linux-lab/data/hst/'+longgal[w]+'/'+psf+'/'\
            +filters[2]+'/final_psf.fits\n')
            text.write('E) 1\n')# change to one when we are using a fine image
            text.write('F) none,none,none\n')

            if model == 'psf':
                text.write('G) \n')
                constraintfile= time+'_'+model+'_'+plate+'_'+togetherness+\
                '_'+imgsize+'_psf'+psf+'/'+'parameters_constraints_'\
                +galaxies[w]+'.txt'

                contraint_text = open(constraintfile, 'w')
                contraint_text.write('# Component/    parameter   constraint\
                Comment\n')
                contraint_text.write('# operation	(see below)   range\n\n')
                contraint_text.write('1              n              4 to 4  \
                # Soft constraint: Constrains the\n')
                contraint_text.write('					# sersic index n to\
                within values\n')
                contraint_text.write('				        # from 0.7 to\
                5.\n\n')
                contraint_text.write('1              re         '\
                +str(tmp[w][q]/2)+' to '+str(tmp[w][q]/2)+'')
                contraint_text.close()

            if model == 'sersic':
                text.write('G) \n')
                constraintfile= time+'_'+model+'_'+plate+'_'+togetherness+'_'\
                +imgsize+'_psf'+psf+'/'+'parameters_constraints_'+str([q])+\
                '_'+galaxies[w]+'.txt'

                contraint_text = open(constraintfile, 'w')
                contraint_text.write('# Component/    parameter   constraint\
                Comment\n')
                contraint_text.write('# operation	(see below)   range\n\n')
                contraint_text.write('1              n              4 to 4\
                # Soft constraint: Constrains the\n')
                contraint_text.write('					# sersic index n to\
                within values\n')
                contraint_text.write('				        # from 0.7 to\
                5.\n\n')
                contraint_text.write('1              re         '\
                +str(tmp[w][q]/2)+' to '+str(tmp[w][q]/2)+'')
                contraint_text.write('1              q           \
                0.99 to 0.99  \n')
                contraint_text.write('1              pa          \
                0 to 0  \n')
                contraint_text.close()

            text.write('G) parameters_constraints_'+str([q])+'_'\
            +galaxies[w]+'.txt\n')

            galcoords = 'galcoords_'+plate+'.dat'
            catalog = ascii.read(galcoords)

            xcoor = str(catalog[w][1])
            ycoor = str(catalog[w][2])
            xcoorlow = str(int(catalog[w][1]-(np.float(imgsize)/2)))
            ycoorlow = str(int(catalog[w][2]-(np.float(imgsize)/2)))
            xcoorhigh = str(int(catalog[w][1]+(np.float(imgsize)/2)))
            ycoorhigh = str(int(catalog[w][2]+(np.float(imgsize)/2)))

            text.write('H) '+xcoorlow+' '+xcoorhigh+' '+ycoorlow+' '\
            +ycoorhigh+'\n')
            text.write('I) 150    150\n')
            text.write('J) 25.027,25.613\n')
            text.write('J) 25.027,25.613,25.946\n')
            text.write('K) '+plate+'  '+plate+'\n')
            text.write('O) regular\n')
            text.write('P) 0\n')
            text.write('U) 0 0.750000 25 4 40\n')
            text.write('V) 0 0 50 0.800000 0.500000 100000\n')
            text.write('W) input,sigma,psf,component,model,residual\n')
            
            # Psf for this code does not work.....
            if model == 'psf':
                text.write(' 0) psf\n')
                if togetherness == 'simultaneous':
                    text.write(' 1) '+xcoor+','+xcoor+'    1,1              \
                    band\n')
                    text.write(' 2) '+ycoor+','+ycoor+'    1,1              \
                    band\n')
                if togetherness == 'semi':
                    text.write(' 1) '+xcoor+','+xcoor+'    1,0              \
                    band\n')
                    text.write(' 2) '+ycoor+','+ycoor+'    1,0              \
                    band\n')
                text.write(' 3) 19.5,19.5     1,1                 band\n')
                text.write(' Z) 0                  #  Skip this model in \
                output image?  (yes=1, no=0)\n')

            #Sersic for the code does work...... 
            if model == 'sersic':
                text.write(' 0) sersic\n')
                if togetherness == 'simultaneous':
                    text.write(' 1) '+xcoor+','+xcoor+','+xcoor+'    1,1,1  \
                    band\n')
                    text.write(' 2) '+ycoor+','+ycoor+','+ycoor+'    1,1,1  \
                    band\n')
                if togetherness == 'semi':
                    text.write(' 1) '+xcoor+','+xcoor+'    1,0              \
                    band\n')
                    text.write(' 2) '+ycoor+','+ycoor+'    1,0              \
                    band\n')
                text.write(' 3) 19.5,19.5,18.8     1,1,1                \
                band\n') 
                text.write(' 4) '+str(tmp[w][q]/2)+','+str(tmp[w][q]/2)+\
                ','+str(tmp[w][q]/2)+'    1,0,0                band\n')
                text.write(' 5) 4.000,4.000,4.000    1,0,0               \
                band\n')
                text.write(' 6) 0,0,0               0,0,0                \
                band\n')
                text.write(' 7) 0,0,0              0,0,0                \
                band\n')
                text.write(' 8) 0,0,0               0,0,0                \
                band\n')
                text.write(' 9) 0.99,0.99,0.99           1,0,0            \
                band\n')
                text.write(' 10) 0,0,0          1,0,0                 band\n')
                text.write(' Z) 0')
            text.close()


#below we create shell scripts, sh, to run galfit on the input files. NOTE! 
#that the only one that are written to produce a proper sh, shell file, have a 
# model of 'sersic' and a togetherness equal to 'independent' or
#'simulataneous'.
if model == 'psf' and togetherness == 'independent':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+\
    '/'+'psf_independent_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        for i in range(0,2):
            text.write('galfitm '+galaxies[w]+'_'+filters[i]+'_'+model+\
            '_input.txt\n')
    text.close()

if model == 'sersic' and togetherness == 'independent':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+\
    '/'+'sersic_independent_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,ngal):
        for i in range(0,2):
            for j in range (0,nmag):
                for q in range (0,nre):
                    text.write('galfitm '+galaxies[w]+'_'+filters2[i]+\
                    '_'+model+'_effective_re'+str(q)+'_magnitude'+str(j)+\
                    '_input.txt\n')
    text.close()

if model == 'psf' and togetherness == 'simultaneous':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+\
    '/'+'psf_simultaneous_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        text.write('galfitm '+galaxies[w]+'_F814WandF475W_'+model+\
        '_input.txt\n')
    text.close()

if model == 'sersic' and togetherness == 'simultaneous':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+\
    '/''sersic_simultaneous_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,ngal):
        for q in range (0,nre):
            text.write('galfitm '+galaxies[w]+'_F814W,F475W,F160W_'+model+\
            '_effective_re'+str(q)+'_input.txt\n')
    text.close()

if model == 'psf' and togetherness == 'semi':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+\
    '/'+'psf_semi_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        text.write('galfitm '+galaxies[w]+'_F814WandF475W_'+model+\
        '_input.txt\n')
    text.close()

if model == 'sersic' and togetherness == 'semi':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+\
    '/'+'sersic_semi_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        text.write('galfitm '+galaxies[w]+'_F814WandF475W_'+model+\
        '_input.txt\n')
    text.close()
