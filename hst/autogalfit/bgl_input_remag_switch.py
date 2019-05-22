# automating the creation and running of galfit input files
# this version is for when we want to explore parameter space for re and magnitude

import os
import numpy as np
from astropy.io import ascii
import datetime
import sys
import glob

#sys.argv[1] is psf or sersic (model), sys.arg[2] is coarse or fine (plate scale), sys.arg[3] is a number (image region size), sys.arg[4] is coarse or fine (psf), sys.arg[5] is simultaneous or independent or semi, sys.arg[6] is either none or s_index or re or s_index, re

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']
longgal = ['J0826+4305','J0901+0314','J0905+5759','J0944+0930','J1107+0417','J1219+0336','J1341-0321','J1506+5402','J1558+3957','J1613+2834','J2116-0634','J2140+1209']
photzeromag = ['25.613','25.027']
kind = ['independent','simultaneous']

model = sys.argv[1]
plate = sys.argv[2]
imgsize = sys.argv[3]
psf = sys.argv[4]
togetherness = sys.argv[5]
constraint = sys.argv[6]

now = datetime.datetime.now()
time = now.strftime("%Y%m%d-%H%M")

bgl_values = ascii.read('bgl_mag_size_table.dat')

re_array = np.array([1/3, 1/2, 1/1.5, 1/1.2, 1, 1.2, 1.5, 2, 3])
mag_array = np.array([-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2])

ngal = len(bgl_values)
nre = len(re_array)
nmag = len(mag_array)

# md1f stands for "magnitude difference of 1st filter" (i.e., F475W) 
md1f = np.zeros((ngal, nmag))
# md2f stands for "magnitude difference of 2nd filter" (i.e., F814W).
md2f = np.zeros((ngal, nmag))

for i in range(0, ngal):
    for j in range(0, nmag):
        md1f[i,j] = bgl_values['m475'][i] + mag_array[j]
        md2f[i,j] = bgl_values['m814'][i] + mag_array[j]

tmp = np.zeros((ngal, nre))
for i in range(0,ngal):
    for j in range(0,nre):
        tmp[i,j] = bgl_values['re'][i] * re_array[j]
        
hstdir = os.environ['HSTDIR']

#INDEPENDENT FITS YO
if togetherness == 'independent':
    for w in range(0,12):
        for i in range(0,2):
            sigma_file = glob.glob(hstdir+'sigma/'+longgal[w]+'_'+filters[i]+'_'+plate+'_sigma.fits')

            if not os.path.exists(time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf):
                os.makedirs(time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf)


            file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+galaxies[w]+'_'+filters[i]+'_'+model+'_input.txt'
            text = open(file,'w')
            text.write('#  Input menu file: '+galaxies[w]+'_'+filters[i]+'_'+plate+'\n') #probably not essential
            text.write('#  Chi^2/nu = ,  Chi^2 = ,  Ndof = \n') #probably not essential
            text.write('# IMAGE and GALFIT CONTROL PARAMETERS\n') #will not change , probably not essential
            text.write('A) '+hstdir+longgal[w]+'/'+plate+'/'+filters[i]+'/final_'+filters[i]+'_drc_sci.fits\n')
            text.write('B) '+galaxies[w]+'_'+filters[i]+'_output.fits\n')
            text.write('C) '+sigma_file[0]+' \n') 
            #text.write('C) none\n') #will not change
            text.write('D) '+hstdir+longgal[w]+'/'+psf+'/'+filters[i]+'/final_psf.fits\n')
            if (plate == 'fine' and psf == 'fine'):
                text.write('E) 1\n') 
            if (plate == 'fine' and psf == 'coarse'):
                text.write('E) 0.5\n') 
            if (plate == 'coarse' and psf == 'fine'):
                text.write('E) 2\n') 
            if (plate == 'coarse' and psf == 'coarse'):
                text.write('E) 1\n')
            text.write('F) none\n') #will not change

            #if constraint == 'none':
            #    text.write('G) \n')
            #if constraint == 's_index':
            #    text.write('G) '+hstdir+'../sersic_index_equaltofour.txt\n')
            #if constraint == 're':
            #     text.write('G) '+hstdir+'../size_constraint_1to1.txt\n')
            #if constraint == 's_index, re':
            #     text.write('G) '+hstdir+'../s_index_and_size.txt\n')    


            if model == 'sersic':
                constraintfile= time+'_'+model+'_'+plate+'_'+togetherness+'_'\
                +imgsize+'_psf'+psf+'/'+'parameters_constraints_'+galaxies[w]+'.txt'

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
                contraint_text.write('1              q           '+\
                str(bgl_values['ba'][w])+' to '+str(bgl_values['ba'][w])+'  \n')
                contraint_text.write('1              pa          '+\
                str(bgl_values['pa'][w])+' to '+str(bgl_values['pa'][w])+'  \n')
                contraint_text.close()

                text.write('G) parameters_constraints_'+galaxies[w]+'.txt\n')
                
            if model == 'psf':
                text.write('G) \n')
                 
            galcoords = 'galcoords_'+plate+'.dat'
            catalog = ascii.read(galcoords)
        
            xcoor = str(catalog[w][1])
            ycoor = str(catalog[w][2])
            xcoorlow = str(int(catalog[w][1]-(np.float(imgsize)/2)))
            ycoorlow = str(int(catalog[w][2]-(np.float(imgsize)/2)))
            xcoorhigh = str(int(catalog[w][1]+(np.float(imgsize)/2)))
            ycoorhigh = str(int(catalog[w][2]+(np.float(imgsize)/2)))
        
            text.write('H) '+xcoorlow+' '+xcoorhigh+' '+ycoorlow+' '+ycoorhigh+'\n')
            text.write('I) 150    150\n') 
            text.write('J) '+photzeromag[i]+'\n')

            if plate == 'coarse':
                text.write('K) 0.050  0.050\n')
            if plate == 'fine':
                text.write('K) 0.025  0.025\n')
            
            text.write('O) regular\n') #will not change
            text.write('P) 0\n') #will not change
            text.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
            text.write('# Component number: 1\n')

            if model == 'psf':
                text.write(' 0) psf\n')
                text.write(' 1) '+xcoor+' '+ycoor+' 1 1  # position x, y        [pixel]\n')
                if i==0:
                    mag=bgl_values[w][1] # the m475 value
                    text.write(' 3) '+str(mag)+'     1\n')
                if i==1:
                    mag=bgl_values[w][2] # the m814 value
                    text.write(' 3) '+str(mag)+'     1\n')
                text.write(' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n')
            if model == 'sersic':
                text.write(' 0) sersic\n')
                text.write(' 1) '+xcoor+' '+ycoor+' 1 1  # position x, y        [pixel]\n')
                if i==0:
                    mag=bgl_values[w][1] # the m475 value
                    text.write(' 3) '+str(mag)+'     1\n')
                if i==1:
                    mag=bgl_values[w][2] # the m814 value
                    text.write(' 3) '+str(mag)+'     1\n')
                text.write(' 4) '+str(bgl_values[w][3])+'      1\n')
                text.write(' 5) 4.0      1\n')
                text.write(' 6) 0      0\n')
                text.write(' 7) 0      0\n')
                text.write(' 8) 0      0\n')
                text.write(' 9) '+str(bgl_values['ba'][w])+'      1\n')
                text.write(' 10) '+str(bgl_values['pa'][w])+'     1\n')
        
                text.write(' Z) 0      \n')
                    
            text.close()



#SIMULTANEOUS FITS OF VARYING DEGREE
else:
    now = datetime.datetime.now()
    time = now.strftime("%Y%m%d-%H%M")
    for w in range(0,ngal):
        galcoords = 'galcoords_'+plate+'.dat'
        catalog = ascii.read(galcoords)
                
        xcoor = str(catalog[w][1])
        ycoor = str(catalog[w][2])
        xcoorlow = str(int(catalog[w][1]-(np.float(imgsize)/2)))
        ycoorlow = str(int(catalog[w][2]-(np.float(imgsize)/2)))
        xcoorhigh = str(int(catalog[w][1]+(np.float(imgsize)/2)))
        ycoorhigh = str(int(catalog[w][2]+(np.float(imgsize)/2)))

        sigma_file_814 = glob.glob(hstdir+'sigma/'+longgal[w]+'_'+filters[1]+'_'+plate+'_sigma.fits')
        sigma_file_475 = glob.glob(hstdir+'sigma/'+longgal[w]+'_'+filters[0]+'_'+plate+'_sigma.fits')
        if not os.path.exists(time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf):
            os.makedirs(time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf)

        for j in range(0, nmag):
            for q in range(0,nre):
                
                file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+galaxies[w]+'_F814WandF475W_'+model+'_re'+str(q)+'_mag'+str(j)+'_input.txt'
                text = open(file,'w')
    
                text.write('#  Input menu file: '+galaxies[w]+'_F814WandF475W_'+model+'_re'+str(q)+'_mag'+str(j)+'_input.txt\n')
                text.write('A) '+hstdir+longgal[w]+'/'+plate+'/'+filters[1]+'/final_'+filters[1]+'_drc_sci.fits,'+hstdir+longgal[w]+'/'+plate+'/'+filters[0]+'/final_'+filters[0]+'_drc_sci.fits\n')
                text.write('A1) V,U\n')
                text.write('A2) 814.000,475.000\n')
                text.write('B) '+galaxies[w]+'_F814WandF475W_'+model+'_re'+str(q)+'_mag'+str(j)+'_output.fits\n')
                text.write('C) '+sigma_file_814[0]+','+sigma_file_475[0]+' \n') 
                text.write('D) '+hstdir+longgal[w]+'/'+psf+'/'+filters[1]+'/final_psf.fits,'+hstdir+longgal[w]+'/'+psf+'/'+filters[0]+'/final_psf.fits\n')
                if (plate == 'fine' and psf == 'fine'):
                    text.write('E) 1\n') 
                if (plate == 'fine' and psf == 'coarse'):
                    text.write('E) 0.5\n') 
                if (plate == 'coarse' and psf == 'fine'):
                    text.write('E) 2\n') 
                if (plate == 'coarse' and psf == 'coarse'):
                    text.write('E) 1\n')
                text.write('F) none,none\n')

                if model == 'sersic':
                    #text.write('G) \n')
                    constraintfile= time+'_'+model+'_'+plate+'_'+togetherness+'_'\
                    +imgsize+'_psf'+psf+'/'+'parameters_constraints_'+galaxies[w]+'_re'+str(q)+'_mag'+str(j)+'.txt'
                
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
                        +str(tmp[w][q])+' to '+str(tmp[w][q])+' \n')
                    contraint_text.write('1              MAG_U     \
                            '+str(md1f[w][j])+' to '+str(md1f[w][j])+' \n')
                    contraint_text.write('1              MAG_V     \
                            '+str(bgl_values['m814'][w])+' to '+str(bgl_values['m814'][w])+' \n')                            
                    contraint_text.write('1              q           '+\
                    str(bgl_values['ba'][w])+' to '+str(bgl_values['ba'][w])+'  \n')
                    contraint_text.write('1              pa          '+\
                    str(bgl_values['pa'][w])+' to '+str(bgl_values['pa'][w])+'  \n')
                    contraint_text.close()
                
                    text.write('G) parameters_constraints_'+galaxies[w]+'_re'+str(q)+'_mag'+str(j)+'.txt\n')
                    
                if model == 'psf':
                    text.write('G) \n')
                
        
                text.write('H) '+xcoorlow+' '+xcoorhigh+' '+ycoorlow+' '+ycoorhigh+'\n')
                text.write('I) 150    150\n')
                text.write('J) 25.027,25.613\n')
                text.write('K) '+plate+'  '+plate+'\n')
                text.write('O) regular\n')
                text.write('P) 0\n')
                text.write('U) 0 0.750000 25 4 40\n')
                text.write('V) 0 0 50 0.800000 0.500000 100000\n')
                text.write('W) input,sigma,psf,component,model,residual\n')

                if model == 'psf':
                    text.write(' 0) psf\n')
                    if togetherness == 'simultaneous':
                        text.write(' 1) '+xcoor+','+xcoor+'    1,1                 band\n')
                        text.write(' 2) '+ycoor+','+ycoor+'    1,1                 band\n')
                    if togetherness == 'semi':
                        text.write(' 1) '+xcoor+','+xcoor+'    1,0                 band\n')
                        text.write(' 2) '+ycoor+','+ycoor+'    1,0                 band\n')
                    text.write(' 3) '+str(bgl_values[w][2])+','+str(bgl_values[w][1])+'     1,1                 band\n')
                    text.write(' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n')
                
                if model == 'sersic':
                    text.write(' 0) sersic\n')
                    if togetherness == 'simultaneous':
                        text.write(' 1) '+xcoor+','+xcoor+'    1,1                 band\n')
                        text.write(' 2) '+ycoor+','+ycoor+'    1,1                 band\n')
                    if togetherness == 'semi':
                        text.write(' 1) '+xcoor+','+xcoor+'    1,0                 band\n')
                        text.write(' 2) '+ycoor+','+ycoor+'    1,0                 band\n')
                    text.write(' 3) '+str(bgl_values['m814'][w])+','+str(md1f[w][j])+'     1,1                 band\n')
                    text.write(' 4)  '+str(tmp[w][q])+',0    1,0                 cheb\n')
                    text.write(' 5) 4.000,0    1,0                 cheb\n')
                    text.write(' 6) 0,0               0,0                 cheb\n')
                    text.write(' 7) 0,0               0,0                 cheb\n')
                    text.write(' 8) 0,0               0,0                 cheb\n')
                    text.write(' 9) '+str(bgl_values['ba'][w])+',0      1,0  cheb\n')
                    text.write(' 10) '+str(bgl_values['pa'][w])+',0     1,0  cheb\n')
                    text.write(' Z) 0')
                text.close()


#below we create shell scripts to run galfit on the input files
if model == 'psf' and togetherness == 'independent':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+'psf_independent_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        for i in range(0,2):
            text.write('galfitm '+galaxies[w]+'_'+filters[i]+'_'+model+'_input.txt\n')
    text.close()
    
if model == 'sersic' and togetherness == 'independent':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+'sersic_independent_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        for i in range(0,2):
            text.write('galfitm '+galaxies[w]+'_'+filters[i]+'_'+model+'_input.txt\n')
    text.close()

if model == 'psf' and togetherness == 'simultaneous':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+'psf_simultaneous_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,ngal):
        text.write('galfitm '+galaxies[w]+'_F814WandF475W_'+model+'_input.txt\n')
    text.close()

if model == 'sersic' and togetherness == 'simultaneous':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+'sersic_simultaneous_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,ngal):
        for j in range(0,nmag):
            for q in range(0,nre):
                text.write('galfitm '+galaxies[w]+'_F814WandF475W_'+model+'_re'+str(q)+'_mag'+str(j)+'_input.txt\n')
    text.close()

if model == 'psf' and togetherness == 'semi':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+'psf_semi_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        text.write('galfitm '+galaxies[w]+'_F814WandF475W_'+model+'_input.txt\n')
    text.close()

if model == 'sersic' and togetherness == 'semi':
    file = time+'_'+model+'_'+plate+'_'+togetherness+'_'+imgsize+'_psf'+psf+'/'+'sersic_semi_run_files.sh'
    text = open(file,'w')
    text.write('shopt -s expand_aliases\n')
    text.write('source ~/.bash_profile\n')
    for w in range(0,12):
        text.write('galfitm '+galaxies[w]+'_F814WandF475W_'+model+'_input.txt\n')
    text.close()


