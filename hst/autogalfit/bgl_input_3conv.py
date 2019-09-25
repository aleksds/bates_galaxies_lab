# automating the creation and running of galfit input files
# this version is for when we perform a 3-band fit with constrained re values

import os
import numpy as np
from astropy.io import ascii
import datetime
import sys
import glob


galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W','F160W']
longgal = ['J0826+4305','J0901+0314','J0905+5759','J0944+0930','J1107+0417','J1219+0336','J1341-0321','J1506+5402','J1558+3957','J1613+2834','J2116-0634','J2140+1209']
photzeromag = ['25.613','25.027','25.946']


#model = sys.argv[1]
#plate = sys.argv[2]
#imgsize = sys.argv[3]
#psf = sys.argv[4]
#togetherness = sys.argv[5]
#constraint = sys.argv[6]

plate = 'coarse'
psf = 'coarse'
imgsize = 200

bgl_values = ascii.read('bgl_mag_size_table.dat')

re_vals = bgl_values['re03']
ba_vals = bgl_values['ba']
pa_vals = bgl_values['pa']

data = ascii.read('bgl_phot.dat')
m475 = data['m475']
m814 = data['m814']
m160 = data['m160']

ngal = len(bgl_values)

hstdir = os.environ['HSTDIR']



now = datetime.datetime.now()
time = now.strftime("%Y%m%d-%H%M")

for w in range(0,ngal):
    sigma_file_814 = glob.glob(hstdir+'sigma/'+longgal[w]+'_'+filters[1]+'_'+plate+'_sigma.fits')
    sigma_file_475 = glob.glob(hstdir+'sigma/'+longgal[w]+'_'+filters[0]+'_'+plate+'_sigma.fits')
    sigma_file_160 = glob.glob(hstdir+'sigma/'+longgal[w]+'_'+filters[2]+'_'+plate+'_sigma.fits')
    if not os.path.exists(time+'_3conv'):
        os.makedirs(time+'_3conv')

    file = time+'_3conv'+'/'+galaxies[w]+'_F814W_F475W_F160W_3conv_input.txt'
    text = open(file,'w')
                
    text.write('#  Input menu file: '+galaxies[w]+'_F814W_F475W_F160W_3conv_input.txt\n') 
    text.write('#  Chi^2/nu = ,  Chi^2 = ,  Ndof = \n') 
    text.write('# IMAGE and GALFIT CONTROL PARAMETERS\n') 
                
    text.write('A) '+hstdir+longgal[w]+'/'+plate+'/'+filters[1]+'/final_'+filters[1]+'_drc_sci.fits,'+
                     #hstdir+longgal[w]+'/'+plate+'/'+filters[1]+'/convolved_image_'+filters[1]+'.fits,'+
                     hstdir+longgal[w]+'/'+plate+'/'+filters[0]+'/final_'+filters[0]+'_drc_sci.fits,'+
                     #hstdir+longgal[w]+'/'+plate+'/'+filters[0]+'/convolved_image_'+filters[0]+'.fits,'+
                     hstdir+longgal[w]+'/'+plate+'/'+filters[2]+'/final_'+filters[2]+'_drz_sci.fits\n')
    text.write('A1) V,U,J\n')
    text.write('A2) 814.000,475.000,160.000\n')
                
    text.write('B) '+galaxies[w]+'_F814W_F475W_F160W_3conv_output.fits\n')
    text.write('C) '+sigma_file_814[0]+','+sigma_file_475[0]+','+sigma_file_160[0]+' \n') 
    text.write('D) '+hstdir+longgal[w]+'/'+psf+'/'+filters[1]+'/final_psf.fits,'+
                     #hstdir+longgal[w]+'/'+psf+'/'+filters[1]+'/psf_convol.fits,'+
                     hstdir+longgal[w]+'/'+psf+'/'+filters[0]+'/final_psf.fits,'+
                     #hstdir+longgal[w]+'/'+psf+'/'+filters[0]+'/psf_convol.fits,'+
                     hstdir+longgal[w]+'/'+psf+'/'+filters[2]+'/final_psf.fits\n')
    text.write('E) 1\n')
    text.write('F) none,none,none\n')
 
    text.write('G) \n')
    constraintfile= time+'_3conv'+'/'+'parameters_constraints_'+galaxies[w]+'.txt'
                
    contraint_text = open(constraintfile, 'w')
    contraint_text.write('# Component/    parameter   constraint Comment\n')
    contraint_text.write('# operation	(see below)   range\n')
    contraint_text.write('1              n              4 to 4 # Soft constraint: Constrains the\n')
    contraint_text.write('					# sersic index n to within values\n')
    contraint_text.write('				        # from 0.7 to 5. \n')
    contraint_text.write('1              re         '+str(re_vals[w]/2)+' to '+str(re_vals[w]/2)+' \n')
    contraint_text.write('1              MAG_U      '+str(m475[w])+' to '+str(m475[w])+' \n')
    contraint_text.write('1              MAG_V      '+str(m814[w])+' to '+str(m814[w])+' \n')
    contraint_text.write('1              MAG_J      '+str(m160[w])+' to '+str(m160[w])+' \n')  
    contraint_text.write('1              q           '+str(ba_vals[w])+' to '+str(ba_vals[w])+'  \n')
    contraint_text.write('1              pa          '+str(pa_vals[w])+' to '+str(pa_vals[w])+'  \n')
    contraint_text.close()
                
    text.write('G) parameters_constraints_'+galaxies[w]+'.txt\n')                
                
    galcoords = 'galcoords_'+plate+'.dat'
    catalog = ascii.read(galcoords)
                    
    xcoor = str(catalog[w][1])
    ycoor = str(catalog[w][2])
    xcoorlow = str(int(catalog[w][1]-(np.float(imgsize)/2)))
    ycoorlow = str(int(catalog[w][2]-(np.float(imgsize)/2)))
    xcoorhigh = str(int(catalog[w][1]+(np.float(imgsize)/2)))
    ycoorhigh = str(int(catalog[w][2]+(np.float(imgsize)/2)))
                
    text.write('H) '+xcoorlow+' '+xcoorhigh+' '+ycoorlow+' '+ycoorhigh+'\n')
    #text.write('I) '+str(np.float(imgsize)/2)+'    '+str(np.float(imgsize)/2)+'\n')
    text.write('I) 150     150   \n')
    text.write('J) 25.027,25.613,25.946\n')
    text.write('K) 0.05  0.05  \n')
    text.write('O) regular\n')
    text.write('P) 0\n')
    #text.write('U) 0 0.750000 25 4 40\n')
    #text.write('V) 0 0 50 0.800000 0.500000 100000\n')
    text.write('W) input,sigma,psf,component,model,residual\n')
                
    text.write(' 0) sersic\n')
    text.write(' 1) '+xcoor+','+xcoor+','+xcoor+'    1,1,1  band\n')
    text.write(' 2) '+ycoor+','+ycoor+','+ycoor+'    1,1,1  band\n')
                    
    text.write(' 3) '+str(m814[w])+','+str(m475[w])+','+str(m160[w])+'     1,1,1    band\n')
    text.write(' 4) '+str(re_vals[w]/2)+','+str(re_vals[w]/2)+','+str(re_vals[w]/2)+'    1,0,0                band\n')
    text.write(' 5) 4.000,4.000,4.000    1,0,0     band\n')
    text.write(' 6) 0,0,0               0,0,0      band\n')
    text.write(' 7) 0,0,0              0,0,0       band\n')
    text.write(' 8) 0,0,0               0,0,0      band\n')
    text.write(' 9) '+str(bgl_values['ba'][w])+',0,0   1,0,0  cheb\n')
    text.write(' 10) '+str(bgl_values['pa'][w])+',0,0  1,0,0  cheb\n')

    text.close()


file = time+'_3conv'+'/'+'run_file.sh'
text = open(file,'w')
text.write('shopt -s expand_aliases\n')
text.write('source ~/.bash_profile\n')
for w in range(0,ngal):
    text.write('galfitm '+galaxies[w]+'_F814W_F475W_F160W_3conv_input.txt\n')
text.close()



