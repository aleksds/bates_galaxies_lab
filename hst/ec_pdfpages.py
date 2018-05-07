#This program creates a pdf with the following elements:
    #color-magnitude plot (sersic and psf)
    #size-size plot (sersic)
    #data/model/residual plot (zoomed in to 5x5 central section) (sersic and psf)

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

#sersic values:
#These values are from independant trials for each filter
re_F475W_indep_s = 0.272
re_F814W_indep_s = 0.432

#(semi-simultaneous)These values are from galfitm where x and y values are allowed to be different between filters
re_F475W_simul_s = 0.337
re_F814W_simul_s = 0.337 

#These values are from galfitm where x and y positions must remain constant
re_F475W_simul2_s = 9.258*10**(-2)
re_F814W_simul2_s = 9.258*10**(-2)

re_diff_s = [(re_F475W_indep_s/re_F814W_indep_s), (re_F475W_simul_s/re_F814W_simul_s), (re_F475W_simul2_s/re_F814W_simul2_s)]
F814W_s = [re_F814W_indep_s, re_F814W_simul_s, re_F814W_simul2_s]

#These values are from independant trials for each filter
m_F475W_indep_s = 19.497
m_F814W_indep_s = 19.089

#(semi-simultaneous) These values are from galfitm where x and y values are allowed to be different between filters
m_F475W_simul_s = 19.537  
m_F814W_simul_s = 19.101

#These values are from galfitm where x and y positions must remain constant
m_F475W_simul2_s = 19.396
m_F814W_simul2_s = 18.966

m_diff_s = [(m_F475W_indep_s - m_F814W_indep_s), (m_F475W_simul_s - m_F814W_simul_s), (m_F475W_simul2_s - m_F814W_simul2_s)]
F814W_s = [m_F814W_indep_s, m_F814W_simul_s, m_F814W_simul2_s]

#psf values:
#These values are from independant trials for each filter
m_F475W_indep_p = 19.626
m_F814W_indep_p = 19.224

#(semi-simultaneous) These values are from galfitm where x and y values are allowed to be different between filters
m_F475W_simul_p = 19.597  
m_F814W_simul_p = 19.203

#These values are from galfitm where x and y positions must remain constant
m_F475W_simul2_p = 19.663
m_F814W_simul2_p = 19.231

m_diff_p = [(m_F475W_indep_p - m_F814W_indep_p), (m_F475W_simul_p - m_F814W_simul_p), (m_F475W_simul2_p - m_F814W_simul2_p)]
F814W_p = [m_F814W_indep_p, m_F814W_simul_p, m_F814W_simul2_p]

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
dx = dy = 5
with PdfPages('ec_pdfpages.pdf') as pdf:
    fig = plt.figure()
    plt.suptitle('J0905 Individual Model')
    alldata = []
    for f in range(0,2):
        for t in range(0,len(type)):
            file = glob.glob('/Volumes/physics/linux-lab/data/galfit/eves_files/sersic_fine/ec_J0905_'+filename[f]+'_fine_sersic_output.fits') 
                             
            multi = fits.open(file[0])
            print(file)
            
            data, data_header = multi[1].data, multi[1].header
            model, res_header = multi[2].data, multi[2].header
            res, res_header = multi[3].data, multi[3].header


            stampdata = data[round(409-dy*10):round(409+dy*10), round(415-dx*10):round(415+dx*10)] 
            stampmodel = model[round(409-dy*10):round(409+dy*10), round(415-dx*10):round(415+dx*10)]
            stampres = res[round(409-dy):round(409+dy), round(415-dx):round(415+dx)]
            
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

    dx = dy = 5
    fig = plt.figure()
    plt.suptitle('J0905 Simultaneous Model')
    alldata = []
    for f in range(0,2):
        for t in range(0,len(type)):
            file = glob.glob('/Volumes/physics/linux-lab/data/galfit/eves_files/sersic_fine/ec_J0905_'+filename[2]+'_sersic_output.fits') 
                             
            multi = fits.open(file[0])
            
            data, data_header = multi[f].data, multi[f].header
            model, res_header = multi[2+f].data, multi[2+f].header
            res, res_header = multi[4+f].data, multi[4+f].header


            stampdata = data[round(99-dy*10):round(99+dy*10), round(95-dx*10):round(95+dx*10)] 
            stampmodel = model[round(99-dy*10):round(99+dy*10), round(95-dx*10):round(95+dx*10)]
            stampres = res[round(99-dy):round(99+dy), round(95-dx):round(95+dx)]
            
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

    fig = plt.figure()
    plt.suptitle('J0905 Difference in Magnitude Between Filters')
    ax = fig.add_subplot(2,1,1)
    plt.xlabel('Integrated Magnitude of F814W')
    plt.ylabel('m_F475W - m_F814W')
    labels = ['sersic independent (r_e, x, y, magnitude allowed to vary)', 'sersic semi-simultaneous (x, y, magnitude allowed to vary)', 'sersic simultaneous (magnitude allowed to vary)','psf independent (r_e, x, y, magnitude allowed to vary)', 'psf semi-simultaneous (x, y, magnitude allowed to vary)', 'psf simultaneous (magnitude allowed to vary)']
    colors = ['lightcoral', 'firebrick', 'darkorange', 'palegreen', 'lightseagreen', 'lightskyblue']
    markers = ['o', 'v','s','o', 'v','s']
    for i in range(0,6):
        if i < 3:
            plt.scatter(F814W_s[i],m_diff_s[i], label = labels[i], color = colors[i], marker = markers[i])
        if i >= 3:
            plt.scatter(F814W_p[i-3],m_diff_p[i-3], label = labels[i], color = colors[i], marker = markers[i])
    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.25))

    fig = plt.figure()
    plt.suptitle('J0905 Difference in Effective Radius Between Filters')
    ax = fig.add_subplot(2,1,1)
    plt.xlabel('Integrated Magnitude of F814W')
    plt.ylabel('re_F475W/reF814W')

    labels = ['sersic independent (r_e, x, y, magnitude allowed to vary)', 'sersic semi-simultaneous (x, y, magnitude allowed to vary)', 'sersic simultaneous (magnitude allowed to vary)', 'psf independent', 'psf semi-simultaneous', 'psf simultaneous']
    colors = ['lightcoral', 'firebrick', 'darkorange', 'palegreen', 'lightseagreen', 'lightskyblue']
    markers = ['o', 'v','s']
    for i in range(0,3):
        plt.scatter(F814W_s[i],re_diff_s[i], label = labels[i], color = colors[i], marker = markers[i])
    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.25))


    pdf.savefig(dpi=1000)
    plt.close()

    pdf.savefig(dpi=1000)
    plt.close()
os.system('open %s &' % 'ec_pdfpages.pdf')




