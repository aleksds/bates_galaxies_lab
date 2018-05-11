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

#read in values from file:
#currently only looks at magnitude: will be issues reading in radius - only in sersic files, sometimes is in scientific notation (more characters)
j = 0
magnitude = np.zeros(12)
chi = np.zeros(8)
filenames = ['F475W_psf','F475W_sersic','F814W_psf','F814W_sersic','F814WandF475W_psf_full','F814WandF475W_psf_semi','F814WandF475W_sersic_full','F814WandF475W_sersic_semi']
for i in range(0,8):
    file = '/Volumes/physics/linux-lab/data/galfit/eves_files/pdf_pages/J0905_'+filenames[i]+'_output.galfit.01.band'
    with open(file) as f:
        content = f.readlines()
        chi[i] = str(content[3][14:19])
        if i < 4:
            magnitude[j] = np.float(content[47][4:10])
            j += 1
        else:
            magnitude[j] = np.float(content[47][4:10])
            j += 1
            magnitude[j] = np.float(content[47][11:17])
            j += 1
        
m_diff_s = [magnitude[1]-magnitude[3], magnitude[9]-magnitude[8], magnitude[11]-magnitude[10]]
m_diff_p = [magnitude[0]-magnitude[2], magnitude[5]-magnitude[4], magnitude[7]-magnitude[6]]

F814W_s = [magnitude[3], magnitude[8], magnitude[10]]
F814W_p = [magnitude[2], magnitude[4], magnitude[6]]


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
    plt.imshow(img_scale.log(rotated, scale_min=5, scale_max=10000), cmap='Greys', interpolation='nearest')


type = ['data','model','residual']
dx = dy = 6
plotnames = ['psf independent', 'sersic independent', 'psf semi-simultaneous', 'psf fully simultaneous', 'sersic semi-simultaneous', 'sersic fully simultaneous']
counter = 0
with PdfPages('ec_pdfpages.pdf') as pdf: #
    for i in range(0, len(plotnames)):
        fig = plt.figure()
        plt.suptitle('J0905 '+plotnames[i])
        for f in range(0,2):
            for t in range(0, len(type)):
                if i == 0: #when i = 0, we want to open filenames[0] when f = 0 and fn[2] when f = 1
                    file = glob.glob('/Volumes/physics/linux-lab/data/galfit/eves_files/pdf_pages/J0905_'+filenames[f*2]+'_output.fits')
                elif i == 1: #whhen i = 1, we want to open filenames[1] when f = 0 and fn[3] when f = 1
                    file = glob.glob('/Volumes/physics/linux-lab/data/galfit/eves_files/pdf_pages/J0905_'+filenames[1+f*2]+'_output.fits')                   
                else:
                    file = glob.glob('/Volumes/physics/linux-lab/data/galfit/eves_files/pdf_pages/J0905_'+filenames[i+2]+'_output.fits')
                multi = fits.open(file[0])
                if i < 2:
                    data, data_header = multi[1].data, multi[1].header
                    model, res_header = multi[2].data, multi[2].header
                    res, res_header = multi[3].data, multi[3].header
                else:
                    data, data_header = multi[1-f].data, multi[1-f].header
                    model, res_header = multi[3-f].data, multi[3-f].header
                    res, res_header = multi[5-f].data, multi[5-f].header            
    
            if i == 0 | 1:
                stampdata = data[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)] 
                stampmodel = model[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)]
                stampres = res[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)]
            elif i == 2 | 3:
                stampdata = data[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)] 
                stampmodel = model[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)]
                stampres = res[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)]
            else:
                stampdata = data[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)] 
                stampmodel = model[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)]
                stampres = res[round(201-dy*10):round(201+dy*10), round(201-dx*10):round(201+dx*10)]
            
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
                if i < 2:
                    fig.text(.5, .52, 'chisq/nu = '+str(chi[counter]), va = 'center', ha = 'center')

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
                if i < 2:
                    fig.text(.5, .08, 'chisq/nu = '+str(chi[2*i+1]), va = 'center', ha = 'center') #when i = 0, we want chi[1], when i = 1, we want chi[3]
                else:
                    fig.text(.5, .05, 'chisq/nu = '+str(chi[i+2]), ha='center')
                    
            counter += 1
        cax = plt.axes([0.9, 0.1, 0.01, 0.8])
        plt.colorbar(cax=cax)
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




