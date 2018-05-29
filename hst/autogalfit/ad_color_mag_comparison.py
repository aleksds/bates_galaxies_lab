# Aleks Diamond-Stanic, May 23, 2018
# code that compares magnitudes, colors, and sizes for two runs of GALFIT
# requires names of two directories as input, current only works for 'independent' runs of galfit
# run ad_color_mag_comparison 20180511-1229_psf_independent 20180511-1229_sersic_independent
import sys
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

dir = os.environ['HSTDIR']

one = sys.argv[1]
two = sys.argv[2]

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']

redshifts = [0.603,0.459,0.712,0.514,0.467,0.451,0.451,0.658,0.608,0.402,0.449,0.728,0.752]
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

#minmax function used to set axis boundaries
def minmax(values):
    axes = np.zeros([len(values)])
    coords = np.zeros(2)
    coords[0] = np.min(values)-0.1
    coords[1] = np.max(values)+0.1
    return coords

# function to add text to plots
def addtext(xvals, yvals):
    med = np.median(yvals)
    mean = np.mean(yvals)
    std = np.std(yvals)
    plt.text(np.min(xvals), np.min(yvals)+(np.max(yvals)-np.min(yvals))*0.6, 'median='+"{:.3f}".format(med))
    plt.text(np.min(xvals), np.min(yvals)+(np.max(yvals)-np.min(yvals))*0.5, 'mean='+"{:.3f}".format(mean))
    plt.text(np.min(xvals), np.min(yvals)+(np.max(yvals)-np.min(yvals))*0.4, 'std='+"{:.3f}".format(std))

model = [str(one),str(two)]
mags = np.zeros([2,12,2]) #order: psf models then sersic models, normal order for galaxies and filters
chi = np.zeros([2,12,2])
sizepix = np.zeros([12,2])
modeltype = ['type1', 'type2']
for m in range(0,len(model)): #loop through each directory
    w = 0
    for file in glob.glob(model[m]+"/*.band"): #loop through each file
        with open(file) as f:
            content = f.readlines()
        print(content[10])
        if content[8][4] == 'V':
            print('yes')
            #file is simultaneous
            #814:
            print(file)
            print(content[47][4:10], m, w)
            mags[m][w][1] = np.float(content[47][4:10])
            chi[m][w][1] = np.float(content[3][13:19])
            if content[44][4:7] == 'psf':
                modeltype[m] = 'psf'
                sizepix[w][0] = 0
            else:
                modeltype[m] = 'sersic'
                sizepix[w][0] = np.float(content[48][4:8])
            #475:
            mags[m][w][0] = np.float(content[47][11:17])
            chi[m][w][0] = np.float(content[3][13:19]) #will be the same as chi[m][w][0]
            if content[44][4:7] == 'psf':
                modeltype[m] = 'psf'
                sizepix[w][1] = 0
            else:
                modeltype[m] = 'sersic'
                sizepix[w][1] = np.float(content[48][4:8])
        else:
            #file is independent
            if np.float(content[10][10:13]) == 814: #note that when i = 0, filter = 814, etc.
                print('814')
                mags[m][w][0] = np.float(content[47][4:10])
                chi[m][w][0] = np.float(content[3][14:19])
                if np.float(content[48][4]) == 0:
                    modeltype[m] = 'psf'
                    sizepix[w][0] = 0
                else:
                    modelype[m] = 'sersic'
                    sizepix[w][0] = np.float(content[48][4:8])
            else:
                print('475')
                mags[m][w][1] = np.float(content[47][4:10])
                chi[m][w][1] = np.float(content[3][14:19])
                if np.float(content[48][4]) == 0:
                   modeltype[m] = 'psf'
                   sizepix[w][1] = 0
                else:
                   modeltype[m] = 'sersic'
                   sizepix[w][1] = np.float(content[48][4:8])
        w = w+1
        print('w =',w)
#color magnitude plot
one_mag_475 = np.zeros(12)
one_mag_814 = np.zeros(12)
one_color = np.zeros(12)
two_mag_475 = np.zeros(12)
two_mag_814 = np.zeros(12)
two_color = np.zeros(12)
for w in range(0,12):
    one_mag_475[w] = mags[0][w][0]
    one_mag_814[w] = mags[0][w][1]
    one_color[w] = mags[0][w][0] - mags[0][w][1]
    two_mag_475[w] = mags[1][w][0]
    two_mag_814[w] = mags[1][w][1]    
    two_color[w] = mags[1][w][0] - mags[1][w][1]

x1 = minmax([one_mag_814, two_mag_814])
y1 = minmax([one_color, two_color])

#name_cm = 'color_vs_mag_'+one+'_'+two+'.pdf'
#with PdfPages(name_cm) as pdf:   

#os.system('open %s &' % name_cm)

#color vs chi-squared plot values
psfyvals = np.zeros(12)
psfxvals_four = np.zeros(12)
sersicyvals = np.zeros(12)
sersicxvals_four = np.zeros(12)
psfxvals_eight = np.zeros(12)
sersicxvals_eight = np.zeros(12)

for w in range(0,12):
    psfxvals_four[w] = chi[0][w][0]
    sersicxvals_four[w] = chi[1][w][0]
    psfyvals[w] = mags[0][w][0] - mags[0][w][1]
    sersicyvals[w] = mags[1][w][0] - mags[1][w][1]
    psfxvals_eight[w] = chi[0][w][1]
    sersicxvals_eight[w] = chi[1][w][1]

x2 = minmax([psfxvals_four, sersicxvals_four, psfxvals_eight, sersicxvals_eight])
y2 = minmax([psfyvals, sersicyvals])


#color vs size for sersic fits
size_color = np.zeros(12)
size_four = np.zeros(12)
size_eight = np.zeros(12)

kpcrad=np.zeros([12,2])
for w in range(0,12):
    arcsecperkpc = cosmo.arcsec_per_kpc_proper(redshifts[w])
    for i in range(0,2):
        kpcrad[w][i] = (0.025*sizepix[w][i])/arcsecperkpc.value

for w in range(0,12):
    size_four[w] = kpcrad[w][0]
    size_color[w] = mags[1][w][0] - mags[1][w][1]
    size_eight[w] = kpcrad[w][1] 

x3 = minmax([size_four, size_eight])
y3 = minmax([size_color])


x4 = minmax([one_mag_475, one_mag_814])
y4 = minmax([one_mag_475-two_mag_475, one_mag_814-two_mag_814])

#comparison plots of mags, colors, chis, and sizes
name_co = 'comparison_'+one+'_'+two+'.pdf'
with PdfPages(name_co) as pdf:
    fig = plt.figure()
    plt.scatter(one_mag_814,one_color, label=one, marker='o', color='orange')
    plt.scatter(two_mag_814,two_color, label=two, marker='^', color='purple')
    plt.xlabel('mag_F814W')
    plt.ylabel('mag_F475W - mag_F814W')
    plt.title('Color Magnitude Plot')
    plt.legend(loc='upper right')
    plt.xlim(x1[0], x1[1])
    plt.ylim(y1[0], y1[1])
    addtext(one_mag_814,one_color)
    addtext(two_mag_814,two_color)
    pdf.savefig()
    plt.close()

    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(one_mag_814[i],one_color[i], label=one, marker='o', color='orange')
        plt.scatter(two_mag_814[i],two_color[i], label=two, marker='^', color='purple')
        plt.xlim(x1[0], x1[1])
        plt.ylim(y1[0], y1[1])
        plt.title(galaxies[i])
        plt.tight_layout()
        
    pdf.savefig()
    plt.close()


    fig = plt.figure()
    
    plt.scatter(psfxvals_four,psfyvals, label=one+' (F475W chi)', marker='o', color='blue')
    plt.scatter(sersicxvals_four,sersicyvals, label=two+' (F475W chi)', marker='^', color='blue')
    plt.scatter(psfxvals_eight,psfyvals, label=one+' (F814W chi)', marker='o', color='green')
    plt.scatter(sersicxvals_eight,sersicyvals, label=two+' (F814W chi)', marker='^', color='green')
    plt.xlim(x2[0],x2[1])
    plt.ylim(y2[0],y2[1])
    plt.xlabel('chi-squared/nu')
    plt.ylabel('mag_F475W - mag_F814W')
    plt.title('Chi Squared vs. Color Plot (with psf and sersic fits)')
    plt.legend(loc='upper right')
    addtext(psfxvals_four,psfyvals)
    pdf.savefig()
    plt.close()

    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(psfxvals_four[i],psfyvals[i], marker='o', color='blue')
        plt.scatter(sersicxvals_four[i],sersicyvals[i], marker='^', color='blue')
        plt.scatter(psfxvals_eight[i],psfyvals[i], marker='o', color='green')
        plt.scatter(sersicxvals_eight[i],sersicyvals[i], marker='^', color='green')
        plt.xlim(x2[0],x2[1])
        plt.ylim(y2[0],y2[1])
        plt.title(galaxies[i])
        plt.tight_layout()
        
    pdf.savefig()
    plt.close()

    fig = plt.figure()
    
    plt.scatter(size_four,size_color, label='F475W size', marker='^', color='blue')
    plt.scatter(size_eight,size_color, label='F814W size', marker='^', color='green')
    plt.xlim(x3[0],x3[1])
    plt.ylim(y3[0],y3[1])
    plt.xlabel('size (kpc)')
    plt.ylabel('magF475W - magF814W')
    plt.title('Color vs Size')
    plt.legend(loc='lower right')
    addtext(size_four,size_color)
    pdf.savefig()
    plt.close()

    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(size_four[i],size_color[i], label='F475W size', marker='^', color='blue')
        plt.scatter(size_eight[i],size_color[i], label='F814W size', marker='^', color='green')
        plt.xlim(x3[0],x3[1])
        plt.ylim(y3[0],y3[1])
        plt.title(galaxies[i])
        plt.tight_layout()
        
    pdf.savefig()
    plt.close()
    
    fig = plt.figure()

    plt.scatter(one_mag_475, one_mag_475-two_mag_475, marker='o', color='blue')
    plt.xlabel('magF475W_'+one)
    plt.ylabel('magF475W - magF475W')    
    plt.title('magF475W comparison')
    plt.xlim(x4[0],x4[1])
    plt.ylim(y4[0],y4[1])
    addtext(one_mag_475, one_mag_475-two_mag_475)
    pdf.savefig()
    plt.close

######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(one_mag_475[i], one_mag_475[i]-two_mag_475[i], marker='o', color='blue')
        plt.xlim(x4[0],x4[1])
        plt.ylim(y4[0],y4[1])
        plt.title(galaxies[i])
        plt.tight_layout()
        
    pdf.savefig()
    plt.close()

    fig = plt.figure()
    
    plt.scatter(one_mag_814, one_mag_814-two_mag_814, marker='o', color='green')
    plt.xlabel('magF814W_'+one)
    plt.ylabel('magF814W - magF814W')    
    plt.title('magF814W comparison')
    plt.xlim(x4[0],x4[1])
    plt.ylim(y4[0],y4[1])
    addtext(one_mag_814, one_mag_814-two_mag_814)
    pdf.savefig()
    plt.close

######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(one_mag_814[i], one_mag_814[i]-two_mag_814[i], marker='o', color='green')
        plt.xlim(x4[0],x4[1])
        plt.ylim(y4[0],y4[1])
        plt.title(galaxies[i])
        plt.tight_layout()

    pdf.savefig()
    plt.close()

    x5 = minmax([one_color])
    y5 = minmax([one_color-two_color])

    fig = plt.figure()

    plt.scatter(one_color, one_color-two_color, marker='o', color='orange')
    plt.xlabel('color_'+one)
    plt.ylabel('difference in color')
    plt.title('color comparison')
    plt.xlim(x5[0],x5[1])
    plt.ylim(y5[0],y5[1])
    addtext(one_color, one_color-two_color)
    pdf.savefig()
    plt.close

######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(one_color[i], one_color[i]-two_color[i], marker='o', color='orange')
        plt.xlim(x5[0],x5[1])
        plt.ylim(y5[0],y5[1])
        plt.title(galaxies[i])
        plt.tight_layout()

    pdf.savefig()
    plt.close()

    x6 = minmax([psfxvals_four, psfxvals_eight])
    y6 = minmax([psfxvals_four/sersicxvals_four, psfxvals_eight/sersicxvals_eight])

    fig = plt.figure()

    plt.scatter(psfxvals_four, psfxvals_four/sersicxvals_four, label='F475W', marker='o', color='blue')
    plt.scatter(psfxvals_eight, psfxvals_eight/sersicxvals_eight, label='F814W', marker='o', color='green')
    plt.xlabel('psf chi square/nu values')
    plt.ylabel('ratio of chi sqr/nu values')
    plt.title('chi sqr/nu comparison')
    plt.legend(loc='lower right')
    plt.xlim(x6[0],x6[1])
    plt.ylim(y6[0],y6[1])
    addtext(psfxvals_eight, psfxvals_eight/sersicxvals_eight)
    pdf.savefig()
    plt.close


######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(psfxvals_four[i], psfxvals_four[i]/sersicxvals_four[i], label='F475W', marker='o', color='blue')
        plt.scatter(psfxvals_eight[i], psfxvals_eight[i]/sersicxvals_eight[i], label='F814W', marker='o', color='green')
        plt.xlim(x6[0],x6[1])
        plt.ylim(y6[0],y6[1])
        plt.title(galaxies[i])
        plt.tight_layout()

    pdf.savefig()
    plt.close()

    x7 = minmax([size_four])
    y7 = minmax([size_four/size_eight])

    fig = plt.figure()

    plt.scatter(size_four, size_four/size_eight, marker='o', color='red')
    plt.xlabel('size_four(kpc)')
    plt.ylabel('ratio of size_four to size_eight')
    plt.title('size comparison')
    plt.xlim(x7[0],x7[1])
    plt.ylim(y7[0],y7[1])
    addtext(size_four, size_four/size_eight)
    pdf.savefig()
    plt.close()

######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(size_four[i], size_four[i]/size_eight[i], marker='o', color='red')
        plt.xlim(x7[0],x7[1])
        plt.ylim(y7[0],y7[1])
        plt.title(galaxies[i])
        plt.tight_layout()

    pdf.savefig()
    plt.close()
    
os.system('open %s &' % name_co)


#goal of this section is to print data/model/residual for each galaxy.
#each galaxy will consist of two pages of plots: first page will be F475W and F814W plots for model one, and second page is same for model two
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
models = [str(one),str(two)]
name_res = 'residuals_'+one+'_'+two+'.pdf'
dx = dy = 25
with PdfPages(name_res) as pdf:
    for i in range(0, len(galaxies)):
        for j in range(0, len(models)):
            fig = plt.figure()
            plt.suptitle(galaxies[i]+' '+models[j]+' model')
            for h in range(0, len(filters)):
                file = glob.glob(models[j]+'/'+galaxies[i]+'_'+filters[h]+'_output.fits')
                multi = fits.open(file[0])
                
                data, data_header = multi[1].data, multi[1].header
                model, res_header = multi[2].data, multi[2].header
                res, res_header = multi[3].data, multi[3].header
                
                stampdata = data[round(101-dy):round(101+dy), round(101-dx):round(101+dx)] 
                stampmodel = model[round(201-dy):round(201+dy), round(201-dx):round(201+dx)]
                stampres = res[round(201-dy):round(201+dy), round(201-dx):round(201+dx)]
                
                if h==0:
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
                    fig.text(.5, .52, 'chisq/nu = '+str(chi[j][i][h]), va = 'center', ha = 'center')

                if h==1: 
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
                    fig.text(.5, .05, 'chisq/nu = '+str(chi[j][i][h]), ha='center')
            cax = plt.axes([0.9, 0.1, 0.01, 0.8])
            plt.colorbar(cax=cax)
            pdf.savefig(dpi=1000)
            plt.close()
os.system('open %s &' % name_res)
                         

