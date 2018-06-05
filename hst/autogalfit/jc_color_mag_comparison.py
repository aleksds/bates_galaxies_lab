# Aleks Diamond-Stanic, May 16, 2018
# goal: start with color_mag_plot.py, make code that compares magnitudes, colors, and sizes for two runs of GALFIT
# requires names of two directories as input, current only works for 'independent' runs of galfit
# run color_mag_comparison.py 20180511-1229_psf_independent 20180511-1229_sersic_independent
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

# check if first directory exists
if os.path.isdir(one):
    print('1st directory exists!')
else:
    print('Uh oh! First directory does not exist!')
	
# check if second directory exists
if os.path.isdir(two):
    print('2nd directory exists!')
else:
    print('Uh oh! Second directory does not exist!')

galaxies = ['J0826', 'J0901', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']
filters = ['F475W','F814W']
dependency = ['dep1', 'dep2']
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

# w = 0;
model = [str(one),str(two)]
mags = np.zeros([2,12,2]) #order: psf models then sersic models, normal order for galaxies and filters
chi = np.zeros([2,12,2])
sizepix = np.zeros([2,12,2])
modeltype = ['type1', 'type2']
for m in range(0,len(model)): #loop through each directory
    for file in glob.glob(model[m]+"/*.band"): #loop through each file
        with open(file) as f:
            content = f.readlines()
        print(content[10])
        print(content[48])
        gal = content[10][3:8]
        w = galaxies.index(gal)
        print(content[10])
        if content[8][4] == 'V':
            print('yes')
            dependency[m] = 'simultaneous'
            #file is simultaneous
            #814:
            mags[m][w][1] = np.float(content[47][4:10])
            chi[m][w][1] = np.float(content[3][14:19])
            test = content[48]
            comma_location = test.find(',')
            size = content[48][4:comma_location]
            if np.float(size) == 0:
            #if np.float(content[48][4]) == 0:
                modeltype[m] = 'psf'
                sizepix[m][w][1] = 0
            else:
                modeltype[m] = 'sersic'
                sizepix[m][w][1] = np.float(content[48][4:8])
            #475:
            mags[m][w][0] = np.float(content[47][11:17])
            chi[m][w][0] = np.float(content[3][14:19]) #will be the same as chi[m][w][0]
            if np.float(size) == 0:
            #if np.float(content[48][4]) == 0:
                modeltype[m] = 'psf'
                sizepix[m][w][0] = 0
            else:
                modeltype[m] = 'sersic'
                sizepix[m][w][0] = np.float(content[48][4:8])
        else:
            #file is independent
            dependency[m] = 'independent'
            if np.float(content[10][10:13]) == 814: #note that when i = 0, filter = 814, etc.
                print('814')
                mags[m][w][1] = np.float(content[47][4:10])
                chi[m][w][1] = np.float(content[3][14:19])
                if np.float(content[48][4]) == 0:
                    modeltype[m] = 'psf'
                    sizepix[m][w][1] = 0
                else:
                    modeltype[m] = 'sersic'
                    sizepix[m][w][1] = np.float(content[48][4:8])
            else:
                print('475')
                mags[m][w][0] = np.float(content[47][4:10])
                chi[m][w][0] = np.float(content[3][14:19])
                if np.float(content[48][4]) == 0:
                   modeltype[m] = 'psf'
                   sizepix[m][w][0] = 0
                else:
                   modeltype[m] = 'sersic'
                   sizepix[m][w][0] = np.float(content[48][4:8])


#color vs size for sersic fits
size_color = np.zeros(12)
size_four_one = np.zeros(12)
size_eight_one = np.zeros(12)
size_four_two = np.zeros(12)
size_eight_two = np.zeros(12)

kpcrad=np.zeros([12,2])
for w in range(0,12):
    arcsecperkpc = cosmo.arcsec_per_kpc_proper(redshifts[w])
    for i in range(0,2):
        kpcrad[w][i] = (0.025*sizepix[1][w][i])/arcsecperkpc.value

for w in range(0,12):
    size_four_one[w] = kpcrad[w][0]
    size_four_two[w] = kpcrad[w][0]
    size_color[w] = mags[1][w][0] - mags[1][w][1]
    size_eight_one[w] = kpcrad[w][1]
    size_eight_two[w] = kpcrad[w][1] 

x3 = minmax([size_four_one, size_eight_one])
x3 = minmax([size_four_two, size_eight_two])
y3 = minmax([size_color])
#if modeltype[0] == 'sersic' or modeltype[1] == 'sersic':
#    name_cs = 'color_v_size_'+one+'_'+two+'.pdf'
#    
#    with PdfPages(name_cs) as pdf:   
#        plt.figure()
#        
#        plt.scatter(size_four,size_color, label='F475W size', marker='^', color='blue')
#        plt.scatter(size_eight,size_color, label='F814W size', marker='^', color='green')
#        plt.xlim(x3[0],x3[1])
#        plt.ylim(y3[0],y3[1])
#        plt.xlabel('size (kpc)')
#        plt.ylabel('magF475W - magF814W')
#        plt.title('Color vs Size')
#        plt.legend(loc='lower right')
#        pdf.savefig()
#        plt.close()
#
#        fig = plt.figure()
#        for i in range(0, len(galaxies)):
#            ax = fig.add_subplot(3,4,i+1)
#            plt.scatter(size_four[i],size_color[i], label='F475W size', marker='^', color='blue')
#            plt.scatter(size_eight[i],size_color[i], label='F814W size', marker='^', color='green')
#            plt.xlim(x3[0],x3[1])
#            plt.ylim(y3[0],y3[1])
#            plt.title(galaxies[i])
#            plt.tight_layout()
#            
#        pdf.savefig()
#        plt.close()
#    os.system('open %s &' % name_cs)
#else:
#        print('No size plots for psf models')
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
#color vs chi-squared plot
oneyvals = np.zeros(12)
onexvals_four = np.zeros(12)
twoyvals = np.zeros(12)
twoxvals_four = np.zeros(12)
onexvals_eight = np.zeros(12)
twoxvals_eight = np.zeros(12)

for w in range(0,12):
    onexvals_four[w] = chi[0][w][0]
    twoxvals_four[w] = chi[1][w][0]
    oneyvals[w] = mags[0][w][0] - mags[0][w][1]
    twoyvals[w] = mags[1][w][0] - mags[1][w][1]
    onexvals_eight[w] = chi[0][w][1]
    twoxvals_eight[w] = chi[1][w][1]

x2 = minmax([onexvals_four, twoxvals_four, onexvals_eight, twoxvals_eight])
y2 = minmax([oneyvals, twoyvals])
x4 = minmax([one_mag_475, one_mag_814])
y4 = minmax([one_mag_475-two_mag_475, one_mag_814-two_mag_814])
#comparison plots of mags, colors, chis, and sizes
name_co = 'comparison_'+one+'_'+two+'.pdf'
with PdfPages(name_co) as pdf:   
    fig = plt.figure()

    plt.scatter(one_mag_475, one_mag_475-two_mag_475, marker='o', color='blue')
    plt.xlabel('magF475W_'+one)
    plt.ylabel('magF475W - magF475W')    
    plt.title('magF475W comparison')
    plt.xlim(x4[0],x4[1])
    plt.ylim(minmax(one_mag_475-two_mag_475))
    #plt.ylim(y4[0],y4[1])
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

    x6 = minmax([onexvals_four, onexvals_eight])
    y6 = minmax([onexvals_four/twoxvals_four, onexvals_eight/twoxvals_eight])

    fig = plt.figure()

    plt.scatter(onexvals_four, onexvals_four/twoxvals_four, label='F475W', marker='o', color='blue')
    plt.scatter(onexvals_eight, onexvals_eight/twoxvals_eight, label='F814W', marker='o', color='green')
    plt.xlabel('chi square/nu values')
    plt.ylabel('ratio of chi sqr/nu values')
    plt.title('chi sqr/nu comparison')
    plt.legend(loc='lower right')
    addtext(onexvals_eight, onexvals_eight/twoxvals_eight)
    plt.xlim(x6[0],x6[1])
    #plt.ylim(minmax(onexvals_four/twoxvals_four))
    #plt.ylim(y6[0],y6[1])
    pdf.savefig()
    plt.close


    ######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(onexvals_four[i], onexvals_four[i]/twoxvals_four[i], label='F475W', marker='o', color='blue')
        plt.scatter(onexvals_eight[i], onexvals_eight[i]/twoxvals_eight[i], label='F814W', marker='o', color='green')
        plt.xlim(x6[0],x6[1])
        #plt.ylim(y6[0],y6[1])
        plt.title(galaxies[i])
        plt.tight_layout()

    pdf.savefig()
    plt.close()

    x7 = minmax([size_four_one])
    x7 = minmax([size_four_two])
    #y7=x7
    if np.min(size_eight_one) == 0.:
    	y7 = minmax([size_four_one/size_eight_one])
    else:
    	y7 = minmax(size_eight_one)
    if np.min(size_eight_one) == 0.:
    	y7 = minmax([size_four_two/size_eight_one])
    else:
    	y7 = minmax(size_eight_one)
    if np.min(size_eight_two) == 0.:
    	y7 = minmax([size_four_one/size_eight_two])
    else:
    	y7 = minmax(size_eight_two)
    if np.min(size_eight_two) == 0.:
    	y7 = minmax([size_four_two/size_eight_two])
    else:
    	y7 = minmax(size_eight_two)

    fig = plt.figure()

    plt.scatter(size_four_one, size_four_one/size_eight_one, marker='o', color='red')
    plt.xlabel('size_four_one(kpc)')
    plt.ylabel('ratio of size_four_one to size_eight_one')
    plt.title('size comparison')
    plt.xlim(x7[0],x7[1])
    #plt.ylim(y7[0],y7[1])
    addtext(size_four_one, size_four_one/size_eight_one)
    pdf.savefig()
    plt.close()

    ######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(size_four_one[i], size_four_one[i]/size_eight_one[i], marker='o', color='red')
        plt.xlim(x7[0],x7[1])
        #plt.ylim(y7[0],y7[1])
        plt.title(galaxies[i])
        plt.tight_layout()     
    pdf.savefig()
    plt.close()

    plt.scatter(size_four_two, size_four_two/size_eight_two, marker='o', color='red')
    plt.xlabel('size_four_two(kpc)')
    plt.ylabel('ratio of size_four_two to size_eight_two')
    plt.title('size comparison')
    plt.xlim(x7[0],x7[1])
    #plt.ylim(y7[0],y7[1])
    addtext(size_four_two, size_four_two/size_eight_two)
    pdf.savefig()
    plt.close()

    ######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(size_four_two[i], size_four_two[i]/size_eight_two[i], marker='o', color='red')
        plt.xlim(x7[0],x7[1])
        #plt.ylim(y7[0],y7[1])
        plt.title(galaxies[i])
        plt.tight_layout()     
    pdf.savefig()
    plt.close()

    plt.scatter(size_four_one, size_four_one/size_eight_two, marker='o', color='red')
    plt.xlabel('size_four_one(kpc)')
    plt.ylabel('ratio of size_four_one to size_eight_two')
    plt.title('size comparison')
    plt.xlim(x7[0],x7[1])
    #plt.ylim(y7[0],y7[1])
    addtext(size_four_one, size_four_two/size_eight_two)
    pdf.savefig()
    plt.close()

    ######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(size_four_one[i], size_four_one[i]/size_eight_two[i], marker='o', color='red')
        plt.xlim(x7[0],x7[1])
        #plt.ylim(y7[0],y7[1])
        plt.title(galaxies[i])
        plt.tight_layout()     
    pdf.savefig()
    plt.close()

    plt.scatter(size_four_two, size_four_two/size_eight_one, marker='o', color='red')
    plt.xlabel('size_four_two(kpc)')
    plt.ylabel('ratio of size_four_two to size_eight_one')
    plt.title('size comparison')
    plt.xlim(x7[0],x7[1])
    #plt.ylim(y7[0],y7[1])
    addtext(size_four_two, size_four_two/size_eight_one)
    pdf.savefig()
    plt.close()

    ######
    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(size_four_two[i], size_four_two[i]/size_eight_one[i], marker='o', color='red')
        plt.xlim(x7[0],x7[1])
        #plt.ylim(y7[0],y7[1])
        plt.title(galaxies[i])
        plt.tight_layout()     
    pdf.savefig()
    plt.close()

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
  
    plt.figure()
    
    plt.scatter(onexvals_four,oneyvals, label=one+' (F475W chi)', marker='o', color='blue')
    plt.scatter(twoxvals_four,twoyvals, label=two+' (F475W chi)', marker='^', color='blue')
    plt.scatter(onexvals_eight,oneyvals, label=one+' (F814W chi)', marker='o', color='green')
    plt.scatter(twoxvals_eight,twoyvals, label=two+' (F814W chi)', marker='^', color='green')
    plt.xlim(x2[0],x2[1])
    plt.ylim(y2[0],y2[1])
    plt.xlabel('chi-squared/nu')
    plt.ylabel('mag_F475W - mag_F814W')
    plt.title('Chi Squared vs. Color Plot (with one and two fits)')
    plt.legend(loc='upper right')
    addtext(onexvals_four,oneyvals)
    pdf.savefig()
    plt.close()

    fig = plt.figure()
    for i in range(0, len(galaxies)):
        ax = fig.add_subplot(3,4,i+1)
        plt.scatter(onexvals_four[i],oneyvals[i], marker='o', color='blue')
        plt.scatter(twoxvals_four[i],twoyvals[i], marker='^', color='blue')
        plt.scatter(onexvals_eight[i],oneyvals[i], marker='o', color='green')
        plt.scatter(twoxvals_eight[i],twoyvals[i], marker='^', color='green')
        plt.xlim(x2[0],x2[1])
        plt.ylim(y2[0],y2[1])
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
                file = glob.glob(models[j]+'/'+galaxies[i]+'_*_output.fits')
                

                if dependency[j] == 'independent':
                    multi = fits.open(file[h])
                    data, data_header = multi[1].data, multi[1].header
                    model, res_header = multi[2].data, multi[2].header
                    res, res_header = multi[3].data, multi[3].header
                    print('independent')
                else:
                    multi = fits.open(file[0])
                    data, data_header = multi[1-h].data, multi[1-h].header
                    model, res_header = multi[3-h].data, multi[3-h].header
                    res, res_header = multi[5-h].data, multi[5-h].header
                    print('not independent', galaxies[i], dependency[j])
                
                stampdata = data#[round(75-dy):round(75+dy), round(75-dx):round(75+dx)] 
                stampmodel = model#[round(75-dy):round(75+dy), round(75-dx):round(75+dx)]
                stampres = res#[round(75-dy):round(75+dy), round(75-dx):round(75+dx)]
                
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
                         

