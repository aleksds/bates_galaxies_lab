#this will be a code that compares the results of PPDs across different prospector runs and photometry values
#made by Kingdell Valdez

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

#saving all the csv files as tables
oldmass = ascii.read('Old_Mass_PPD_Results.csv')
oldage = ascii.read('Old_Age_PPD_Results.csv')
olddust = ascii.read('Old_Dust_PPD_Results.csv')
newmass = ascii.read('New_Mass_PPD_Results.csv')
newage = ascii.read('New_Age_PPD_Results.csv')
newdust = ascii.read('New_Dust_PPD_Results.csv')
penmass = ascii.read('Penultimate_Mass_PPD.csv')
penage = ascii.read('Penultimate_Age_PPD.csv')
pendust2 = ascii.read('Penultimate_Dust2_PPD.csv')
pendust1 = ascii.read('Penultimate_Dust1_PPD.csv')
totalmass = ascii.read('/Users/kvaldez/github/bates_galaxies_lab/hst/total_PPD_Mass.csv')
totalage = ascii.read('/Users/kvaldez/github/bates_galaxies_lab/hst/total_PPD_Age.csv')
totaldust1 = ascii.read('/Users/kvaldez/github/bates_galaxies_lab/hst/total_PPD_Dust1.csv')
totaldust2 = ascii.read('/Users/kvaldez/github/bates_galaxies_lab/hst/total_PPD_Dust2.csv')

# this gives the names of the columns
oldmass.colnames

###MASS###
#linear old mass
old_ssp_best_mass = np.array(oldmass['SSP Best'])
old_ssp_up_mass = np.array(oldmass['SSP Up Unc'])
old_ssp_lo_mass = np.array(oldmass['SSP Lo Unc'])

#log of old mass variable
old_log_ssp_best_mass = np.log10(old_ssp_best_mass)
old_log_ssp_up_mass = np.log10(old_ssp_up_mass)
old_log_ssp_lo_mass = np.log10(old_ssp_lo_mass)

#linear of new mass
new_ssp_best_mass = np.array(newmass['SSP Best Mass'])
new_ssp_up_mass = np.array(newmass['SSP Up Unc'])
new_ssp_lo_mass = np.array(newmass['SSP Lo Unc'])

#log of new mass
new_log_ssp_best_mass = np.log10(new_ssp_best_mass)
new_log_ssp_up_mass = np.log10(new_ssp_up_mass)
new_log_ssp_lo_mass = np.log10(new_ssp_lo_mass)

#penultimate mass
pen_best_mass = np.array(penmass['best'])
pen_up_mass = np.array(penmass['up unc'])
pen_lo_mass = np.array(penmass['lo unc'])

#total SED mass
net_best_mass = np.array(totalmass['best'])
net_up_mass = np.array(totalmass['up unc'])
net_lo_mass = np.array(totalmass['lo unc'])

###AGE###
#old age variables
old_ssp_best_age = np.array(oldage['SSP Best'])
old_ssp_up_age = np.array(oldage['SSP Up Unc'])
old_ssp_lo_age = np.array(oldage['SSP Lo Unc'])

#new age variables
new_ssp_best_age = np.array(newage['SSP Best Mass'])
new_ssp_up_age = np.array(newage['SSP Up Unc'])
new_ssp_lo_age = np.array(newage['SSP Lo Unc'])

#penultimate age variables
pen_best_age = np.array(penage['best'])
pen_up_age = np.array(penage['up unc'])
pen_lo_age = np.array(penage['lo unc'])

#total SED age
net_best_age = np.array(totalage['best'])
net_up_age = np.array(totalage['up unc'])
net_lo_age = np.array(totalage['lo unc'])

###DUST###
#linear of old dust variables
old_ssp_best_dust = np.array(olddust['SSP Best'])
old_ssp_up_dust = np.array(olddust['SSP Up Unc'])
old_ssp_lo_dust = np.array(olddust['SSP Lo Unc'])

#linear of new age variables
new_ssp_best_dust = np.array(newdust['SSP Best Mass'])
new_ssp_up_dust = np.array(newdust['SSP Up Unc'])
new_ssp_lo_dust = np.array(newdust['SSP Lo Unc'])

#pen dust1
pen_best_dust1 = np.array(pendust1['best'])
pen_up_dust1 = np.array(pendust1['up unc'])
pen_lo_dust1 = np.array(pendust1['lo unc'])

#pen dust2
pen_best_dust2 = np.array(pendust2['best'])
pen_up_dust2 = np.array(pendust2['up unc'])
pen_lo_dust2 = np.array(pendust2['lo unc'])

#total SED dust1
net_best_dust1 = np.array(totaldust1['best'])
net_up_dust1 = np.array(totaldust1['up unc'])
net_lo_dust1 = np.array(totaldust1['lo unc'])

#total SED dust2
net_best_dust2 = np.array(totaldust2['best'])
net_up_dust2 = np.array(totaldust2['up unc'])
net_lo_dust2 = np.array(totaldust2['lo unc'])

filename = 'kv_PPD_comparison.pdf'

with PdfPages(filename) as pdf:
    #mass comparison figure in linear
    fig = plt.figure()

    plt.scatter(old_ssp_best_mass,new_ssp_best_mass)
    plt.errorbar(x=old_ssp_best_mass,y=new_ssp_best_mass,xerr=[old_ssp_lo_mass,old_ssp_up_mass],
                 yerr=[new_ssp_lo_mass,new_ssp_up_mass], ls='none')
    plt.xlabel('Old SSP Mass')
    plt.ylabel('New SSP Mass')
    plt.title('SSP Mass Comparison (Linear)')
    plt.xlim(0, 0.8e11)
    plt.ylim(0, 0.8e11)
    a = np.arange(150)/75.+9
    b = 10**(a)

    print(len(b))
    plt.plot(b,b,linestyle=':',color='red',linewidth=1, label='x = y')

    pdf.savefig()
    plt.close()

    #mass comparison figure in log
    fig = plt.figure()

    log_uperror = old_log_ssp_best_mass-np.log10(old_ssp_best_mass - old_ssp_up_mass)
    log_loerror = old_log_ssp_best_mass-np.log10(old_ssp_best_mass - old_ssp_lo_mass)
    log_uperror2 = new_log_ssp_best_mass-np.log10(abs(new_ssp_best_mass-new_ssp_up_mass))
    log_loerror2 = new_log_ssp_best_mass-np.log10(abs(new_ssp_best_mass - new_ssp_lo_mass))

    plt.scatter(old_log_ssp_best_mass, new_log_ssp_best_mass)
    plt.errorbar(x=old_log_ssp_best_mass, y=new_log_ssp_best_mass, xerr=[log_loerror,log_uperror],
                 yerr=[log_loerror2,log_uperror2], ls='none')
    plt.xlabel('log(Old SSP Mass)')
    plt.ylabel('log(New SSP Mass)')
    plt.title('SSP Mass Comparison (Log)')
    plt.xlim(9, 11.5)
    plt.ylim(9, 11.5)
    a = np.arange(15)
    b = a

    print(len(b))
    plt.plot(a, b, linestyle=':', color='red', linewidth=1, label='x = y')

    pdf.savefig()
    plt.close()

    #age comparison figure
    fig = plt.figure()

    plt.scatter(old_ssp_best_age, new_ssp_best_age)
    plt.errorbar(x=old_ssp_best_age, y=new_ssp_best_age,
                 xerr=[old_ssp_lo_age, old_ssp_up_age],
                 yerr=[new_ssp_lo_age, new_ssp_up_age],
                 ls='none')
    plt.xlabel('Old SSP Age')
    plt.ylabel('New SSP Age')
    plt.title('SSP Age Comparison')
    plt.xlim(0, 0.09)
    plt.ylim(0, 0.09)
    a = np.arange(2)
    b = a
    print(len(b))
    plt.plot(b, b, linestyle=':', color='red', linewidth=1, label='x = y')

    pdf.savefig()
    plt.close()

    # dust comparison figure
    fig = plt.figure()

    plt.scatter(old_ssp_best_dust, new_ssp_best_dust)
    plt.errorbar(x=old_ssp_best_dust, y=new_ssp_best_dust,
                 xerr=[old_ssp_lo_dust,old_ssp_up_dust],
                 yerr=[new_ssp_lo_dust,new_ssp_up_dust],
                 ls='none')
    plt.xlabel('Old SSP Dust')
    plt.ylabel('New SSP Dust')
    plt.title('SSP Dust Comparison')
    plt.xlim(0, 2)
    plt.ylim(0, 2)
    a = np.arange(4)
    b = a

    print(len(b))
    plt.plot(b, b, linestyle=':', color='red', linewidth=1, label='x = y')

    pdf.savefig()
    plt.close()

    #now comparing new prospector to old prospector with new photometry
    # mass comparison figure in log
    fig = plt.figure()

    plt.scatter(new_log_ssp_best_mass, pen_best_mass)
    plt.errorbar(x=new_log_ssp_best_mass, y=pen_best_mass, xerr=[log_loerror2,log_uperror2],
                 yerr=[pen_lo_mass,pen_up_mass], ls='none')
    plt.xlabel('log(New SSP Mass)')
    plt.ylabel('Penultimate SSP Mass')
    plt.title('Penultimate vs New SSP Mass Comparison')
    plt.xlim(9, 11)
    plt.ylim(9, 11)
    a = np.arange(15)
    b = a

    print(len(b))
    plt.plot(a, b, linestyle=':', color='red', linewidth=1, label='x = y')

    pdf.savefig()
    plt.close()

    # this figure shows the relationship between the centroid mass and total galactic mass
    fig = plt.figure()

    net_ratio = pen_best_mass/net_best_mass
    cen_err = (pen_up_mass+pen_lo_mass)/2
    net_err = (net_up_mass+net_lo_mass)/2
    net_ratio_err = net_ratio*np.sqrt((cen_err/pen_best_mass)**2 + (net_err/net_best_mass)**2)
    plt.scatter(net_best_mass,net_ratio)
    plt.errorbar(x=net_best_mass, y=net_ratio, xerr=[net_lo_mass, net_up_mass],
                 yerr=[net_ratio_err,net_ratio_err], ls='none')
    plt.xlabel('Total Galactic Stellar Mass')
    plt.ylabel('Ratio of Centroid Stellar Mass to Total Galactic Stellar Mass')
    plt.title('Centroid-Galactic Stellar Mass Comparison')
    plt.ylim(0,1)

    pdf.savefig()
    plt.close()

os.system('open %s &' % filename)