# Kwamae Delva
# Early March, 2018
#
# https://github.com/linetools/linetools/blob/master/docs/examples/Voigt_examples.ipynb

import sys
sys.path.append('/Users/kdelva/github/linetools/')

# suppress warnings for these examples
import warnings
warnings.filterwarnings('ignore')


### import
import astropy.units as u
import linetools
from linetools.spectralline import AbsLine, SpectralLine
from linetools.spectra import io as lsio
from linetools.analysis import voigt as lav
from linetools import spectralline as ltsp
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.isgm import abscomponent as lt_abscomp
from linetools.lists.linelist import LineList
import numpy as np
from pkg_resources import resource_filename
from scipy import integrate
from matplotlib.ticker import AutoMinorLocator
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
import os

# plots
from matplotlib import pyplot as plt
import pylab
pylab.rcParams['figure.figsize'] = (8.0, 6.0)
import matplotlib
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)


# define the data directory
dir = os.environ['HIRESDIR']


## from hires_plot_nh.py

gal_info = ascii.read(dir+'gal_info.txt')
gal = gal_info['gal']
zem = gal_info['zem']
vcen = gal_info['vcen']  # Could possible be abslin.attrib['b']

xspec = XSpectrum1D.from_file(resource_filename('linetools','/spectra/tests/files/UM184_nF.fits'))


# wavelengths of relevant absorption lines
mgii2803 = 2803.5314853
mgii2796 = 2796.3542699
mgi2852 = 2852.96328

feii2600 = 2600.1724835
feii2586 = 2586.6495659
feii2382 = 2382.7641781
feii2374 = 2374.4603294
feii2344 = 2344.2129601

#oscillator strengths (fosc)
f2803 = 0.3058
f2796 = 0.6155
f2852 = 1.83

f2600 = 0.2394
f2586 = 0.069126
f2382 = 0.320
f2374 = 0.0313
f2344 = 0.1142

# array of names, lines of interest, oscillator strength:
lines = [mgii2796, mgii2803, mgi2852, feii2586, feii2600, feii2374, feii2382, feii2344]
names = ['MgII 2796', 'MgII 2803', 'MgI 2852', 'Fe II 2586', 'Fe II 2600', 'Fe II 2374', 'Fe II 2382', 'Fe II 2344']
fosc = [f2796, f2803, f2852, f2586, f2600, f2374, f2382, f2344]

minorLocator = AutoMinorLocator()
filename = 'Mg_Tau_Flux_Column_Comparison.pdf'
with PdfPages(filename) as pdf:
    # for h in range(0, len(gal)):
    for h in range(0, 1):
        datafile = dir+gal[h]+'/'+gal[h]+'_stitched_v1.txt'
        data = ascii.read(datafile)
        wave = data['wv'] 
        flux = data['norm']
        fx = data['fx']
        var = data['var']


        # Define a plotting function
        def plt_line(spec):
            #plt.clf()
            #plt.figure(dpi=1200)
            plt.plot(spec.wavelength.value, spec.flux.value, 'k-', drawstyle='steps-mid', lw=1.5)
            plt.xlim(2500, 13000)
            plt.ylabel('Normalized Flux', fontsize=20.)
            plt.xlabel('Wavelength', fontsize=20.)
            #ax = plt.gca()
            #ax.xaxis.set_major_locator(plt.MultipleLocator(2.))
            #ax.get_xaxis().get_major_formatter().set_useOffset(False)
            plt.ylim(0., 1.1)
            plt.show()
            plt.close()

        #Defines linetools plotting fucntion
        def linetools(wv, ln, N, b, z):
            ln.analy['spec'] = xspec
            vmodel = ln.generate_voigt()
            plt_line(vmodel)

        def column(vel, col_dens):
            cor_col = np.array([])
            for i in range(0, len(vel)):
                if (vel[i] >= -3000 and vel[i] <= 500):
                    cor_col = np.append(cor_col, col_dens[i])
            return cor_col

        ism = LineList('ISM')
        
        #Velocity Calculations

        vel_kms = np.zeros([len(lines),len(wave)])
        # define the velocity scale [km / s]
        for i in range(0, len(lines)):
            vel_kms[i] = ((wave-lines[i]*(1+zem[h]))/(lines[i]*(1+zem[h]))) * 3E5

        #Tau Calculations

        tau = np.zeros([len(flux)])
        # loop over each spectral line-tau is an 8 by 50515 array, with 50515 values of tau for each spectral line

        ## Error Bar Calculations

        sigma = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma[i] = flux[i] * np.sqrt(var[i])/fx[i]
        blah = np.log(1/flux)
        tau = blah

        #Error in Column Density Calculations

        sigma_coldens = np.zeros([len(flux)])
        for i in range(0, len(flux)):
            sigma_coldens[i] = (sigma[i]/flux[i])

#COLUMN Info        
        col_2796 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2803 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        col_2852 = column(vel_kms[2],tau/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))

        col_2586 = column(vel_kms[0],tau/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        col_2600 = column(vel_kms[1],tau/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        col_2374 = column(vel_kms[2],tau/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))
        col_2382 = column(vel_kms[3],tau/(2.654E-15*fosc[3]**2 *(wave/(1+zem[h]))))
        col_2344 = column(vel_kms[4],tau/(2.654E-15*fosc[4]**2 *(wave/(1+zem[h]))))

        

        vel_2796 = np.linspace(-3000,500, num = len(col_2796), endpoint = 'True')
        vel_2803 = np.linspace(-3000,500, num = len(col_2803), endpoint = 'True')
        vel_2852 = np.linspace(-3000,500, num = len(col_2852), endpoint = 'True')

        vel_2586 = np.linspace(-3000,500, num = len(col_2586), endpoint = 'True')
        vel_2600 = np.linspace(-3000,500, num = len(col_2600), endpoint = 'True')
        vel_2374 = np.linspace(-3000,500, num = len(col_2374), endpoint = 'True')
        vel_2382 = np.linspace(-3000,500, num = len(col_2382), endpoint = 'True')
        vel_2344 = np.linspace(-3000,500, num = len(col_2344), endpoint = 'True')
    
#Error Limits        
        sigma_coldens2796 = column(vel_kms[0], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2803 = column(vel_kms[1], sigma_coldens/(2.654E-15*fosc[1]**2 *(wave/(1+zem[h]))))
        sigma_coldens2852 = column(vel_kms[2], sigma_coldens/(2.654E-15*fosc[2]**2 *(wave/(1+zem[h]))))

        sigma_coldens2586 = column(vel_kms[3], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2600 = column(vel_kms[4], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2374 = column(vel_kms[5], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2382 = column(vel_kms[6], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        sigma_coldens2344 = column(vel_kms[7], sigma_coldens/(2.654E-15*fosc[0]**2 *(wave/(1+zem[h]))))
        

        # print('This is max coldens for MgII2796', np.max(sigma_coldens2796)))
        # print('This is max coldens for MgII2803', np.max(sigma_coldens2803)))
        # print('This is max coldens for MgI2852', np.max(sigma_coldens2852)))

        # print('This is max coldens for MgII2586', np.max(sigma_coldens2586)))
        # print('This is max coldens for MgII2600', np.max(sigma_coldens2600)))
        # print('This is max coldens for MgII2374', np.max(sigma_coldens2374)))
        # print('This is max coldens for MgII2382', np.max(sigma_coldens2382)))
        # print('This is max coldens for MgII2344', np.max(sigma_coldens2344)))



        wave = np.linspace(2500, 13000,100)*u.AA   ### Need to modify wave so that it focuses on +/- 100 or 400 nm above and below wavelength of flux min, depending on if multi component
        
        #Generates Voigt Profile

        mgi2852 = 2852.96328 * u.AA
        abslin = AbsLine(2852.96328 * u.AA,z=zem[h], linelist=ism)
        abslin.attrib['N'] = np.max(sigma_coldens2852)/u.cm**2  # log N
        abslin.attrib['b'] = 25.*u.km/u.s   ### Still gotta figure out how to get fwhm
        abslin.attrib['z'] = zem[h]
        # abslin.analy['spec'] = xspec
        # vmodel = abslin.generate_voigt()
        # plt_line(vmodel)  # Plot

        def linetools(wv, ln, N, b, z):   ##### LETS GOOOO MY FUNCTION WORKS! #####
            ln.analy['spec'] = xspec
            vmodel = ln.generate_voigt()
            plt_line(vmodel)

        linetools(mgi2852, abslin, abslin.attrib['N'], abslin.attrib['b'], abslin.attrib['z'])



        # # Mg 2796 2803 Multi Component
        # abslin1 = AbsLine(2803.5314853 * u.AA,z=zem[h], linelist=ism)
        # abslin1.attrib['N'] = np.max(sigma_coldens2803)/u.cm**2  # log N
        # abslin1.attrib['b'] = 25.*u.km/u.s
        # abslin.attrib['z'] = zem[h]
        # abslin1.analy['spec'] = xspec
        # vmodel = abslin.generate_voigt()

        # abslin2 = AbsLine(2796.3542699 * u.AA,z=zem[h], linelist=ism)
        # abslin2.attrib['N'] = np.max(sigma_coldens2796)/u.cm**2  # log N
        # abslin2.attrib['b'] = 25.*u.km/u.s
        # abslin.attrib['z'] = zem[h]
        # abslin2.analy['spec'] = xspec
        # vmodel2 = lav.voigt_from_abslines(wave,[abslin1,abslin2])
        # plt_line(vmodel2)  # Plot

        

        # # Fe 2586 and 2600 Multi Comp
        
        # abslin3 = AbsLine(2600.1724835 * u.AA,z=zem[h], linelist=ism)
        # abslin3.attrib['N'] = np.max(sigma_coldens2600)/u.cm**2  # log N
        # abslin3.attrib['b'] = 25.*u.km/u.s
        # abslin3.attrib['z'] = zem[h]
        # abslin3.analy['spec'] = xspec
        # vmodel3 = abslin.generate_voigt()

        # abslin4 = AbsLine(2586.6495659 * u.AA,z=zem[h], linelist=ism)
        # abslin4.attrib['N'] = np.max(sigma_coldens2586)/u.cm**2  # log N
        # abslin4.attrib['b'] = 25.*u.km/u.s
        # abslin.attrib['z'] = zem[h]
        # abslin4.analy['spec'] = xspec
        # vmodel4 = lav.voigt_from_abslines(wave,[abslin3,abslin4])
        # plt_line(vmodel4)  # Plot



        
        # # Fe 2374 and 2382 Multi Comp
        # abslin5 = AbsLine(2382.7641781 * u.AA,z=zem[h], linelist=ism)
        # abslin5.attrib['N'] = np.max(sigma_coldens2382)/u.cm**2  # log N
        # abslin5.attrib['b'] = 25.*u.km/u.s
        # abslin.attrib['z'] = zem[h]
        # abslin5.analy['spec'] = xspec
        # vmodel5 = abslin.generate_voigt()

        # abslin6 = AbsLine(2374.4603294 * u.AA,z=zem[h], linelist=ism)
        # abslin6.attrib['N'] = np.max(sigma_coldens2374)/u.cm**2  # log N
        # abslin6.attrib['b'] = 25.*u.km/u.s
        # abslin.attrib['z'] = zem[h]
        # abslin6.analy['spec'] = xspec
        # vmodel6 = lav.voigt_from_abslines(wave,[abslin5,abslin6])
        # plt_line(vmodel6)  # Plot
        


        
        # #Fe 2344 Comp
        # abslin7 = AbsLine(2344.2129601 * u.AA,z=zem[h], linelist=ism)
        # abslin7.attrib['N'] = np.max(sigma_coldens2344)/u.cm**2  # log N
        # abslin7.attrib['b'] = 25.*u.km/u.s
        # abslin.attrib['z'] = zem[h]
        # abslin7.analy['spec'] = xspec
        # vmodel7 = abslin.generate_voigt()
        # plt_line(vmodel7)  # Plot

        

        # Next steps, figuring out how to get b w/ kd_linetools_b_solver.py
        #Modifying wave based on each component's min flux
                    
