"""
Aleks Diamond-Stanic 20171013
code to test out and run some pydis functions

can install pydis with the following:
    pip install git+https://github.com/jradavenport/pydis.git
"""

import numpy as np
import os
import pydis

"""
note that pydis code that has lines using np.loadtxt like the following

    np.loadtxt(biaslist,dtype='string')

that should be replaced with something like

    files = np.genfromtxt(biaslist,dtype='str')
    
 """
    
dir = os.environ['APODIR']+'UT170422/'
#rbias_file = dir + 'rbias.txt'
#bias = pydis.biascombine(rbias_file), trim=True)
#bias = pydis.biascombine(rbias_file, trim=False)


objlist = dir + 'objlist.r.txt'
flatlist = dir + 'flatlist.r.txt'
biaslist = dir + 'biaslist.r.txt'
arc = dir + 'arc0905.0030r.fits'

pydis.autoreduce(objlist, flatlist, biaslist, arc, stdstar='spec50cal/fiege34.dat', trim=False)
