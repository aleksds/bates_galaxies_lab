#KingValdez170615 
#Code to test out and run some pydis functions
# if pyDIS isn't in the currnet working directory, add to path
import sys
sys.path.append('/home/kvaldez/github/pydis/')
import numpy as np
# must import, of course
import pydis
    
# reduce and extract the data with the fancy autoreduce script
#wave, flux, err = pydis.autoreduce('objlist.r.txt', 'flatlist.r.txt', 'biaslist.r.txt', 
                                   #'HeNeAr.0005r.fits', stdstar='spec50cal/fiege34.dat')
dir = '/usr/netapp/physics/linux-lab/data/apo/UT170422/'
rbias_file = dir + 'rbias.txt'
bias = pydis.biascombine(rbias_file, trim=True)
