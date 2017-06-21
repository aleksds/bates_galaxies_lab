# Sophia C W Gottlieb I
# 20170620
#
# This code is to read in the flux values I just made....

# importing relevant code, i hope
import os
import numpy as np
from astropy.io import fits
import glob
import math
datatable = []
# find file directory? Not necessary here. Stop bereting me!
dir = os.environ['HSTDIR']
print('HELLO')
f = open("sg_photontable.txt","r")
for line in f:
    words = line.split()
    datatable.append(words)
    #print(words)
    #print(line)
#f.read()
print('GOODBYE')
