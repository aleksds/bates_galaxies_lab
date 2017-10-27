# Kwamae Delva
# 10/27/17
# New Voigt Profile Code Using Linetools and Barak with main goal being to include covering fraction

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import sympy
from sympy import mpmath as mp


def voigt_slow(a, u):
      with mp.workdps(20):
            z = mp.mpc(u, a)
            result = mp.exp(-z*z) * mp.erfc(-1j*z)
      return result.real

    # a : float
    #   Ratio of Lorentzian to Gaussian linewidths (see below).
    # u : array of floats, shape (N,)
    #   The frequency or velocity offsets from the line centre, in units
    #   of the FWHM of the Gaussian broadening (see below).
     # workdps(n)(f):
     #    returns a decorated version of the function f
     #    that sets the decimal precision to n before execution,
     #    and restores the precision afterwards.
     
