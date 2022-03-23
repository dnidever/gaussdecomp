#!/usr/bin/env python

"""UTILS.PY - Utility functions

"""

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20200122'  # yyyymmdd                                                                                                                           
import os
import numpy as np
import warnings
from scipy import sparse
from scipy.interpolate import interp1d
from dlnpyutils import utils as dln

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


def computenoise(spec):
    # compute the noise in the spectrum

    # Use the "blank" velocity ranges
    if hasattr(spec,'blankrange'):
        if spec.blankrange is not None:
            sp = np.array([])
            for r in spec.blankrange:
                gd,ngd = dln.where((spec.vel>=r[0]) & (spec.vel<=r[1]))
                if ngd>0: sp=np.hstack((sp,spec.flux[gd]))
            if len(sp)>0:
                noise = dln.mad(sp)
                return noise
    # Use values near zero with outlier rejection 
    sig = dln.mad(spec.flux)
    gd,ngd = dln.where(np.abs(spec.flux) <= 3*sig)
    if ngd>20:
        sig = dln.mad(spec.flux[gd])
        gd,ngd = dln.where(np.abs(spec.flux) <= 3*sig)        
        if ngd>20:
            noise = dln.mad(spec.flux[gd])
            return noise
    # Last resort, subtract smoothed version of spectrum
    smspec = dln.smooth(spec.flux,50)
    noise = dln.mad(spec.flux-smspec)
    return noise
    

def gaussmodel(gpars,vel):
    # compute the model spectrum with a bunch of Gaussians and the velocity array
    pass
