#!/usr/bin/env python

"""SPECTRUM.PY - Spectrum model class

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

class Spectrum:
    # A class for the radio spectra
        
    # Initialize the object
    def __init__(self,flux,vel,blankvrange=None)
        self.flux = flux
        self.vel = vel
        self.blankvrange = blankvrange
        self._noise = None
        self._model = None
        self._gpars = None
        
    @property
    def noise(self):
        if self._noise is not None:
            return self._noise
        noise = utils.computenoise(self.flux)
        self._noise = noise
        return noise

    @property
    def model(self):
        if self._model is not None:
            return self._model
        if self._gpars is None:
            return None
        model = utils.gaussmodel(self._gpars,self.vel)
        self._model = model
        return model

    
