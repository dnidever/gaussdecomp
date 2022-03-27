#!/usr/bin/env python

"""CUBE.PY - Datacube model class

"""

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@montana.edu>'
__version__ = '20220326'  # yyyymmdd                                                                                                                           
import os
import numpy as np
import warnings
import copy
from scipy import sparse
from scipy.interpolate import interp1d
from astropy.wcs import WCS
from dlnpyutils import utils as dln
from . import utils

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

class Cube:
    # A class for the radio data cube
        
    # Initialize the object
    def __init__(self,data=None,header=None,getfunction=None,veldim=None,vel=None):
        self._getfunction = getfunction
        self._data = data
        if data is not None:       
            self._shape = data.shape
        else:
            self._shape = None
        self._header = header
        if self._header is not None:
            self._wcs = WCS(header)
        if veldim is not None:
            self._veldim = veldim
        if header is not None and veldim is None:
            for i in range(3):
                ctype = header.get('ctype'+str(i+1))
                if 'vel' in ctype.lower() or 'vrad' in ctype.lower() or 'vlsr' in ctype.lower():
                    self._veldim = i+1
        if self._veldim is not None:
            naxis = header.get('naxis'+str(self._veldim+1))
            cdelt = header.get('cdelt'+str(self._veldim+1))
            crval = header.get('crval'+str(self._veldim+1))
            crpix = header.get('crpix'+str(self._veldim+1))
            vel = crval + cdelt * (np.arange(naxis)+1-crpix)
            self._vel = vel
        else:
            self._vel = None
        # X and Y dimensions and sizes
        if self._veldim is not None:
            left = np.arange(3)
            left = np.delete(left,self._veldim)
            self._xdim = left[0]
            naxis1 = header.get('naxis'+str(self._xdim+1))            
            ctype1 = header.get('ctype'+str(self._xdim+1))
            cdelt1 = header.get('cdelt'+str(self._xdim+1))
            crval1 = header.get('crval'+str(self._xdim+1))
            crpix1 = header.get('crpix'+str(self._xdim+1))
            x = crval1 + cdelt1 * (np.arange(naxis1)+1-crpix1)
            self._x = x
            self._nx = len(x)
            self._xtype = ctype1
            self._ydim = left[1]            
            naxis2 = header.get('naxis'+str(self._ydim+1))
            ctype2 = header.get('ctype'+str(self._ydim+1))            
            cdelt2 = header.get('cdelt'+str(self._ydim+1))
            crval2 = header.get('crval'+str(self._ydim+1))
            crpix2 = header.get('crpix'+str(self._ydim+1))
            y = crval2 + cdelt2 * (np.arange(naxis2)+1-crpix2)            
            self._y = y
            self._ny = len(y)
            self._ytype = ctype2
        else:
            self._xdim = None
            self._ydim = None
            self._nx = None
            self._ny = None
            self._x = None
            self._y = None
        
    def __repr__(self):
        out = self.__class__.__name__ + '('
        out += 'N=%d, %.2f < V < %.2f\)n' % \
               (self.n,self.vrange[0],self.vrange[1])
        return out

    def __str__(self):
        out = self.__class__.__name__ + '('
        out += 'N=%d, %.2f < V < %.2f)\n' % \
               (self.n,self.vrange[0],self.vrange[1])
        return out

    def __call__(self,x,y):
        """ Return the spectrum at a given X/Y position."""
        if self._getfunction is not None:
            return self._getfunction(x,y)
        else:
            if self._veldim == 0:
                return np.copy(self._vel), np.copy(self._data[:,x,y])
            elif self._veldim == 1:
                return np.copy(self._vel), np.copy(self._data[x,:,y])                
            elif self._veldim == 2:
                return np.copy(self._vel), np.copy(self._data[x,y,:])                
            else:
                print('not understood')

    def coords(self,x,y):
        """ Get the coordinates for this X/Y position."""
        
        return self._wcs.pixel_to_world(x,y)

    def copy(self):
        """ Create a copy of the cube."""
        return copy.deepcopy(self)
                    
    def write(self,outfile):
        """ Write cube to a file."""
        if self._cube is None:
            print('No data to write out')
            return
        hdu = fits.HDUList()
        hdu.append(fits.PrimaryHDU(self._cube,self._header))
        # add values to header
        if self._veldim is not None:
            hdu[0].header['veldim'] = veldim
        hdu.writeto(outfile,overwrite=True)

    @classmethod
    def read(cls,infile):
        """ Read in a cube from a file."""
        cube,head = fits.getdata(infile,header=True)
        # get information from header?
        return Cube(cube,header=head)
