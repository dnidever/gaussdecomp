#!/usr/bin/env python

"""UTILS.PY - Utility functions

"""

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20200122'  # yyyymmdd                                                                                                                           
import os
import numpy as np
import warnings
from scipy.signal import argrelextrema
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

def garea(par,npow=1,column=False):
    """
    This function returns the area of a gaussian given its parameters 
    If multiple gaussians are given then it returns an array of the 
    areas of the individual gaussians (not combined). 
    The area of a gaussian is A = ht*wid*sqrt(2*pi) 
    """

    if len(par)<3 or par[0]<0:
        raise ValueError('Need at least 3 parameters and height must be non-negative')
     
    npar = len(par) 
    ngauss = npar//3 
    area = (par[0:ngauss*3:3]**npow)*(par[2:ngauss*3:3]/np.sqrt(npow))*np.sqrt(2.0*np.pi) 
     
    if ngauss==1: 
        area = area[0]
     
    # Output the column density 
    # N = 1.83e18 * Total(T_B(K))*dvel = 1.83e18 * area
    if column:
        area *= 1.83e18 
     
    return area 
    

def parcheck(x,y,par):
    """
    This program checks that the gaussian parameters 
    are okay. 
     
    INPUT 
       x      Array of x values 
       y      Array of y values 
       par    Array of gaussian parameters 
     
    OUTPUT 
       flag   Flag to see if some of the gaussian parameters 
              are bad.  flag = 1   Bad 
                        flag = 0   Okay 
     
    Created by David Nidever April 2005 
    """
 
    npts = len(x) 
    ny = len(y) 
    npar =  len(par) 
    ngauss = npar//3 
     
    flag = 0 
     
    # Calculate dX    
    dx = np.median(np.diff(x[0:np.minimum(10,npts)]))
     
    # Checking for bad stuff 
    if npts < 3: 
        flag = 1 
    if npts < npar: 
        flag = 1 
    if npar < 3: 
        flag = 1 
     
    xmin = np.min(x) 
    xmax = np.max(x) 
    ymin = np.min(y) 
    ymax = np.max(y) 
     
    # Checking the parameters 
    for i in range(ngauss): 
        # Checking height 
        if par[3*i] <= 0.01: 
            flag = 1 
        if par[3*i] > np.max(y)*1.1: 
            flag = 1 
     
        # Checking center 
        if par[3*i+1] < np.min(x) or par[3*i+1] > np.max(x): 
            flag = 1 
     
        # Checking width 
        if par[3*i+2] <= 0.5*dx: 
            flag = 1 
        if par[3*i+2] > (np.max(x)-np.min(x))*0.5: 
            flag = 1 
 
    # Checking for Infinite or NAN
    nfinite = np.sum(~np.isfinite(par))
    if nfinite > 0: 
        flag = 1 

    return flag
        
def gremove(par,x,y):
    """
    This program removes gaussians that have bad parameters. 
     
    INPUT 
       par    Array of gaussian parameters 
       x      Array of x values 
       y      Array of y values 
     
    OUTPUT 
       par    Array of gaussian parameters with bad gaussians removed 
     
    Created by David Nidever April 2005 
    """

    nx = len(x) 
    ny = len(y) 
    npar = len(par) 
    ngauss = npar//3 
     
    orig_par = par 
    fpar = None
     
    for i in range(ngauss): 
        ipar = par[3*i:3*i+3] 
     
        # Checking for infinite or NAN 
        nfinite = np.sum(np.isfinite(ipar))
        if nfinite > 0 : 
            flag = 1 
     
        # Getting the area 
        area = garea(ipar) 
     
        # Using PARCHECK.PRO
        flag = 0
        flag = parcheck(x,y,ipar)
    
        # If everything's fine add them to the final parameters
        if flag==0:
            if fpar is not None:
                fpar = np.hstack((fpar,ipar))
            else:
                fpar = ipar

    return fpar

def gremdup(pararr,v):
    """
    This program removes gaussian parameters that 
    are very similar to another one 
     
    INPUT 
       pararr  Array of parameters [ngauss*3] 
       v       Array of velocities 
     
    OUTPUT 
       pararr  Array of parameters with duplicate removed 
     
    Created by David Nidever April 2005 
    """
 
    npar = len(pararr) 
    nv = len(v) 
     
    orig_pararr = pararr 
     
    n = len(pararr)/3 
    vmin = min(v) 
    vmax = max(v) 
    dv = v(1)-v(0) 
    dum = closest(0,v,ind=vindcen) 
     
    # Removing "duplicates" 
    endflag = 0 
    i = 0 
 
    # Checking all N(N-1)/2 possible combinations 
    # First loop 
    while (endflag==0):
        n = len(pararr)//3 
        if i >= n-2: 
            continue
        ipar = pararr[3*i:3*i+3]
        ig = gfunc(v,ipar) 
     
        endflag2 = 0 
        j = i + 1 
 
        # Second loop 
        while (endflag2 == 0): 
            jpar = pararr[3*j:3*j+3]
            jg = gfunc(v,jpar) 
     
            # Use the narrowest range possible
            vlo = np.maximum(vmin,np.minimum(ipar[1]-ipar[2]*3,jpar[1]-jpar[2]*3)) * 0.5   # just to be safe
            vhi = np.minimum(vmax,np.maximum(ipar[1]+ipar[2]*3, jpar[1]+jpar[2]*3)) * 0.5  # just to be safe
            vran = vhi-vlo
            if ~np.isfinite(vran/dv+1.):
                raise ValueError('Velocity range problem')
            ind = np.arange(int(np.round(vran/dv+1)))+vindcen+vlo/dv
     
            # Calculating the similarity 
            sim = np.sum( (ig[ind]-jg[ind])**2 )*dv / ( garea(ipar,npow=2) + garea(jpar,npow=2) ) 
     
            n = len(pararr)/3 
            if j == n-1: 
                endflag2 = 1 
 
            # Similar gaussians, replace by single average 
            if sim < 0.1: 
                n = len(pararr)/3
                torem = np.arange(3)+3*j
                pararr = np.delete(pararr,torem)
                # don't increment
            else: 
                # Increment to next 
                j += 1 
 
        n = len(pararr)/3 
        if i >= n-2: 
            endflag = 1 
 
        i += 1 

    return pararr

def gfunc(x,par,dp,noderiv=True):
    
    """
    This function makes gaussians 
     
    each gaussian has 3 parameters height, center, width 
    one extra for a constant offset 
     
    INPUT 
       x         Array of X values 
       par       Gaussian parameters 
       /noderiv  Don't return the derivative 
     
    OUTPUT 
       th        Theoretical Y values 
       dp        Derivative of gaussian with the parameters (optional) 
     
    When there are problems in this program it returns: 
      th = x*0.+999999. 
      dp = x*0.+999999. 
     
    Created by David Nidever April 2005 
    """
 
    npar = len(par) 
    npts = len(x) 
    ngauss = npar//3 
    th = 0. 
     
    # Looping through gaussians
    if i in range(ngauss):
        ipar = par[i*3:i*3+3]
        g = ipar[0]*np.exp(-0.5*((x-ipar[1])/ipar[2])**2.) 
        th += g
        
    # Adding power terms 
    npow = npar-3*ngauss 
    if npow>0: 
        for i in range(npow): 
            th += par[i+3]*x**float(i) 
 
    # Computing partial derivatives with respect to the parameters 
    if noderiv == False:
        dp = np.zeros((npts,npar),float)
        for i in range(ngauss):
            ipar = par[i*3:i*3+3]
            # Height 
            dp[:,i*3] = np.exp(-0.5*((x-ipar[1])/ipar[2])**2) 
            # Center 
            dp[:,i*3+1] = dp[:,i*3]*ipar[0]*(x-ipar[1])/ipar[2]**2
            # Width
            dp[:,i*3+2] = dp[:,i*3]*ipar[0]*((x-ipar[1])**2)/ipar[2]**3 
        return th,dp
            
    return th 

def printgpar(par,sigpar=None,rms=None,noise=None,chisq=None,niter=None,status=None):
    """
    Printing the gaussian parameters 
     
    Parameters
    ----------
       par      Array of gaussian parameters 
       sigpar   Array of errors of gaussian parameters 
       rms      RMS of gaussian fit 
       noise    Noise level 
       chisq    Chi squared of gaussian fit 
       niter    Number of iterations 
       status   Status of fitting procedure 
     
    Return
    ------
       Results printed to the screen.
     
    Created by David Nidever April 2005 
    """
     
    npar = len(par)
    ngauss = len(par)/3
    if sigpar is None:
        sigpar = np.zeros(npar)
     
    # Printing the results 
    print('----------------------------------------------------------')
    print('#       Height         Center         Width       Area' )
    print('----------------------------------------------------------') 
    for i in range(ngauss): 
        ipar = par[i*3:i*3+3]
        isigpar = sigpar[i*3:i*3+3]
        iarea = garea(ipar)
        print(' %d  %6.2f (%4.2f) %6.2f (%4.2f) %6.2f (%4.2f) %6.2' % (i+1,ipar[0],isigpar[0],ipar[1],isigpar[1],ipar[2],isigpar[2],iarea))
        print('----------------------------------------------------------')
        if niter is not None:
            print('N Iter = ',str(niter))
        if chisq is not None:
            print('ChiSq = %.3f ' % chisq)
        if rms is not None:
            print('RMS = %.3f' % rms)
        if noise is not None:
            print('Noise = %.3f' % noise)
        if status is not None:
            if status<0:
                print('Status = '+str(status)+'   NOT SUCCESSFUL')
        print('')

def gpeak1(spec,top):     
    # Gets peaks 
     
    # getting maxima
    maxind, = argrelextrema(spec, np.greater)

    if len(maxind)>0:
        gd, = np.where(spec[maxind] > top)
        ngd = len(gd)
        if ngd>0:
            maxarr = maxarr[gd] 
        else:
            maxarr = np.array([])
     
    return maxarr

        
def gest(spec,ind,par,nsig=3):
    """
    This program estimates the gaussian parameters 
    peaking at v0. 
    par = [ht,cen,wid] 
     
    Parameters
    ----------
       v     Array of velocities 
       spec  Array of spectrum 
       ind   Index of the maximum to fit 
       nsig  Number of sigmas for gind (nsig=3 by default) 
     
    Returns
    -------
       par   Gaussian parameters of best fit to maximum 
       gind  The indices of the gaussian out to 3 sigma 
     
    If there is some kind of problem this program returns: 
      par = [999999., 999999., 999999.] 
      gind = -1. 
     
    Created by David Nidever April 2005 
    """
 
    gind = -1 
    par = np.zeros(3,float)+999999.    # assume bad unless proven okay 
     
    npar = len(par) 
    npts = len(v) 
    nspec = len(spec) 
     
    # Bad Input Values 
    if (npar == 0) or (npts == 0) or (npts != nspec) or (len(ind) == 0):
        raise ValueError('Bad inputs')
     
    if (ind < 0) or (ind > npts-1): 
        raise ValueError('Index out of range')
     
    cen = v[ind]
     
    # Using derivatives to get the parameters 
    dv = v[1]-v[0]
    dsdv = np.diff(spec)/dv
    vv0 = v-cen 
    r10 = dsdv/spec 
     
    # Finding maxima closest to IND
    mina, = argrelextrema(r10, np.less)
    maxa, = argrelextrema(r10, np.greater)
    if len(mina)==0 or len(maxa)==0:
        return np.array([]),np.array([])
    maxind, = np.where(maxa < ind)
    maxind = maxind[-1]
    lo = maxa[maxind]
    minind, = np.where(mina > ind)
    minind = minind[0]
    hi = mina[minind]
     
    # Finding the slope for r = (d(spec)/dv)/spec 
    # dsdv/spec = -vv0/sig^2. 
    m = (r10[ind+1]-r10[ind])/dv 
    if (-1./m < 0.): 
        return np.array([]),np.array([])
    sig0 = np.sqrt(-1./m) 
    par = np.array([spec[ind],cen,sig0])
     
    # Finding the indices for this gaussian 
    #  3 sigma on each side
    iwid = sig0/dv 
    lo = np.floor(np.maximum(ind-nsig*iwid,0)) 
    hi = np.ceil(np.minimum(ind+nsig*iwid,npts-1))
    ngind = hi-lo+1
     
    # Need at least 6 points, 3 at each side 
    if ngind<6: 
        lo = np.floor(np.maximum(ind-3,0))
        hi = np.ceil(np.minimum(ind+3,npts-1)) 
        ngind = hi-lo+1
     
    # Need at least 6 points, 5 at each side 
    if ngind<6:
        lo = np.floor(np.maximum(ind-5,0))
        hi = np.ceil(np.minimum(ind+5,npts-1))
        ngind = hi-lo+1 

    gind = np.arange(ngdind)+lo

    return par,gind
                     

def setlimits(x,y,par):
    """
    This program sets limits for the gaussian parameters. 
     
    Parameters
    ----------
       x        Array of x values 
       y        Array of y values 
       par      Gaussian parameters 
     
    Returns
    -------
       parinfo  Structure of parameters limits for gfit.pro 
       ngauss   Number of gaussians 
     
    When there are problems in the this program it returns: 
      parinfo = repliceate({limited:[0,0],limits:[0.,0.]},npar) 
      ngauss = 0. 
    Not sure if these are the best things to return when there 
      are problems. 
     
    Created by David Nidever April 2005 
    """

    nx = len(x) 
    ny = len(y) 
    npar = len(par) 
    ngauss = npar//3 
     
    # Bad Input Values 
    if (ngauss == 0) or (npar == 0) or (nx == 0) or (ny == 0) or (nx != ny):
        raise ValueError('Problems with input values')
     
    # Array limits 
    xmin = np.min(x) 
    xmax = np.max(x) 
    ymin = np.min(y) 
    ymax = np.max(y) 

    # Initialize the bounds
    bounds = (np.zeros(npar,float)-np.inf,
              np.zeros(npar,float)+np.fin)
    
    # Setting limits 
    for i in range(ngauss): 
        # Height limits, must be greater than 0, le than height 
        bounds[0][3*i] = 0.01
        bounds[1][3*i] = ymax*2
        # Center limits, must be within velocity bound, must be close to center 
        bounds[0][3*i+1] = xmin
        bounds[1][3*i+1] = xmax
        # Width limits, must be greater than 0, less than half the xscale 
        bounds[0][3*i+2] = 0.5
        bounds[1][3*i+2] = (xmax-xmin)*0.5
     
    # Constant offset, must be great than 0 and less than max height 
    if npar > ngauss*3: 
        bounds[0][ngauss*3] = 0.
        bounds[1][ngauss*3] = ymax+0.1

    return bounds


def gfit(x,y,par,bounds=None,noise=None):
    """
    This program fits gaussians to a spectrum, given initial estimates 
     
    Parameters
    ----------
       x         Array of X values 
       y         Array of Y values 
       par       Initial guess parameters 
       bounds    boundaries for the parameters
       noise     Noise level 
     
    Returns
    ------
       fpar      Final parameters 
       perror    Errors of final parameters 
       rms       RMS of final fit 
       chisq     Chi squared of final fit 
       resid     Residuals of spectrum - final fit 
       noise     Noise level 
       rtime     Run time in seconds 

     
    If there are any problems this program returns: 
        fpar = par 
        perror = par*0.+999999. 
        rms = 999999. 
        chisq = 999999. 
        iter = 0 
        resid = y*0. 
        weights = x*0.+1. 
     
    Written by D. Nidever, April 2005 
    """
                     
    t0 = time.time() 
    npts = len(x) 
    npar = len(par)
     
    # function to use 
    if not keyword_set(func) : 
        func='gdev' 
     
    # Getting the noise level
    if noise is None:
        noise = utils.computenoise(spec)
     
    # Setting the limits
    if bounds is None:
        bounds = utils.setlimits(x,y,par)
    badflag = utils.parcheck(x,y,par)
                     
    # Initial parameters are okay 
    if badflag == 0:
        sigma = np.ones(npts,float)

        initpar = np.copy(par)
        try:
            fpar,cov = curve_fit(gfunc, x, y, p0=initpar, sigma=sigma, bounds=bounds)
            perror = np.sqrt(np.diag(cov))
            dof = npts-npar
            sigpar = perror * np.sqrt(chisq/dof) 
            # total resid 
            result = np.gfunc(x,fpar) 
            resid = y-result 
            rms = np.sqrt(np.sum((resid/sigma)**2.)/(npts-1))  # using Haud's method, pg. 92, eg.6 
        except:
            print('Problems in curve_fit')             
            fpar = par 
            perror = par*0.+999999. 
            rms = 999999. 
            chisq = 999999. 
            resid = np.copy(y)*0. 
        # Problems with the initial parameters 
    else: 
        fpar = par 
        perror = par*0.+999999. 
        rms = 999999. 
        chisq = 999999. 
        resid = y*0. 
     
    rtime = time.time()-t0 

    return fpar,perror,rms,chisq,resid,noise,rtime


def gplot(v,y,par,save=False,outfile=None,ylime=None,xlim=None,
          tit=None,noresid=False,xtit=None,ytit=None,noannot=False):
    """
    Plots the gaussian fit.  You can plot just the gaussians 
    by typing: gplot,v,0,par  (or anything for spec) 
     
    Parameters
    ----------
       v        Array of velocities 
       y        Array of spectral points 
       par      Array of gaussian parameters 
       /color   Plot in color 
       /save    Save to postscript file 
       file     Name of postscript file 
       xrange   Range of x-values 
       ylim   Range of y-values 
       tit      Plot title 
       /noresid Don't overplot the residuals 
     
    Returns
    -------
       None 
     
    Created by David Nidever April 2005 
    """

    nv = len(v) 
    ny = len(y) 
    npar = len(par) 
     
    # In case y is not set you 
    #  can still plot the gaussians by themselves 
    if (ny < nv): 
        y = np.copy(v)*0. 
    ny = len(y) 
     
    # Bad Input Values 
    if (nv == 0) or (ny == 0):
        raise ValueError('Problem with inputs')
     
    if outfile is None:
        outfile = 'gaussfit'
     
    ngauss = len(par)//3
    if (npar > 0): 
        result = gfunc(v,par) 
    else: 
        result = np.copy(y)*0.
     
    # Plotting
    if tit is None:
        tit = 'Gaussian Analysis of HI Spectrum'
    if xtit is None:
        xtit = 'V_{LSR} (km s^{-1})'
    if ytit is None:
        ytit = 'T_B (K)' 
     
    # Setting the ranges 
    if (np.std(y) != 0.):
        if ylim is not None:
            yr = ylim 
        else: 
            yr = [np.min(np.hstack((y,result)))-3.0,np.max(np.hstack((y,result)))*1.1]
                         
    else: 
        if ylim is not None:
            yr = ylim 
        else:
            yr = [np.min(np.hstack((y,result)))-0.1,np.max(np.hstack((y,result)))*1.1] 
    if xlim is not None:
        xr = xlim 
    else: 
        xr = [np.min(v),np.max(v)] 

    plt.plot(xr,[0,0])
    if (ny > 0): 
        plt.plot(v,y,thick=thick)
    if (npar > 0): 
        plt.plot(v,result+0.02,c='r',thick=thick)
             
    # Looping through the gaussians 
    if ngauss > 1: 
        for i in range(ngauss): 
            ipar = par[i*3:i*3+3] 
            d, = np.where((v >= (ipar[1]-3.5*ipar[2])) & (v <= (ipar[1]+3.5*ipar[2])))
            nd = len(d)
            g = gfunc(v[d],ipar) 
            plt.plot(v[d],g,linestyle=i/7,thick=thick)
         
    # Overplotting the residuals 
    if (ny > 0) and (npar > 0) and (not keyword_set(noresid)): 
        resid = y-result 
        if np.std(y) != 0.: 
            plt.plot(v,resid-2.5,thick=thick) 
            plt.plot(v,v*0.-2.5,thick=thick)
            plt.annotate('Residua.s',xy=[xr[0]+0.1*dln.valrange(xr),-1.5])
         
    # Annotations
    if noannot is False:
        plt.annotate('Observed HI',xy=[xr[0]+0.1*dln.valrange(xr),yr[0]+0.9*dln.valrange(yr)])
        plt.annotate('Sum of Gaussians',xy=[xr[0]+0.1*dln.valrange(xr),yr[0]+0.86*dln.valrange(yr)])

    # Save the plot
    if save:
        plt.savefig(outfile)



def gaussmodel(gpars,vel):
    # compute the model spectrum with a bunch of Gaussians and the velocity array
    pass
