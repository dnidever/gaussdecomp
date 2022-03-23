#!/usr/bin/env python

import os
import time
import numpy as np
from . import utils
from dlnpyutils import utils as dln

def gaussfitter(l,b,par0,sigpar,rms,noise,v,spec,resid,noplot=True,ngauss=None,noprint=False,
                color=False,vmax=None,vmin=None,inpar=None,inv=None,inspec=None,debug=False):
                
    """
    This program tries to do gaussian analysis 
    on HI spectra. 
     
    Parameters
    ----------
       l         Galactic longitude 
       b         Galactic latitude 
       inpar     Initial guess of parameters 
       inv       Use input velocity array 
       inspec    Use input spectrum array 
       /noplot   Don't plot anything 
       /noprint  Don't print anything 
       /debug    Diagnostic printing and plotting 
       /color    Use color 
       vmin      The minimum velocity to use 
       vmax      The maximum velocity to use 
     
    Returns
    -------
       par0      Final parameters 
       v         Array of velocities 
       spec      Spectrum 
       resid     The residuals 
       rms       The rms of the residuals 
     
    Example
    -------

    par,v,spec,resid,rms = gaussfitter(1,1)
     
    Written by D. Nidever, April 2005 
    Translated to python by D. Nidever, March 20222
    """
 
    # Start time 
    gt0 = time.time() 
     
    # Default position 
    if len(l) == 0 or len(b) == 0: 
        l = 295. 
        b = -60. 
     
    # Original printing
    if noprint==False:
        print('Fitting Gaussians to the HI spectrum at (',str(l),',',str(b),')')
     
    # Getting the HI spectrum 
    if len(inv) > 0 and len(inspec) > 0: 
        v = inv 
        spec = inspec 
    else: 
        v,spec = rdhispec(l,b)
     
    # Estimating the noise 
    noise = hinoise(v,spec)
    noise_orig = noise 
     
    # Velocity range 
    if vmin is None:
        vmin = np.min(v)
    if vmax is None:
        vmax = np.max(v)
    gd, = np.where((v >= vmin) & (v <= vmax))
    ngd = len(gd)
    old_v = np.copy(v)
    old_spec = np.copy(spec)
    spec = spec[gd]
    v = v[gd] 
    dv = v[1]-v[0] 
    vmin = np.min(v)
    vmax = np.max(v)
    dum = dln.closest(0,v,ind=vindcen) 
     
    # Initalizing Some Parameters
    if par0 is not None:
        par0 = 0 
    npts = len(v) 
    rms = 999999. 
    rms0 = 999999. 
    count = 0 
    npar0 = 0 
    y = spec 
    endflag = 0 
    addt0 = time.time() 
 
    # TRYING THE FIRST GUESS PARAMETERS 
    if (len(inpar) > 1): 
     
        print('Using First Guess Parameters')
        orig_inpar = inpar 
     
        inpar = utils.gremove(inpar,v,y) # removing bad ones 
        inpar = utils.gremdup(inpar,v)   # removing the duplicates 
     
        # Letting everything float 
        finpar,sigpar,rms,chisq,niter,residuals,status,weights = utils.gfit(v,y,inpar,noise=noise,rt=rt1)
     
        finpar = utils.gremove(finpar,v,y)  # removing bad ones 
        finpar = utils.gremdup(finpar,v)    # removing the duplicates 

        
    # ADDING GAUSSIANS 
 
    # One gaussian added per loop 
    while (flag != 1): 

        if noprint is False:
            print('count = ',str(count))
     
        # Using the initial guess 
        if len(finpar) > 1 : 
            par0 = finpar
        if par0 is not None:
            npar0 = len(par0)
     
        # Getting the Residuals 
        if par0 is not None:
            th = utils.gfunc(v,par0) 
            if th[0] != 999999.: 
                resid = y-th 
            else: 
                resid = y 
        else:
            resid = y 
         
        # Smoothing the data 
        smresid4 = dln.savgolsm(resid, [4,4,2]) 
        smresid16 = dln.savgolsm(resid, [16,16,2]) 
        if (npts > 61): 
            smresid30 = dln.savgolsm(resid, [30,30,2]) 
        else: 
            smresid30 = np.zeros(npts,float)
        if (npts > 101): 
            smresid50 = dln.savgolsm(resid, [50,50,2]) 
        else: 
            smresid50 = np.zeros(npts,float)
        if (npts > 201): 
            smresid100 = dln.savgolsm(resid, [100,100,2]) 
        else: 
            smresid100 = np.zeros(npts,float)
        maxy = dln.max(resid) 

         
        # Setting up some arrays 
        rmsarr = np.zeros(200,float)+999999. 
        pararr = np.zeros((200,3),float)+999999. 
        gcount = 0 
         
        # Looping through the different smoothed spectra 
        for s in range(5):
         
            # Getting the right smoothed spectrum 
            if s==0:
                smresid = smresid4
            elif s==1:
                smresid = smresid16
            elif s==2:
                smresid = smresid30
            elif s==3:
                smresid = smresid50
            elif s==4:
                smresid = smresid100 
         
            # Getting maxima 
            maxarr = utils.gpeak1(smresid,np.maximum(5*noise,0.5*maxy))
            if len(maxarr)==0:
                maxarr = utils.gpeak1(smresid,np.maximum(5*noise,0.1*maxy))
            if len(maxarr)==0:
                maxarr = utils.gpeak1(smresid,5*noise)
            if len(maxarr)==0:
                maxarr = utils.gpeak1(smresid,3*noise)
            if len(maxarr)==0:
                maxarr = utils.gpeak1(smresid,2*noise)
            ngd = len(maxarr) 
         
            # If there are any peaks check them out 
            if ngd > 0:
                # Looping through the peaks and fitting them 
                # to the proper residual spectrum 
                # only one gaussian, plus a line 
                for i in range(ngd): 
                    maxi = maxarr[i]
                    # Getting guess and fitting it 
                    par2,gind = utils.gest(v,smresid,maxi,par2)
                    # Only do it if there is a good range 
                    if gind[0] != -1: 
                        par2s = [par2,0.01,0.01] 
                        fpar2s,sigpar2,rms2,chisq2,niter,residuals,status,weights = utils.gfit(v[gind],smresid[gind],par2s,noise=noise,rt=rt2,func='gdev1')
                        fpar2 = fpar2s[0:3]
                        rms = np.std(resid-gfunc(v,fpar2))
                    else: 
                        rms = 999999. 
                        fpar2 = par2*0.+999999. 
             
                    # Adding to the rms and par arrays 
                    rmsarr[gcount] = rms 
                    pararr[gcount,:] = fpar2 
             
                    gcount += 1
 
        # Only taking the good rms 
        gdrms , = np.where(rmsarr != 999999.) 
        ngdrms = len(gdrms)
        
        # No good ones 
        if ngdrms == 0: 
            if (count == 0): 
                flag = 1 
                par0 = -1 
                continue
            else: 
                flag = 1 
                continue 
        rmsarr = rmsarr[gdrms]
        pararr = pararr[gdrms,:]
        pararr = pararr.flatten()

        # Removing bad gaussians 
        pararr = utils.gremove(pararr,v,y)
 
        # Removing the duplicates 
        pararr = utils.gremdup(pararr,v)
 
        # No good ones 
        if len(pararr) < 2: 
            flag = 1 
            continue
 
 
        # Setting up some arrays 
        n = len(pararr)/3 
        rmsarr2 = np.zeros(n,float)
        if par0 is not None:
            pararr2 = np.zeros((n,len(par0)+3),float)
        else:
            pararr2 = np.zeros((n,3),float)
            sigpararr2 = np.copy(pararr2) 
 
        if noprint is False:
            print(str(n,2),' gaussian candidates')
 
 
        # Looping through the candidates that are left 
        #  allowing everything to float 
        for i in range(n): 
            ipar = pararr[3*i:3*i+3] 
            if par0 is not None:
                par = np.hstack((par0,ipar))
            else:
                par = np.copy(ipar)
 
            # Letting everything float 
            fpar,sigpar,rms,chisq,niter,residuals,status,weights = utils.gfit(v,y,par,noise=noise,rt=rt5)
            npar = len(fpar)
            rmsarr2[i] = rms 
            pararr2[i,:] = fpar 
            sigpararr2[i,:] = sigpar 
 
            # Giagnostic plotting and printing
            if debug:
                utils.gplot(v,y,fpar)
                utils.printgpar(fpar,sigpar,len(fpar)/3,rms,noise,chisq,niter,status)
                sleep(1.5)
 
        # Adding the gaussian with the lowest rms 
        gdi = where((rmsarr2 == np.min(rmsarr2)) & (rmsarr2 != 999999))[0][0]
        ngdi = len(gdi)
        if ngdi > 0: 
            par0 = pararr2[gdi,:]
            npar0 = len(par0) 
            sigpar = sigpararr2[gdi,:]
            rms = rmsarr2[gdi]
            sigpar00 = sigpar 
            rms00 = rms 
 
        # Plot what we've got so far 
        #if not keyword_set(noplot) then gplot,v,spec,par0 
 
        # Printing the parameters 
        #ngauss = n_elements(par0)/3
        if noprint is False and rms != 999999: 
            utils.printgpar(par0,sigpar,ngauss,rms,noise,chisq,niter)
 
        # Ending Criteria 
        drms = (rms0-rms)/rms0 
        if (rms <= noise) : 
            flag = 1 
        if (rms >= rms0) : 
            flag = 1 
        if (drms < 0.02) : 
            flag = 1 
 
        # Saving the rms 
        rms0 = rms 
 
        count += 1

      
    if noprint is False:
        print(' Addition of Gaussians ',time.time()-addt0,' sec')
 
    # Saving the results of the fit 
    sigpar00 = sigpar 
    rms00 = rms 
 
    # Plot what we've got so far
    if noplot is False:
        utils.gplot(v,spec,par0)
 
    # Printing the parameters 
    ngauss = len(par0)/3
    if noprint is False:
        utils.printgpar(par0,sigpar00,ngauss,rms,noise)
 
    # REMOVING GAUSSIANS 
 
    # Removing bad gaussians 
    par0 = utils.gremove(par0,v,y)
    
    # Removing the duplicates 
    par0 = utils.gremdup(par0,v)
 
    # None left - End 
    if par0(0) == -1 : 
        return None
 
    # Sorting the gaussians by area 
    orig_par0 = par0 
    ngauss = len(par0)/3 
    par0 = np.argsort(par0) 
 
    # See if removing the smallest (in area) gaussian changes much 
 
    count = 0 
    ngauss = len(par0)/3 
 
    orig_par0 = par0 
    nrefit = np.ceil(ngauss/2) 
 
    weights = np.ones(npts,float)
    rms0 = np.sqrt(np.sum(weights*(spec-gfunc(v,par0))**2.)/(npts-1)) 
 
    # Looping through the smallest 1/2 
    for i in range(nrefit): 
        ig = ngauss-i-1 
        tpar = par0 
 
        # Remove the gaussian in question
        tpar = np.delete(tpar,np.arange(3)+ig*3)
 
        # Finding the best fit 
        fpar,sigpar,rms,chisq,iter,residuals,status,weights = utils.gfit(v,y,tpar,weights=weights,noise=noise,rt=rt5)
 
        drms = rms-rms0 
        # remove the gaussian 
        if (drms/rms0 < 0.02): 
            par0 = fpar 
            rms0 = np.sqrt(np.sum(weights*(spec-utils.gfunc(v,par0))**2)/(npts-1))
            # saving the results of the fit 
            #  without this gaussian 
            sigpar00 = sigpar 
            rms00 = rms 
            chisq00 = chisq 
            iter00 = iter 
            status00 = status 
 
 
    # Removing bad gaussians 
    par0 = utils.gremove(par0,v,y)
 
    # Run it one last time for the final values 
    fpar,sigpar,rms,chisq,niter,residuals,status = utils.gfit(v,y,par0,weights=weights,noise=noise,rt=rt5)
    par0 = fpar 
    ngauss = len(par0)/3 
 
    # Final Results 
    ngauss = len(par0)/3
    if noplot is False:
        utils.gplot(v,spec,par0)
    if noprint is False:
        utils.printgpar(par0,sigpar,ngauss,rms,noise,chisq)
 
    # Printing the total time
    if noprint is False:
        print('Total time = ',str(time.time()-gt0))
 
    #BOMB: 
 
    # No gaussians found 
    if par0[0] == -1: 
        sigpar = par0*0.-1 
        rms = np.std(spec) 
              
    spec = old_spec 
    v = old_v 
 
    return par0
