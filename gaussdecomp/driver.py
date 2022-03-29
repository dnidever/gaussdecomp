#!/usr/bin/env python

import os
import time
import numpy as np
from datetime import datetime
from dlnpyutils import utils as dln
from . import utils,fitter
from .cube import Cube

# Tracking lists
BTRACK = {'data':[],'count':0,'x':np.zeros(n,int)-1,'y':np.zeros(n,int)-1,'ngauss':np.zeros(n,int)-1}
n = 100000
GSTRUC = {'data':[],'count':0,'x':np.zeros(n,int)-1,'y':np.zeros(n,int)-1,'ngauss':np.zeros(n,int)-1}

def gstruc_add(tstr):
    """ Add to the GSTRUC tracking structure."""

    # data:   large list of data, one element per position (with padding)
    # count:  the number of elements in DATA we're using currently, also
    #           the index of the next one to start with.
    # x/y:    the x/y position for the gaussians
    # ngauss: the number of gaussians

    global BTRACK, GSTRUC
    
    # Add new elements
    if len(tstr)+GSTRUC['count'] > len(GSTRUC['x']):
        print('Adding more elements to GSTRUC')
        for n in ['x','y','ngauss']:
            GSTRUC[n] = np.hstack((GSTRUC[n],np.zeros(100000,int)-1))
    # Stuff in the new data
    count = GSTRUC['count']
    GSTRUC['data'] += [tstr]
    GSTRUC['x'][count] = tstr['x']
    GSTRUC['y'][count] = tstr['y']
    if tstr['par'] is not None:
        GSTRUC['ngauss'][count] = tstr['par']//3
    else:
        GSTRUC['ngauss'][count] = 0
    GSTRUC['count'] += 1

def gstruc_replace(tstr):
    """ Replace an old decomposition with a newer and better one."""
    # Double-check that the positions match
    ind, = np.where((GSTRUC['x']==tstr['x']) & (GSTRUC['y']==tstr['y']))
    if len(ind)==0:
        raise ValueError('No position for (%d,%d) found in GSTRUC' % (tstr['x'],tstr['y']))
    ind = ind[0]
    GSTRUC['data'][ind] = [tstr]
    GSTRUC['ngauss'] = len(tstr['par'])
    
    
def btrack_add(track):
    """ Add to the BTRACK tracking structure."""

    # data:   large list of data, one element per position (with padding)
    # count:  the number of elements in DATA we're using currently, also
    #           the index of the next one to start with.
    # x/y:    the x/y position for the gaussians
    # ngauss: the number of gaussians

    global BTRACK, GSTRUC
    
    # Add new elements
    if BTRACK['count']+1 > len(BTRACK['x']):
        print('Adding more elements to BTRACK')
        for n in ['x','y','ngauss']:
            BTRACK[n] = np.hstack((BTRACK[n],np.zeros(100000,int)-1))
    # Stuff in the new data
    count = BTRACK['count']
    BTRACK['data'] += [track]
    BTRACK['x'][count] = track['x']
    BTRACK['y'][count] = track['y']
    if track['par'] is not None:
        BTRACK['ngauss'][count] = track['par']//3
    else:
        BTRACK['ngauss'][count] = 0
    BTRACK['count'] += 1

    
def gincrement(x,y,xr,yr,xsgn=1,ysgn=1,nstep=1,p2=False):
    """
    This program increments the position 
     
    Parameters
    ----------
    x : int
      X of current position 
    y : int
      Y of current position 
    xr : list/array
      Two element array of x limits.
    yr  : list/array
      Two element array of y limits.
    xsgn : int, optional
      Sign of x increment.  Default is 1.
    ysgn : int, optional
      Sign of y increment, Default is 1.
    nstep : int, optional
      Number of steps to increment.  Default is 1.
    p2 : boolean, optional
      Increment in y rather than in x. Default is False.
     
    Returns
    -------
    newx : int,
      X of new position 
    newy : int
      Y of new position 
     
    When this program has problems it returns: 
    newx = None
    newy = None
     
    Created by David Nidever April 2005 
    Translated to python by D. Nidever, March 2022
    """
 
    step = 1
     
    # Bad until proven okay 
    newx = None
    newy = None
     
    # X and y limits 
    x0 = xr[0]
    x1 = xr[1] 
    y0 = yr[0] 
    y1 = yr[1] 
     
    # Are we in the box? 
    if (x < x0) or (x > x1) or (y < y0) or (y > y1):
        return None,None

    # Wrong signs 
    if np.abs(xsgn) != 1: 
        return None,None
    if np.abs(ysgn) != 1:
        return None,None        
     
    # figuring out the case 
    #bit = ((ysgn+1.)*0.5) + 2.*((xsgn+1.)*0.5) 
    #bit = int(bit) 
    # 
    # bit (x,y) 
    # 0 - (-1,-1) 
    # 1 - (-1,+1) 
    # 2 - (+1,-1) 
    # 3 - (+1,+1) 

    tx = x
    ty = y
     
    # Looping through all the steps 
    for i in range(nstep): 
        
        # p2, incrementing vertically 
        if p2:
            # UP, at the end 
            if (ysgn == 1) and (ty == y1):
                return None,None
            # DOWN, at the end 
            if (ysgn == -1) and (ty == y0): 
                return None,None
            # Not at end, normal increment 
            newx = tx 
            newy = ty + ysgn * step 
            
        # Incrementing Sideways 
        else: 
            # RIGHT, xsgn = +1 
            if (xsgn == 1): 
                # UP, the very end 
                if (ysgn == 1) and (tx == x1) and (ty == y1): 
                    return None,None
                # DOWN, the very end 
                if (ysgn == -1) and (tx == x1) and (ty == y0): 
                    return None,None
                # At end of x, increment y 
                if (tx == x1): 
                    newx = x0 
                    newy = ty + ysgn * step 
                # Normal increment 
                if (tx != x1): 
                    newx = tx + xsgn * step 
                    newy = ty 
         
            # LEFT, xsgn = -1 
            if (xsgn == -1): 
                # UP, the very end 
                if (ysgn == 1) and (tx == x0) and (ty == y1): 
                    return None,None
                # DOWN, the very end 
                if (ysgn == -1) and (tx == x0) and (ty == y0): 
                    return None,None
                # At end of x, increment y 
                if (tx == x0): 
                    newx = x1 
                    newy = ty + ysgn * step 
                # Normal increment 
                if (tx != x0): 
                    newx = tx + xsgn * step 
                    newy = ty 
     
        # In case we're looping 
        tx = newx 
        ty = newy 
     
    # Final answer 
    newx = tx 
    newy = ty 

    return newx,newy


def gredo(x,y,guessx,guessy,guesspar):
    """
    This function checks wether we can redo this location again.
    The "redo" is denied if the new guess is essentially the same
    as a previous guess.
    
    Parameters
    ----------
    x : int
      X of current position 
    y : int
      Y of current position 
    guessx : int
      X of guess position 
    guessy : int
      Y of guess position 
    guesspar : list/array
      Guess parameters 
     
    Parameters
    ----------
    flag : boolean,
      Is this redo okay.
        True - Redo okay. This guess has not been done before 
        False - Redo NOT okay. This guess has been done before 
     
    When there are any problems in this program it returns: 
    flag = -1 
     
    Created by David Nidever April 2005 
    Translated to python by D. Nidever, March 2022
    """

    global BTRACK, GSTRUC
    
    flag = True  # do it unless proven wrong 
     
    # FROM **ANY** PREVIOUS POSITION 
    prev, = np.where((BTRACK['x']==x) & (BTRACK['y']==y))
    nprev = len(prev)

    if guesspar is None:
        return False
    tguesspar = guesspar
    nguesspar = len(tguesspar) 
    ngg = nguesspar//3 
     
    # FROM **ANY** PREVIOUS POSITION 
    # we have used a guess from this position before 
    #  but have the parameters changed sufficiently 
    if (nprev > 0): 
         
        # Looping through the previous ones 
        for i in range(nprev):
            guesspar2 = BTRACK['data'][prev[i]]['guesspar']
         
            # Some gaussians found 
            if (guesspar2 is not None):
                tpar = guesspar2
                ntpar = len(tpar) 
                ntg = ntpar//3      # number of gaussians in this guess 
             
                # Same number of gaussians 
                if (ntpar == nguesspar): 
                    # Sorting, largest first 
                    tpar2 = utils.gsort(tpar) 
                    tguesspar2 = utils.gsort(tguesspar) 
                    
                    # Gixing possible zeros that could ruin the ratio 
                    dum = np.copy(tpar2)
                    bd, = np.where(dum == 0.) 
                    if len(bd) > 0:
                        dum[bd] = 1e-5 
                    diff = np.abs(tpar2 - tguesspar2) 
                    ratio = diff/np.abs(dum) 
                 
                    # These differences are too small, NO redo 
                    if (np.max(ratio) < 0.01): 
                        return False
 
    return flag 


def gbetter(res1,res2):
    """
    This function tries to figure out if one gaussian 
    analysis is better than another. 
     
    Parameters
    ----------
    res1 : dict
      Results from position 1 with par, rms, noise.
    res2 : dict
      Results from position 2 with par, rms, noise.
     
    Returns
    -------
    better : int
       The function value is either: 
         1 - Second position better than the first 
         0 - First position better than the second 
        -1 - Any problems 
     
    If the first one is better then it returns 0 
    and if the second one is better it returns 1 
     
    When this program has any problems is return: 
      better = -1 
     
    Created by David Nidever April 2005 
    Translated to python by D. Nidever, March 2022
    """
 
    better = -1   # default unless proven wrong 

    rms1,noise1,par1 = res1['rms'],res1['noise'],res1['par']
    rms2,noise2,par2 = res2['rms'],res2['noise'],res2['par']    
    
    # In case either one is -1 (bad)
    if par1 is not None and par2 is not None:
        if (rms1 == -1) and (rms2 != -1): 
            better = 1
        if (rms1 != -1) and (rms2 == -1): 
            better = 0 
        if (rms1 == -1) and (rms2 == -1): 
            better = -1 
        if (rms1 == -1) or (rms2 == -1): 
            return better
        if (len(par1) < 3) and (len(par2) >= 3): 
            better = 1 
        if (len(par2) < 3) and (len(par1) >= 3): 
            better = 0 
        if (len(par1) < 3) or (len(par2) < 3): 
            return better

    if par2 is None:
        return better
     
    drms1 = rms1-noise1 
    drms2 = rms2-noise2 
    n1 = len(par1)/3 
    n2 = len(par2)/3 
     
    # Clear cut, rms better, n equal or less 
    if (drms1 < drms2) and (n1 <= n2): 
        better = 0 
    if (drms1 > drms2) and (n1 >= n2): 
        better = 1 
     
    # rms same, n different 
    if (drms1 == drms2) and (n1 <= n2): 
        better = 0 
    if (drms1 == drms2) and (n1 > n2): 
        better = 1 
     
    # mixed bag, lower rms but higher n 
    if (drms1 < drms2) and (n1 > n2): 
        ddrms = drms2-drms1 
        rdrms = ddrms/drms2   # ratio compared to worse one 
        dn = n1-n2 
         
        better = 1    # default 
        if (dn == 1) and (rdrms > 0.2) : 
            better = 0 
        if (dn == 2) and (rdrms > 0.5) : 
            better = 0 
        if (dn == 3) and (rdrms > 1.0) : 
            better = 0 
        if (dn >= 4) and (rdrms > 2.0) : 
            better = 0 
     
    if (drms2 < drms1) and (n2 > n1): 
        ddrms = drms1-drms2 
        rdrms = ddrms/drms1    # ratio compared to worse one 
        dn = n2-n1 
         
        better = 0   # default 
        if (dn == 1) and (rdrms > 0.2) : 
            better = 1 
        if (dn == 2) and (rdrms > 0.5) : 
            better = 1 
        if (dn == 3) and (rdrms > 1.0) : 
            better = 1 
        if (dn >= 4) and (rdrms > 2.0) : 
            better = 1 
     
    return better 


def gfind(x,y,xr=None,yr=None):
    """
    This function helps find a y and x in 
    the gaussian components structure. 
     
    Parameters
    ----------
    x : int
      X to search for.
    y : int
      Y to search for.
    xr : list/array, optional
      Two element array of x limits, xr=[xmin,xmax].
    yr : list/array, optional
      Two element array of y limits, yr=[ymin,ymax].
     
    Returns
    -------
    flag : 
      The function value is either 0 or 1: 
         1 - the position exists in the structure 
         0 - the position does NOT exist in the structure 
        -1 - any problems 
    results : dict
      Dictionary of resuts with pind, rms, noise, par.
        pind    Index of position in GSTRUC. 
        rms     RMS of gaussian fit at the desired position 
        noise   Noise level at desired position 
        par     Parameters of gaussians in GSTRUC with the desired position 
     
    When there are any problems in this program it returns: 
     flag = None
     rms = None
     noise = None 
     par = None
     pind = None
     
    Created by David Nidever April 2005 
    Translated to python by D. Nidever, March 2022
    """

    global BTRACK, GSTRUC
    
    # Assume bad until proven otherwise 
    flag,rms,noise,par,pind = None,None,None,None,None
    results = {'pind':pind,'rms':rms,'noise':noise,'par':par}  # initial bad values

    if x is None or y is None:
        return 0,results
    
    # Setting the ranges
    if xr is not None:
        x0 = xr[0] 
        x1 = xr[1] 
    else: 
        x0 = 0. 
        x1 = 359.5
    if yr is not None:
        y0 = yr[0] 
        y1 = yr[1] 
    else: 
        y0 = -90. 
        y1 = 90. 
     
    if (x < x0) or (x > x1) or (y < y0) or (y > y1): 
        flag = -1 
        return flag,results
     
    # No GSTRUC yet, first position
    try:
        dum = len(GSTRUC)
    except:
        return 0,results
    
    # Looking for the position 
    t0 = time.time() 
    # XSTART/YSTART has a value for each position, faster searching 
    #  use NGAUSS and INDSTART to get the indices into DATA
    pind, = np.where((GSTRUC['x']==x) & (GSTRUC['y']==y))
    #print('find ',time.time()-t0)
    
    # Found something, getting the values 
    if len(pind) > 0:
        tstr = GSTRUC['data'][pind[0]]
        rms = tstr['rms']
        noise = tstr['noise']
        par = tstr['par']
        flag = 1 
         
    # Nothing found 
    else:
        pind,rms,noise,par = None,None,None,None
        flag = 0 

    results = {'pind':pind,'rms':rms,'noise':noise,'par':par}
    return flag,results


def gguess(x,y,xsgn=1,ysgn=1,xr=None,yr=None):
    """
    This program finds the best guess for the new profile 
    The one from the backwards positions. 
    If it can't find one then it returns 999999.'s 
     
    Parameters
    ----------
    x      X of current position 
    y      Y of current position 
    xsgn   Sign of x increment 
    ysgn   Sign of y increment 
    xr     Two element array of x limits 
    yr     Two element array of y limits 
     
    Returns
    -------
    guesspar  The first guess gaussian parameters 
    guessx  The x of the position where the 
                   guess parameters came from 
    guessy  The y of the position where the 
                   guess parameters came from 
     
    When this program has problems it returns: 
      guesspar = 999999. 
      guessx = 999999. 
      guessy = 999999. 
     
    Created by David Nidever April 2005 
    Translated to python by D. Nidever, March 2022
    """
     
    # Bad until proven okay 
    guesspar = None
    guessx = None 
    guessy = None
     
    # Making sure it's the right structure 
    #tags = tag_names(*(!gstruc.data)) 
    #if (len(tags) != 6) : 
    #    return guesspar,guessx,guessy 
    #comp = (tags == ['X','Y','RMS','NOISE','PAR','SIGPAR']) 
    #if ((where(comp != 1))(0) != -1) :
    #    return guesspar,guessx,guessy         
     
    # Saving the originals 
    orig_x = x 
    orig_y = y 

    if xr is None:
        xr = [0.,359.5]
    if yr is None:
        yr = [-90.,90.] 
     
    # Is the x range continuous?? 
    if (xr[0] == 0.) and (xr[1] == 359.5): 
        cont = 1 
    else: 
        cont = 0 
     
    # getting the p3 and p4 positions 
    # P3 back in x (l-0.5), same y 
    # P4 back in y (b-0.5), same x 
    x3,y3 = gincrement(x,y,xr,yr,xsgn=-xsgn,ysgn=-ysgn)
    x4,y4 = gincrement(x,y,xr,yr,xsgn=xsgn,ysgn=-ysgn,p2=True)
     
    # CHECKING OUT THE EDGES 
    # AT THE LEFT EDGE, and continuous, Moving RIGHT, use neighbor on other side 
    # Use it for the guess, but never move to it directly 
    if (x == xr[0]) and (xsgn == 1) and (cont == 1): 
        y3 = y 
        x3 = xr[1] 
     
    # AT THE RIGHT EDGE, and continuous, Moving LEFT 
    if (x == xr[1]) and (xsgn == -1) and (cont == 1): 
        y3 = y 
        x3 = xr[0] 
     
    # At the edge, NOT continuous, Moving RIGHT, NO P3 NEIGHBOR 
    if (x == xr[0]) and (xsgn == 1) and (cont == 0): 
        x3 = None 
        y3 = None 
     
    # At the edge, NOT continuous, Moving LEFT, NO P3 NEIGHBOR 
    if (x == xr[1]) and (xsgn == -1) and (cont == 0): 
        x3 = None
        y3 = None
     
    # Have they been visited before? 
    p3,res3 = gfind(x3,y3)
    p4,res4 = gfind(x4,y4)
     
    # Comparing the solutions 
    b34 = gbetter(res3,res4)
     
    # selecting the best guess 
    if (b34 == 0):# using P3 
        guesspar = res3['par']
        guessx = x3
        guessy = y3 
    if (b34 == 1):# using P4 
        guesspar = res4['par']
        guessx = x4 
        guessy = y4 
    if (b34 == -1): 
        guesspar = None
        guessx = None 
        guessy = None
     
    # Putting the originals back 
    x = orig_x 
    y = orig_y 

    return guesspar,guessx,guessy  


def nextmove(x,y,xr,yr,count,xsgn=1,ysgn=1,redo=False,redo_fail=False,back=False,
             noback=False,backret=True,wander=True):
    """
    Figure out the next move, the next position to decompose.

    Parameters
    ----------
    x : int
      X of current position 
    y : int
      Y of current position 
    xr : list/array
      Two element array of x limits.
    yr  : list/array
      Two element array of y limits.
    count : int
      The current iteration count.
    xsgn : int, optional
      Sign of x increment.  Default is 1.
    ysgn : int, optional
      Sign of y increment, Default is 1.
    redo : boolean, optional
      Was the current position a "redo"?  Default is False.
    redo_fail : boolean, optional
      If the current position was a "redo", was the redo successful
        (better than the previous attempt)?  Default is False.
    back : boolean, optional
      Was the current position a "backwards" step?  Default is False.
    noback : boolean, optional
      The program is not allowed to go backwards.  Default is False.
    backret : boolean, optional
      Any backwards motion must return to the position it 
         came from.  Default is True.
    wander : boolean, optional
      Allow backwards motion. Haud's algorithm.  Default is False.


    Returns
    -------
    newx : int
      X of new position 
    newy : int
      Y of new position 
    guessx : int
      X position for guess parameters.
    guessy : int
      X position for guess parameters.
    guesspar : list/array
      Guess parameter array.
    redo : boolean
      New position is a redo.
    skip : boolean
      Skip the next position.

    Example
    -------

    newx,newy,guessx,guessy,guesspar,redo,skip = nextmove(x,y,xr,yr,count,redo=redo,redo_fail=redo_fail,back=back)

    """

    global BTRACK, GSTRUC


    # If back, redo and BACKRET=1 then return to pre-redo position
    #=============================================================
    # This is done separately from the normal algorithm 
    if backret and back and redo: 
        back = False
        nbtrack = len(BTRACK)
        newx = BTRACK[-1]['data']['lastx']
        newy = BTRACK[-1]['data']['lasty']
        lastx = BTRACK[-1]['data']['x']
        lasty = BTRACK[-1]['data']['y']        

        # p0 is the redo position, p5 is the pre-redo position 
        p0,res0 = gfind(lastx,lasty,xr=xr,yr=yr)
        p5,res5 = gfind(newx,newy,xr=xr,yr=yr)

        b = gbetter(res0,res5)
        redo = gredo(newx,newy,lastx,lasty,par0) 

        # back position better, redo pre-redo position 
        if (b==0) and redo: 
            # Getting the guess 
            guesspar,guessx,guessy = gguess(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
            redo = True
            skip = False
        # back position worse, or can't redo pre-position, skip 
        else:
            redo = False
            skip = True

        return newx,newy,guessx,guessy,guesspar,redo,skip


    # Redo Failed! Return to original position
    #   If we went back and backret=1 then return to pre-redo positi
    #   if we went forward then don't do anything, should continue forward 
    if redo and redo_fail and back: 
        # Go back to pre-redo position and skip
        newx = BTRACK[-1]['data']['lastx']
        newy = BTRACK[-1]['data']['lasty']            
        return newx,newy,None,None,None,False,True
    

    # Some default values
    skip = False
    guesspar,guessx,guessy = None,None,None

    
    # Positions
    #
    # ^ Y    P2     
    # |   P3 P0 P1
    # |      P4
    #   --------> X
    
    # Get the positions, THIS IS THE PROPER WAY TO DO IT!!!!! 
    x1,y1 = gincrement(x,y,xr,yr,xsgn=xsgn,ysgn=ysgn)
    x2,y2 = gincrement(x,y,xr,yr,xsgn=xsgn,ysgn=ysgn,p2=True) 
    x3,y3 = gincrement(x,y,xr,yr,xsgn=-xsgn,ysgn=-ysgn)
    x4,y4 = gincrement(x,y,xr,yr,xsgn=xsgn,ysgn=-ysgn,p2=True)

    # Have they been visited before? 
    p0,res0 = gfind(x,y,xr=xr,yr=yr)
    par0 = res0['par']
    p1,res1 = gfind(x1,y1,xr=xr,yr=yr)
    p2,res2 = gfind(x2,y2,xr=xr,yr=yr)
    p3,res3 = gfind(x3,y3,xr=xr,yr=yr)
    p4,res4 = gfind(x4,y4,xr=xr,yr=yr)

    # PRINTING OUT SOME RELEVANT INFORMATION HERE 
    # comparing the solutions at neighboring positions
    b1 = gbetter(res0,res1)
    b2 = gbetter(res0,res2)
    b3 = gbetter(res0,res3)
    b4 = gbetter(res0,res4)

    # Do we need to redo? 
    if (p1==1) and (b1==0): 
        red1 = True
    else: 
        red1 = False
    if (p2==1) and (b2==0): 
        red2 = True
    else: 
        red2 = False
    if (p3==1) and (b3==0): 
        red3 = True
    else: 
        red3 = False
    if (p4==1) and (b4==0): 
        red4 = True
    else: 
        red4 = False

    xx = [x1,x2,x3,x4]
    yy = [y1,y2,y3,y4]
    pp = [p1,p2,p3,p4]
    bb = [b1,b2,b3,b4]
    rr = [red1,red2,red3,red4]
    
    # Printing out the info 
    print('Count = %d' % count)
    print('Last/Current Position = (%d,%d)' %(x,y))
    print('Neighbors (position)  visited  better  redo')
    for i in range(4):
        if xx[i] is not None:
            strx = '%5d' % xx[i]
        else:
            strx = '-----'
        if yy[i] is not None:
            stry = '%5d' % yy[i]
        else:
            stry = '-----'            
        print('P%1d (%5s,%5s)  %7d  %7d  %7s' % (i+1,strx,stry,pp[i],bb[i],str(rr[i])))    
    print('')

    # If P3 or P4 worse than P0 then move back to worst decomp 
    # If P3 and P4 better than P0 then move forward,
    #   -if both have been visited before then do the worst decomp 
    #   -if neither has been visited before then move to P1. 

    

    # Starting Normal Algorithm
    #  (not redo+back+backred)
    #==========================



    #==============================                    
    #---- CHECKING BACKWARDS ----
    #==============================                
    if ((p3==1) or (p4==1)) and (noback==False): 

        # Only P3 visited before
        #=======================
        if (p3==1) and (p4==0): 
            b3 = gbetter(res0,res3)
            # Checking to see if this has been done before 
            #   getting P3 position 
            redo = gredo(x3,y3,x,y,par0)

            # P3 worse than P0, moving back
            #------------------------------
            if (b3==0) and redo: 
                newx,newy = x3,y3
                back = True   # moving backwards 
                guesspar = par0
                guessx = x 
                guessy = y 
            else: 
                back = False
                redo = False

        # Only P4 visited before
        #=======================
        if (p4==1) and (p3==0): 
            b4 = gbetter(res0,res4)
            # Checking to see if this has been done before 
            #   getting P4 position 
            redo = gredo(x4,y4,x,y,par0)

            # P4 worse than P0, moving back
            #------------------------------
            if (b4==0) and redo: 
                newx,newy = x4,y4
                back = True   # moving backwards
                guessx,guessy,guesspar = x,y,par0
            else: 
                back = False
                redo = False

        # Both visited before
        #====================
        if (p3==1) and (p4==1): 
            b3 = gbetter(res0,res3)
            b4 = gbetter(res0,res4)
            redo = True    # redo unless proven otherwise 
            # Checking to see if this has been done before 
            #   getting P3 position 
            redo3 = gredo(x3,y3,x,y,par0) 
            # Checking to see if this has been done before 
            #   getting P4 position 
            redo4 = gredo(x4,y4,x,y,par0) 

            # P3 worse than P0, but P4 better than P0 (b3=0 and b4=1)
            #----------------------------------------
            if (b3==0) and (b4==1): 
                # We can redo it, moving back to P3 
                if redo3: 
                    newx,newy = x3,y3
                # Can't redo, move forward 
                else:
                    redo = False
                    back = False

            # P4 worse than P0, but P3 better than P0 (b3=1 and b4=0)
            #----------------------------------------
            if (b3==1) and (b4==1): 
                # We can redo it, moving back to P4 
                if redo4: 
                    newx,newy = x4,y4
                # Can't redo, move forward 
                else: 
                    redo = False
                    back = False

            # Both P3 and P4 are worse than P0
            #---------------------------------
            if (b3==0) and (b4==0): 
                # Can redo either one 
                if redo3 and redo4: 
                    b34 = gbetter(res3,res4)
                    # Moving back to P3 (P3 worse than P4) 
                    if (b34==1):   # to P3 
                        newx,newy = x3,y3
                    # Moving back to P4 (P4 worse than P3)
                    if (b34==0):   # to P4 
                        newx,newy = x4,y4
                # Can't redo P4, go to P3 
                if redo3 and (redo4 == False): 
                    newx,newy = x3,y3   # to P3
                # Can't redo P3, go to P4 
                if (redo3 == False) and redo4: 
                    newx,newy = x4,y4   # to P4 
                # Can't do either, move forward 
                if (redo3 == False) and (redo4 == False): 
                    redo = False 
                    back = False

            # Both are better than P0, move forward
            #--------------------------------------
            if (b3==1) and (b4==1): 
                back = False
                redo = False

            # One is worse than P0
            #---------------------
            if redo: 
                back = True  # moving backwards 
                guessx,guessy,guesspar = x,y,par0


    #==============================
    # ---- CHECKING FORWARD ----
    #==============================
    if ((p3==0) and (p4==0)) or (back == False) or noback: 

        # This is the very end 
        if (x1 is None): 
            flag = 1
            print('We have reached the very end')
            #goto, BOMB
            import pdb; pdb.set_trace()

        back = False  # moving forward 

        # Only P1 has been visited before
        #================================
        if (p1==1) and (p2==0): 
            b1 = gbetter(res0,res1)
            redo = True
            # Checking to see if this has been done before 
            #   getting P1 position 
            redo1 = gredo(x1,y1,x,y,par0) 

            # Moving to P1 (P1 worse than P0) 
            if (b1==0) and redo1: 
                newx,newy = x1,y1
                # getting the guess 
                guesspar,guessx,guessy = gguess(x,y,xr,yr,xsgn=xsgn,ysgn=ysgn)
            # Can't redo P1, or P1 better than P0, move another step ahead 
            else: 
                newx,newy = x1,y1
                redo = False
                skip = True  # don't fit this one 

        # Only P2 has been visited before, THIS SHOULD NEVER HAPPEN
        #================================                    
        if (p2==1) and (p1==0): 
            print('This should never happen!!')
            import pdb; pdb.set_trace() 

        # Both have been visited before
        #==============================                  
        if (p1==1) and (p2==1): 
            b1 = gbetter(res0,res1)
            b2 = gbetter(res0,res2)
            redo = True   # redo unless proven otherwise 
            # Checking to see if this has been done before 
            #   getting P1 position 
            redo1 = gredo(x1,y1,x,y,par0) 
            # Checking to see if this has been done before 
            #   getting P2 position 
            redo2 = gredo(x2,y2,x,y,par0) 
            if (redo1 == False) and (redo2 == False):  # no redo 
                redo = False

            # P1 worse than P0, and P2 better than P0 (b1=0 and b2=1)
            #----------------------------------------
            if (b1==0) and (b2==1): 
                # Can redo, moving to P1 
                if redo1: 
                    newx,newy = x1,y1
                # Can't redo, increment and skip 
                else: 
                    newx,newy = x1,y1  # to P1 
                    redo = False
                    skip = True

            # P2 worse than P0, and P1 better than P0 (b1=1 and b2=0)
            #----------------------------------------
            if (b1==1) and (b2==0): 
                # Can redo, moving to P2 
                if redo2: 
                    newx,newy = x2,y2
                # Can't redo, increment to P1 and skip 
                else: 
                    newx,newy = x1,y1  # to P1 
                    redo = False
                    skip = True


            # Both worse than P0
            #-------------------
            if (b1==0) and (b2==0):  # both bad, find worst 
                # Can redo either one 
                if redo1 and redo2: 
                    b12 = gbetter(res1,res2)
                    # Moving to P1 (P1 worse than P2) 
                    if (b12==1):  # to P1 
                        newx,newy = x1,y1
                    # Moving to P2 (P2 worse than P1) 
                    if (b12==0):  # to P2
                        newx,newy = x2,y2

                # Can't redo P2, go to P1 
                if redo1 and (redo2 == False): 
                    newx,newy = x1,y1  # to P1 
                # Can't redo P1, go to P2 
                if (redo1 == False) and redo2: 
                    newx,newy = x2,y2   # to P2 
                # Can't do either, increment to P1 and skip 
                if (redo1 == False) and (redo2 == False): 
                    newx,newy = x1,y1  # to P1 
                    redo = False
                    skip = True 

            # Both better, increment to P1 and skip
            #--------------------------------------
            if (b1==1) and (b2==1): 
                newx,newy = x1,y1    # to P1 
                redo = False
                skip = True 

            # Getting the guess 
            if redo:  # redo 
                # Getting the new guess from backward positions 
                gguesspar,guessx,guessy = gguess(newx,newy,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)


        # Neither has been visited before, increment
        #===========================================
        if (p1==0) and (p2==0): 
            # Increment to P1
            newx,newy = x1,y1
            # Getting the guess 
            guesspar,guessx,guessy = gguess(newx,newy,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)


    return newx,newy,guessx,guessy,guesspar,redo,skip



def driver(datacube,xstart=0,ystart=0,outfile=None,silent=False,noplot=False,
           plotxr=None,xsgn=1,ysgn=1,xr=None,yr=None,trackplot=False,
           noback=False,backret=True,gstruc=None,btrack=None,savestep=None,
           gassnum=None,subcube=None,wander=False,clobber=False):
    """
    This program runs the gaussian fitting program 
    on a large part of the HI all sky survey 
     
    This program can be run in three different modes. 
    backret: The program is allowed to go backwards, but MUST return to 
             the position it came from.  This seems to be the best mode 
             to run the program in for most purposes. 
             This is now the DEFAULT mode. 
    noback: The program is not allowed to go backwards.  It will always 
             move in a forward direction until it is finished.  Therefore 
             there is essentially no re-decomposition of positions. 
             To use this mode set noback=True on the command line. 
    wanter: The program is allowed to go backwards and forwards (Haud's 
             original algorithm).  If it goes back in y and only 
             the forward x position has been visited before then 
             it will need to go through all x positions before 
             returning to where it came from.  So either use strips 
             narrow in x or run this only once an initial solution 
             for every position has been found. 
             To use this mode set wander=True on the command line. 
     
    If you input gstruc and btrack from a previous run, then the program 
    will start where it left off.  This is useful when you ran it with 
    backret and then want to use the wander mode.  But currently it 
    start where it left off, so if you want it to start at the very 
    beginning again you need to do something (such as adding btrack(0) 
    to the end of btrack: btrack = [btrack,btrack(0)]. 
    This can also be used for debugging purposes. 
     
    Parameters
    ----------
    datacube : array or st
      Cube object or filename.
    xstart : int, optional
       The x to start with.  Default is 0.
    ystart : int, optional
       The y to start with.  Default is 0.
    xr : list/array
       Two-element X range.
    yr : list/array
       Two-element Y range.
    xsgn : int, optional
       Direction of x increments (-1 or 1).  Default is 1.
    ysgn : int, optional
       Direction of y increments (-1 or 1). Default is 1.
    outfile : str, optional
       File to save the structures to.
    plotxr : list
       Plotting xrange.
    trackplot : boolean, optional
      Track the progression visually 
    noplot : boolean, optional
      Don't plot anything.
    silent : boolean, optional
      Don't print anything 
    noback : boolean, optional
      The program is not allowed to go backwards.  Default is False.
    backret : boolean, optional
       Any backwards motion must return to the position it 
         came from.  Default is True.
    wander : boolean, optional
       Allow backwards motion. Haud's algorithm.  Default is False.
    gstruc : list, optional
       List of gaussians to start with.
    btrack : list, optional
       Tracking structure to start with.
     
    Returns
    -------
    gstruc : list
      List of all the gaussians found.
    btrack : list
      List that keeps track of every move. 

    Example
    -------

    gstruc,btrack = driver(cube)

    Created by David Nidever April 2005 
    Translated to python by D. Nidever, March 2022
    """

    global BTRACK, GSTRUC
    
    flag = 0 
    count = 0
    tstart = time.time() 


    # Load the cube
    if type(datacube) is str:
        datacubefile = datacube
        print('Loading '+datacubefile)
        datacube = Cube.read(datacubefile)
    
    # Setting parameters
    if xr is None:
        xr = [0,datacube.nx-1]
    if yr is None:
        yr = [0,datacube.ny-1]  
    if xstart is None:
        if xsgn == 1: 
            xstart = xr[0] 
        else: 
            xstart = xr[1]
    if ystart is None:
        if ysgn == 1 : 
            ystart = yr[0] 
        else: 
            ystart = yr[1]
    if wander:
        backret = False
    if noback:
        backret = False

    # No mode selected, using default mode (backret) 
    if (backret == False) and (noback == False) and (wander == False): 
        print('' )
        print('!!! WARNING !!!  NO MODE SELECTED  ->  USING DEFAULT (BACKRET) MODE' )
        print('')
        sleep(3) 
        backret = True


      

        
    # Restore file 
    #restore_file = repstr(outfile,'.fits','_restore.sav') 
 
    # Checking the file
    if outfile is None:
        logtime = datetime.now().strftime("%Y%m%d%H%M%S") 
        outfile = 'gaussdecomp_'+logtime+'.fits' 
    #dum = findfile(outfile) 
    #if dum != '': 
    #    print('THE FILE ',outfile,' EXISTS ALREADY !!!' )
    #    print('DO YOU WANT TO CONTINUE?')
    #    quest='' 
    #    read,quest 
    #    if quest != 'y' and quest != 'yes' and quest != 'YES' and     quest != 'yes' and quest != 'Yes' : 
    #    return
 
    # Printing out the inputs 
    print(' RUNNING GAUSSIAN ANALYSIS WITH THE FOLLOWING PARAMETERS')
    print('-----------------------------------------------------------')
    print(' STARTING POSITION = (%d,%d)' % (xstart,ystart))
    print(' X RANGE = [%d,%d]' % (xr[0],xr[1]))
    print(' Y RANGE = [%d,%d]' % (yr[0],yr[1]))
    print(' X DIRECTION = '+str(xsgn))
    print(' Y DIRECTION = '+str(ysgn))
    print(' OUTFILE = '+outfile)
    print('-----------------------------------------------------------')
    if (backret == 1) : 
        print(' USING (BACKRET) MODE')
    if (noback == 1) : 
        print(' USING (NOBACK) MODE')
    if (wander == 1) : 
        print(' USING (WANDER) MODE')
    print('-----------------------------------------------------------')
    print('')
 
    # Initializing some parameters 
    redo_fail = False 
    redo = False
    back = False
    lastx = None
    lasty = None
    p0 = False
    p1 = False
    p2 = False
    p3 = False
    p4 = False
     
    # Where are we starting 
    x = xstart 
    y = ystart 

    track_dict = {'count':None,'x':None,'y':None,'rms':None,'noise':None,'par':None,
                  'guesspar':None,'guessx':None,'guessy':None,'back':None,'redo':None,
                  'redo_fail':None,'skip':None,'lastx':None,'lasty':None}
    gstruc_dict = {'x':None,'y':None,'rms':None,'noise':None,'par':None,
                   'sigpar':None,'long':None,'lat':None}
     
    # STARTING THE LARGE LOOP 
    while (flag == 0): 
         
        t00 = time.time() 
         
        # P0 is the current position 
        # P1 forward in x (l+0.5), same y 
        # P2 forward in y (b+0.5), same x 
        # P3 back in x (l-0.5), same y 
        # P4 back in y (b-0.5), same x 
        # 
        # Move forward in x if possible 
         
        tstr,tstr1,tstr2,skip,guessx,guessy,guesspar = None,None,None,False,None,None,None        
         
        # STARTING WITH BTRACK, RESTORING THE LAST STATE 
        if (count == 0) and (gstruc is not None and btrack is not None):
            import pdb; pdb.set_trace()
            nbtrack = len(btrack) 
            count = btrack[nbtrack-1]['count']
            x = btrack[nbtrack-1]['x']
            y = btrack[nbtrack-1]['y']
            rms = btrack[nbtrack-1]['rms']
            noise = btrack[nbtrack-1]['noise']
            par = btrack[nbtrack-1]['par']
            guesspar = btrack[nbtrack-1]['guesspar']
            guessx = btrack[nbtrack-1]['guessx']
            guessy = btrack[nbtrack-1]['guessy']
            back = btrack[nbtrack-1]['back']
            redo = btrack[nbtrack-1]['redo']
            redo_fail = btrack[nbtrack-1]['redo_fail']
            skip = btrack[nbtrack-1]['skip']
            lastx = btrack[nbtrack-1]['lastx']
            lasty = btrack[nbtrack-1]['lasty']
            btrack_add(btrack)
            gstruc_add(gstruc)
             
            count += 1 
            lastx = x 
            lastlast = y 
         
         
        # FIGURE OUT THE NEXT MOVE 
        #------------------------- 
        if (count > 0):
            lastx,lasty = x,y
            x,y,guessx,guessy,guesspar,redo,skip = nextmove(x,y,xr,yr,count,xsgn,ysgn,backret=backret,noback=noback,
                                                            wander=wander,redo=redo,back=back,redo_fail=redo_fail)


        # MAYBE KEEP THE BTRACK AND GSTRUC *IN* THE CUBE OBJECT!!!!

        #import pdb; pdb.set_trace()
        
 
        # Starting the tracking structure, bad until proven good
        track = track_dict.copy()
        track['count'] = count 
        track['x'] = x 
        track['y'] = y 
        track['lastx'] = lastx 
        track['lasty'] = lasty 
        track['guesspar'] = guesspar 
        track['guessx'] = guessx 
        track['guessy'] = guessy 
        track['back'] = back 
        track['redo'] = redo 
        track['skip'] = skip 

        # Minimal structure, in case we skip
        tstr = {'x':x,'y':y,'rms':np.inf,'noise':None,'par':None,
                'sigpar':None,'lon':None,'lat':None}
        
        # Some bug checking 
        if x is None: 
            import pdb; pdb.set_trace() 
        if (x == lastx) and (y == lasty): 
            import pdb; pdb.set_trace() 
        #if count != 0: 
        #    if (red1+red2+red3+red4 == 0) and redo: 
        #        import pdb; pdb.set_trace() 
 
 
        if skip: 
            print('SKIP')
 
        # FITTING THE SPECTRUM, UNLESS WE'RE SKIPPING IT 
        #------------------------------------------------ 
        if skip == False: 
            t0 = time.time() 
            
            # Initial Printing 
            print('Fitting Gaussians to the HI spectrum at (%d,%d)' % (x,y))
            strout = ''
            if redo:
                strout = strout+'REDO '
            if back:
                strout = strout+'BACK'
            if back is False:
                strout = strout+'FORWARD' 
            print(strout) 
                          
            # Getting the HI spectrum
            spec = datacube(x,y)  # Get the new spectrum
            # No good spectrum 
            if spec is None or np.sum(spec.flux)==0:
                print('No spectrum to fit')
                rms = None
                noise = None
                skip = True
                count += 1
                gstruc_add(tstr)
                btrack_add(track)
                continue

            noise = spec.noise
            npts = spec.n
            
            # Zero-velocity region INCLUDED
            #====================================            
            if np.min(spec.vel) < 0:
                # GETTIING THE VELOCITY RANGE around the zero-velocity MW peak
                print('Zero-velocity region INCLUDED.  Fitting it separately')
                smspec = dln.savgol(spec.flux,21,2) 
                dum,vindcen = dln.closest(spec.vel,0)
            
                # Finding the vel. low point 
                flag = 0 
                i = vindcen
                lo = 0
                while (flag == 0): 
                    if smspec[i] <= noise: 
                        lo = i 
                    if smspec[i] <= noise: 
                        flag = 1
                    i -= 1 
                    if i < 0: 
                        flag = 1 
                lo = np.maximum(0,(lo-20))
 
                # Finding the vel. high point 
                flag = 0 
                i = vindcen
                hi = npts-1
                while (flag == 0): 
                    if smspec[i] <= noise : 
                        hi = i 
                    if smspec[i] <= noise : 
                        flag = 1 
                    i += 1 
                    if i > npts-1: 
                        flag = 1 
                hi = np.minimum((npts-1),(hi+20))
 
                vmin = spec.vel[lo] 
                vmax = spec.vel[hi] 
 
                # RUNNING GAUSSFITTER ON ZERO VELOCITY REGION, WITH GUESS 
                results = fitter.gaussfitter(spec,vmin=vmin,vmax=vmax,initpar=guesspar,silent=True,noplot=True)            
 
                # FIT WITH NO GUESS (if first time and previous fit above with guess) 
                tp0 = gfind(x,y,xr=xr,yr=yr) 
                if (tp0 == 0) and (guesspar is not None):
                    results2 = fitter.gaussfitter(spec,vmin=vmin,vmax=vmax,silent=True,noplot=True)
                    #tpar0,tsigpar0,trms,noise,v2,spec2,resid2 = fitter.gaussfitter(spec,vmin=vmin,vmax=vmax,silent=True,noplot=True)
                    #b = gbetter({'par':par0,'rms':rms,'noise':noise},{'par':tpar0,'rms':trms,'noise':noise})
                    b = gbetter(results,results2)
                    # The fit without the guess is better 
                    if (b == 1): 
                        par0 = tpar0 
                        sigpar0 = tsigpar0 
                        rms = trms 
 
                # ADDING THE BEST RESULTS TO THE STRUCTURE, TSTR1 
                if (par0(0) != -1): 
                    ngauss = len(par0)/3
                    tstr1 = replicate(gstruc_schema,ngauss) 
                    for i in range(ngauss): 
                        tstr1[i].par = par0[3*i:3*i+3]
                        tstr1[i].sigpar = sigpar0[3*i:3*i+3]
                    tstr1['x'] = x 
                    tstr1['y '] = y 
                    tstr1['lon'] = glon 
                    tstr1['lat'] = gy 
                    #tstr1['rms'] = rms 
                    tstr1['noise'] = noise 
 
                # REMOVING ZERO-VELOCITY parameters and spectrum 
                if par0[0] != -1: 
                    th = gfunc(v,par0) 
                    inspec = spec-th 
                    inv = v 
                    npts = len(v) 
                    if guesspar is not None:
                        inpar1 = guesspar 
                        inpar2 = guesspar 
                        inpar1 = gremove(inpar1,v[0:lo],spec[0:lo]) 
                        inpar2 = gremove(inpar2,v[hi:npts],spec[hi:npts])
                        guesspar2 = np.hstack((inpar1,inpar2))
                        gd, = np.where(guesspar2 != -1,ngd) 
                        if (len(gd) == 0):  # no guess 
                            guesspar2 = None
                            guesspar2 = guesspar2[gd]
                else: 
                    inspec = spec 
                    inv = v 
 
                
                # RUNNING GAUSSFITTER ON EVERYTHING WITHOUT THE ZERO-VELOCITY REGION, WITH GUESS 
                #par0,sigpar0,rms,noise,v3,spec3,resid3 = fitter.gaussfitter(inv,inspec,initpar=guesspar2,silent=True,noplot=True)
                results3 = fitter.gaussfitter(inv,inspec,initpar=guesspar2,noplot=True,silent=True)
            
 
                # FIT WITH NO GUESS (if first time and previous fit above with guess) 
                if (tp0 == 0) and (len(guesspar) > 1):
                    #tpar0,tsigpar0,trms,noise,v3,spec3,resid3 = fitter.gaussfitter(inv,inspec,silent=True,noplot=True)
                    results4 = fitter.gaussfitter(inv,inspec,silent=True,noplot=True)                    
                    #b = gbetter({'par':par0,'rms':rms,'noise':noise},{'par':tpar0,'rms':trms,'noise':noise})
                    b = gbetter(results3,results4)
                    # The fit without the guess is better 
                    if (b == 1):
                        results = results4.copy()
                    else:
                        results = results3.copy()


 
                # ADDING THE RESULTS TO THE STRUCTURE, TSTR2 
                if par0(0) != -1: 
                    ngauss = len(par0)/3 
                    tstr2 = replicate(gstruc_schema,ngauss) 
                    for i in range(ngauss): 
                        tstr2[i].par = par0[3*i:3*i+3]
                        tstr2[i].sigpar = sigpar0[3*i:3*i+3]
                    tstr2['x'] = x 
                    tstr2['y'] = y 
                    tstr2['lon'] = lon 
                    tstr2['lat'] = lat 
                    tstr2['noise'] = noise 
 
                # ADDING THE STRUCTURES TOGETHER, TSTR = [TSTR1,TSTR2]
                if tstr1 is not None and tstr2 is not None:
                    tstr = [tstr1,tstr2]
                if tstr1 is not None and tstr2 is None:
                    tstr = tstr1
                if tstr1 is None and tstr2 is not None:
                    tstr = tstr2
                if tstr1 is None and tstr2 is None:  # no gaussians
                    tstr = gstruc_dict.copy()
                    tstr['x'] = x 
                    tstr['y'] = y 
                    tstr['lon'] = lon 
                    tstr['lat'] = lat 
                    tstr['rms'] = rms 
                    tstr['noise'] = noise 

                        
            # Does NOT cover zero-velocity region
            #====================================
            else:
                print('Zero-velocity NOT covered')
                invspec = spec.copy()

                # RUNNING GAUSSFITTER ON EVERYTHING WITH GUESS 
                results = fitter.gaussfitter(invspec,initpar=guesspar,noplot=True)    #silent=True)            
            
                # FIT WITH NO GUESS (if first time and previous fit above with guess)
                tp0 = gfind(x,y,xr=xr,yr=yr)
                if (tp0 == 0) and (guesspar is not None):
                    results2 = fitter.gaussfitter(invspec,silent=True,noplot=True)                    
                    #b = gbetter({'par':par0,'rms':rms,'noise':noise},{'par':tpar0,'rms':trms,'noise':noise})
                    b = gbetter(results1,results2)
                    # The fit without the guess is better 
                    if (b == 1): 
                        results = results2.copy()                        
                         
                # Creating the structure with the results
                if results['par'] is not None:
                    ngauss = len(results['par'])//3
                    tstr = gstruc_dict.copy()
                    tstr['x'] = x
                    tstr['y'] = y
                    for n in ['par','sigpar','rms','noise']:
                        tstr[n] = results[n]
                    #tstr['lon'] = lon 
                    #tstr['lat'] = lat 
                else:
                    tstr = None
                    
                
            print('fitting ',time.time()-t0)
 
            # PLOTTING/PRINTING, IF THERE WAS A FIT 
            if tstr is not None:
                # Getting the rms of all the components of the whole spectrum 
                th = gfunc(v,tstr['par'])
                rms = np.std(spec-th) 
                tstr['rms'] = rms 
 
                # Printing and plotting
                if noplot == False:
                    utils.gplot(v,spec,tstr.par,xlim=plotxr)
                if silent == False:
                    utils.printgpar(tstr['par'],tstr['sigpar'],
                                    len(tstr['par'])//3,tstr['rms'],tstr['noise'])
                if trackplot:
                    utils.gtrackplot(x,y,lastx,lasty,redo, count,xr=xr,yr=yr,pstr=pstr,xstr=xstr,ystr=ystr)
            else:
                if silent == False:
                    print('No gaussians found at this position!')

 
                
            # ADDING SOLUTION TO GSTRUC
            if tstr is not None:
                if count == 0: 
                    gstruc_add(tstr)
                if count > 0: 
                    old,res1 = gfind(x,y,xr=xr,yr=yr)
 
                    # This is a re-decomposition 
                    if (old==1) and redo: 
                        # Checking the two decompositions 
                        #par2 = tstr['par']  # new one 
                        #rms2 = tstr['rms']
                        #b = gbetter({'par':par2,'rms':rms2,'noise':noise2},{'par':par1,'rms':rms1,'noise':noise1})
                        b = gbetter(tstr,res1)
                        # New one is better 
                        if (b == False): 
                            gstruc_replace(tstr)  # replacing the solution
                            t1 = time.time() 
                            print(time.time()-t1)
                            redo_fail = False
                        else: # re-decomposition failed 
                            redo_fail = True
                            print('REDO FAILED!')
 
                    # This is NOT a re-decomposition, add it 
                    if (old==0) or (redo == False): 
                        t1 = time.time() 
                        gstruc_add(tstr)
                        print('gstruc ',time.time()-t1)
                        redo_fail = False
 
        # SKIP FITTING PART
        else: 
            # Creating a dummy structure 
            tstr = None
            redo_fail = False
            redo = False
            back = False
 
            if trackplot:
                utils.gtrackplot(x,y,lastx,lasty,redo,count,xr=xr,yr=yr,pstr=pstr,xstr=xstr,ystr=ystr)
 
 
        # FINISHING UP THE TRACKING STRUCTURE
        if tstr is not None:
            npar = len(tstr['par'])
            track['par'] = tstr['par']            
        else:
            npar = 0
        track['rms'] = rms 
        track['noise'] = noise 
        track['redo_fail'] = redo_fail 
 
        # UPDATING THE TRACKING STRUCTURE
        btrack_add(track)
 
        count += 1 
 
        # Saving the last position 
        lastx = x 
        lasty = y 
 
        print('This iteration ',time.time()-t00)
 
        # SAVING THE STRUCTURES, periodically
        if savestep == False:
            savestep = 50 
            nsave = savestep 
            if (int(count)//int(nsave) == int(count)/float(nsave)): 
                print('SAVING DATA!')
                #MWRFITS,(*(!gstruc.data))[0:!gstruc.count-1],file,/create 
                #MWRFITS,(*(!btrack.data))[0:!btrack.count-1],file,/silent# append 
                #gstruc = !gstruc & btrack = !btrack 
                #SAVE,gstruc,btrack,file=restore_file 
                #undefine,gstruc,btrack 

                # Pickle the btrack structure
                         
    # FINAL SAVE 
    print(str(len(gstruc),2),' final Gaussians')
    print('Saving data to ',file)
    #MWRFITS,(*(!gstruc.data))[0:!gstruc.count-1],file,/create 
    #MWRFITS,(*(!btrack.data))[0:!btrack.count-1],file,/silent# append 
    #gstruc = !gstruc & btrack = !btrack 
    #SAVE,gstruc,btrack,file=restore_file 
    #undefine,gstruc,btrack 
 
    print('dt = ',str(time.time()-tstart,2),' sec.')


