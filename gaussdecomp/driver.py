#!/usr/bin/env python

import os
import time
import numpy as np
from datetime import datetime
from dlnpyutils import utils as dln
from . import utils,fitter
from .cube import Cube

# Tracking lists
BTRACK = []
GSTRUC = []

def gincrement(x,y,xr=None,yr=None,
               xsgn=1,ysgn=1,nstep=1,p2=False):
    """
    This program increments the position 
     
    Parameters
    ----------
    x      X of current position 
    y      Y of current position 
    xsgn   Sign of x increment 
    ysgn   Sign of y increment 
    xr     Two element array of x limits 
    yr     Two element array of y limits 
    nstep  Number of steps to increment 
    p2     Increment in y rather than in x
     
    Returns
    -------
    newx   X of new position 
    newy   Y of new position 
     
    When this program has problems it returns: 
    newx = 999999. 
    newy = 999999. 
     
    Created by David Nidever April 2005 
    """
 
    step = 1
     
    # Bad until proven okay 
    newx = None
    newy = None

    if xr is None:
        xr = [0.,2000.] 
    if yr is None: 
        yr = [0.,2000.] 
     
    # X and y limits 
    x0 = float(xr[0]) 
    x1 = float(xr[1]) 
    y0 = float(yr[0]) 
    y1 = float(yr[1]) 
     
    # Are we in the box? 
    if (x < x0) or (x > x1) or (y < y0) or (y > y1):
        return None,None

    # Wrong signs 
    if np.abs(xsgn) != 1 : 
        return None,None
    if np.abs(ysgn) != 1 :
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
     
    tx = float(x) 
    ty = float(y) 
     
    # Looping through all the steps 
    for i in range(nstep): 
     
        # p2, incrementing vertically 
        if p2:
            # UP, at the end 
            if (ysgn == 1) and (ty == y1):
                newx,newy = None,None
                continue
            # DOWN, at the end 
            if (ysgn == -1) and (ty == y0): 
                newx,newy = None,None
                continue 
         
            # Not at end, normal increment 
            newx = tx 
            newy = ty + ysgn * step 
         
        # Incrementing Sideways 
        else: 
            # RIGHT, xsgn = +1 
            if (xsgn == 1): 
                # UP, the very end 
                if (ysgn == 1) and (tx == x1) and (ty == y1): 
                    newx,newy = None,None
                    continue
                # DOWN, the very end 
                if (ysgn == -1) and (tx == x1) and (ty == y0): 
                    newx,newy = None,None
                    continue 
             
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
                    newx,newy = None,None
                    continue
                # DOWN, the very end 
                if (ysgn == -1) and (tx == x0) and (ty == y0): 
                    newx,newy = None,None
                    continue 
             
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
    This function checks wether we can redo this location again 
    returns True if this redo is okay 
    returns False if this redo is NOT okay 
    
    Parameters
    ----------
    x         X of current position 
    y         Y of current position 
    guessx    X of guess position 
    guessy    Y of guess position 
    guesspar  Guess parameters 
     
    Parameters
    ----------
    flag      Function value is either: 
                  True - Redo okay. This guess has not been done before 
                  False - Redo NOT okay. This guess has been done before 
     
    When there are any problems in this program it returns: 
    flag = -1 
     
    Created by David Nidever April 2005 
    """
     
    # Bad to start out with 
    flag = False
     
    # Making sure it's the right structure 
    #tags = tag_names(*(!btrack.data)) 
    if (len(tags) != 15) : 
        return -1 
    btags = ['COUNT','X','Y','RMS','NOISE','PAR','GUESSPAR','GUESSX','GUESSY',
             'BACK','REDO','REDO_FAIL','SKIP','LASTX','LASTY'] 
    comp = (tags == btags) 
    if ((where(comp != 1))(0) != -1) : 
        return -1 
     
    flag = True # do it unless proven wrong 
     
    # FROM **ANY** PREVIOUS POSITION 
    #prev, = np.where((*(!btrack.data)).x == x and (*(!btrack.data)).y == y)
    nprev = len(prev)
     
    gd1, = np.where(guesspar != 999999.) 
    if (len(gd1) == 0): 
        return False 
    tguesspar = guesspar[gd1]
    nguesspar = len(tguesspar) 
    ngg = nguesspar//3 
     
    # FROM **ANY** PREVIOUS POSITION 
    # we have used a guess from this position before 
    #  but have the parameters changed sufficiently 
    if (nprev > 0): 
         
        # Looping through the previous ones 
        for i in range(nprev): 
            #guesspar2 = (*(!btrack.data))[prev[i]].guesspar 
            gd2, = np.where(guesspar2 != 999999.,ngd2) 
         
            # Some gaussians found 
            if (len(gd2) > 0): 
                tpar = guesspar2[gd2]
                ntpar = len(tpar) 
                ntg = ntpar//3      # of gaussians in this guess 
             
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


def gbetter(par1,rms1,noise1,par2,rms2,noise2):
    """
    This function tries to figure out if one gaussian 
    analysis is better than another. 
     
    Parameters
    ----------
    par1    Gaussian Parameters of position 1 
    rms1    RMS of position 1 
    noise1  Noise level of position 1 
    par2    Gaussian Parameters of position 2 
    rms2    RMS of position 2 
    noise2  Noise level of position 2 
     
    Returns
    -------
    better  The function value is either: 
                1 - Second position better than the first 
                0 - First position better than the second 
               -1 - Any problems 
     
    If the first one is better then it returns 0 
    and if the second one is better it returns 1 
     
    When this program has any problems is return: 
      better = -1 
     
    Created by David Nidever April 2005 
    """
 
    better = -1   # default unless proven wrong 
     
    # In case either one is -1 (bad) 
    if len(par1) > 0 and len(par2) > 0: 
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
    x     X to search for 
    y     Y to search for 
    xr    Two element array of x limits, xr=[xmin,xmax] 
    yr    Two element array of y limits, yr=[ymin,ymax] 
     
    Returns
    -------
    flag    The function value is either 0 or 1: 
                 1 - the position exists in the structure 
                 0 - the position does NOT exist in the structure 
                -1 - any problems 
    pind     Index of position in GSTRUC. 
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
    """

    # Assume bad until proven otherwise 
    flag,rms,noise,par,pind = None,None,None,None,None
     
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
        rms = -1 
        noise = -1 
        par = -1 
        pind = -1 
        return flag,pind,rms,noise,par
     
    # No GSTRUC yet, first position
    try:
        dum = len(GSTRUC)
    except:
        return 0,None,None,None,None
    
    # Looking for the position 
    t0 = time.time() 
    # XSTART/YSTART has a value for each position, faster searching 
    #  use NGAUSS and INDSTART to get the indices into DATA 
    #pind, = np.where(*(!gstruc['xstart']) == x and *(!gstruc['ystart']) == y)
    npind = len(pind)
    print('find ',time.time()-t0)
     
    # Found something, getting the values 
    if npind > 0: 
        #ind = l64indgen((*(!gstruc['ngauss'])[pind[0]])+(*(!gstruc['indstart'])[pind[0]] 
        #rms = first_el((*(!gstruc['data'])[ind].rms) 
        #noise = first_el((*(!gstruc['data'])[ind].noise) 
        #par = ((*(!gstruc.data))[ind].par)(*) 
        flag = 1 
         
    # Nothing found 
    else: 
        rms,noise,par = None,None,None
        flag = 0 
     
    return flag,pind,rms,noise,par


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
    x3,y3 = gincrement(x,y,xsgn=-xsgn,ysgn=-ysgn,xr=xr,yr=yr)
    x4,y4 = gincrement(x,y,xsgn=xsgn,ysgn=-ysgn,xr=xr,yr=yr,p2=True)
     
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
    p3 = gfind(x3,y3,ind=ind3,rms=rms3,noise=noise3,par=par3) 
    p4 = gfind(x4,y4,ind=ind4,rms=rms4,noise=noise4,par=par4) 
     
    # Comparing the solutions 
    b34 = gbetter(par3,rms3,noise3,par4,rms4,noise4) 
     
    # selecting the best guess 
    if (b34 == 0):# using P3 
        guesspar = par3 
        guessx = x3 
        guessy = y3 
    if (b34 == 1):# using P4 
        guesspar = par4 
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


def nextmove(x,y,xsgn,ysgn,xr,yr):
    """
    Figure out the next move, the next position to decompose.

    """

    # Get the positions, THIS IS THE PROPER WAY TO DO IT!!!!! 
    x1,y1 = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
    x2,y2 = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr,p2=True) 
    x3,y3 = gincrement(x,y,xsgn=-xsgn,ysgn=-ysgn,xr=xr,yr=yr)
    x4,y4 = gincrement(x,y,xsgn=xsgn,ysgn=-ysgn,xr=xr,yr=yr,p2=True)

    # Have they been visited before? 
    p0 = gfind(x,y,rms=rms0,noise=noise0,par=par0,xr=xr,yr=yr) 
    p1 = gfind(x1,y1,rms=rms1,noise=noise1,par=par1,xr=xr,yr=yr) 
    p2 = gfind(x2,y2,rms=rms2,noise=noise2,par=par2,xr=xr,yr=yr) 
    p3 = gfind(x3,y3,rms=rms3,noise=noise3,par=par3,xr=xr,yr=yr) 
    p4 = gfind(x4,y4,rms=rms4,noise=noise4,par=par4,xr=xr,yr=yr) 

    # PRINTING OUT SOME RELEVANT INFORMATION HERE 
    # comparing them 
    strb1 = gbetter(par0,rms0,noise0,par1,rms1,noise1) 
    strb2 = gbetter(par0,rms0,noise0,par2,rms2,noise2) 
    strb3 = gbetter(par0,rms0,noise0,par3,rms3,noise3) 
    strb4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 

    # Do we need to redo? 
    if p1  and (strb1 == False): 
        red1 = True
    else: 
        red1 = False
    if p2 and (strb2 == False): 
        red2 = True
    else: 
        red2 = False
    if p3 and (strb3 == False): 
        red3 = True
    else: 
        red3 = 0 
    if p4 and (strb4 == False): 
        red4 = True
    else: 
        red4 = False

    strx1 = '%5.1f' % x1
    stry1 = '%5.1f' % y1
    strx2 = '%5.1f' % x2
    stry2 = '%5.1f' % y2
    strx3 = '%5.1f' % x3
    stry3 = '%5.1f' % y3
    strx4 = '%5.1f' % x4
    stry4 = '%5.1f' % y4
    if (x1 == 999999.): strx1 = '-----'
    if (y1 == 999999.): stry1 = '-----' 
    if (x2 == 999999.): strx2 = '-----'
    if (y2 == 999999.): stry2 = '-----' 
    if (x3 == 999999.): strx3 = '-----'
    if (y3 == 999999.): stry3 = '-----' 
    if (x4 == 999999.): strx4 = '-----'
    if (y4 == 999999.): stry4 = '-----' 

    # Printing out the info 
    print('Count = %d' % count)
    print('Last/Current Position = (%.1f,%.2f)' %(x,y))
    print('Neighbors (position)  visited  better  redo')
    print('P1  (',strx1,',',stry1,')  ',p1, strb1, red1) 
    print('P2  (',strx2,',',stry2,')  ',p2, strb2, red2)
    print('P3  (',strx3,',',stry3,')  ',p3, strb3, red3)
    print('P4  (',strx4,',',stry4,')  ',p4, strb4, red4)
    print('')

    # if P3 or P4 worse than P0 then move back, to worst decomp 
    # if P3 and P4 better than P0 then move forward, if both 
    #  have been visited before then the worst decomp 
    #  if neither has been visited before then move to P1. 


    # If back redo and BACKRET=1 then return to pre-redo position 
    # This is done separately from the normal algorithm 
    if backret and back and redo: 
        back = False
        nbtrack = len(BTRACK)
        newx = BTRACK[-1]['data']['lastx']
        newy = BTRACK[-1]['data']['lasty']
        lastx = BTRACK[-1]['data']['x']
        lasty = BTRACK[-1]['data']['y']        
        #newx = (*(!btrack.data))[nbtrack-1].lastx 
        #newy = (*(!btrack.data))[nbtrack-1].lasty 
        #lastx = (*(!btrack.data))[nbtrack-1].x 
        #lasty = (*(!btrack.data))[nbtrack-1].y 

        # p0 is the redo position, p5 is the pre-redo position 
        p0 = gfind(lastx,lasty,rms=rms0,noise=noise0,par=par0,xr=xr,yr=yr) 
        p5 = gfind(newx,newy,rms=rms5,noise=noise5,par=par5,xr=xr,yr=yr) 

        b = gbetter(par0,rms0,noise0,par5,rms5,noise5) 
        redo = gredo(newx,newy,lastx,lasty,par0) 
        x = newx 
        y = newy 

        # back position better, redo pre-redo position 
        if (b == False) and redo: 
            # Getting the guess 
            guesspar,guessx,guessy = gguess(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
            redo = True
            skip = False
        # back position worse, or can't redo pre-position, skip 
        else:
            redo = False
            skip = True

    # NOT (redo back and /backret) 
    # Starting Normal Algorithm 
    else: 
        # Redo Failed! 
        # if we went forward then don't do anything, should continue forward 
        # If we went back and backret=1 then return to pre-redo position 
        if redo and redo_fail and back: 
            # Go back to pre-redo position
            x = BTRACK[-1]['data']['lastx']
            y = BTRACK[-1]['data']['lasty']            
            #nbtrack = btrack.count 
            #x = btrack['data'][nbtrack-1].lastx 
            #y = btrack['data'][nbtrack-1].lasty 


        #==============================                    
        #---- CHECKING BACKWARDS ----
        #==============================                
        if ((p3 == 1) or (p4 == 1)) and (noback == 0): 

            # Only P3 visited before
            #=======================
            if p3 and (p4 == False): 
                b3 = gbetter(par0,rms0,noise0,par3,rms3,noise3) 
                # Checking to see if this has been done before 
                #   getting P3 position 
                tnewx,tnewy = gincrement(x,y,xsgn=-xsgn,ysgn=-ysgn,xr=xr,yr=yr)
                redo = gredo(tnewx,tnewy,x,y,par0)

                # P3 worse than P0, moving back
                #------------------------------
                if (b3 == False) and redo: 
                    newx,newy = gincrement(x,y,xsgn=-xsgn,ysgn=-ysgn,xr=xr,yr=yr)           
                    back = True   # moving backwards 
                    guesspar = par0 
                    guessx = x 
                    guessy = y 
                    x = newx 
                    y = newy 
                else: 
                    back = False
                    redo = False

            # Only P4 visited before
            #=======================
            if p4 and (p3 == False): 
                b4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 
                # Checking to see if this has been done before 
                #   getting P4 position 
                tnewx,tnewy = gincrement(x,y,xsgn=xsgn,ysgn=-ysgn,xr=xr,yr=yr,p2=True)
                redo = gredo(tnewx,tnewy,x,y,par0)

                # P4 worse than P0, moving back
                #------------------------------
                if (b4 == False) and redo: 
                    newx,newy = gincrement(x,y,xsgn=xsgn,ysgn=-ysgn,xr=xr,yr=yr,p2=True)
                    back = True  # moving backwards 
                    guesspar = par0 
                    guessx = x 
                    guessy = y 
                    x = newx 
                    y = newy 
                else: 
                    back = False
                    redo = False

            # Both visited before
            #====================
            if p3 and p4: 
                b3 = gbetter(par0,rms0,noise0,par3,rms3,noise3) 
                b4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 
                redo = 1 # redo unless proven otherwise 
                # Checking to see if this has been done before 
                #   getting P3 position 
                tnewx3,tnewy3 = gincrement(x,y,xsgn=-xsgn,ysgn=-ysgn,xr=xr,yr=yr)
                redo3 = gredo(tnewx3,tnewy3,x,y,par0) 
                # Checking to see if this has been done before 
                #   getting P4 position 
                tnewx4,tnewy4 = gincrement(x,y,xsgn=xsgn,ysgn=-ysgn,xr=xr,yr=yr,p2=True)
                redo4 = gredo(tnewx4,tnewy4,x,y,par0) 

                # P3 worse than P0, but P4 better than P0 (b3=0 and b4=1)
                #----------------------------------------
                if (b3 == False) and b4: 
                    # We can redo it, moving back to P3 
                    if (redo3 == 1): 
                        newx,newy = gincrement(x,y,xsgn=-xsgn,ysgn=-ysgn,xr=xr,yr=yr)
                    # Can't redo, move forward 
                    else: 
                        redo = False
                        back = False

                # P4 worse than P0, but P3 better than P0 (b3=1 and b4=0)
                #----------------------------------------
                if b3 and b4: 
                    # We can redo it, moving back to P4 
                    if redo4: 
                        newx,newy = gincrement(x,y,xsgn=xsgn,ysgn=-ysgn,xr=xr,yr=yr,p2=True)
                    # Can't redo, move forward 
                    else: 
                        redo = False
                        back = False

                # Both P3 and P4 are bad
                #-----------------------
                if (b3 == False) and (b4 == False): 
                    # Can redo either one 
                    if redo3 and redo4: 
                        b34 = gbetter(par3,rms3,noise3,par4,rms4,noise4)
                        # moving back to P3 (P3 worse than P4) 
                        if b34:   # to P3 
                            newx,newy = gincrement(x,y,xsgn=-xsgn,
                                                       ysgn=-ysgn,xr=xr,yr=yr)
                        # moving back to P4 (P4 worse than P3) 
                        if (b34 == False) :# to P4 
                            newx,newy = gincrement(x,y,xsgn=xsgn,
                                                       ysgn=-ysgn,xr=xr,yr=yr,p2=True)
                    # Can't redo P4, goto P3 
                    if redo3 and (redo4 == False): 
                        newx,newy = gincrement(x,y,xsgn=-xsgn,    # to P3 
                                                   ysgn=-ysgn,xr=xr,yr=yr) 
                    # Can't redo P3, goto P4 
                    if (redo3 == False) and redo4: 
                        newx,newy = gincrement(x,y,newx,newy,xsgn=xsgn,   # to P4 
                                                   ysgn=-ysgn,xr=xr,yr=yr,p2=True)
                    # Can't do either, move forward 
                    if (redo3 == 0) and (redo4 == 0): 
                        redo = False 
                        back = False

                # Both are better than P0, move forward
                #--------------------------------------
                if (b3 == 1) and (b4 == 1): 
                    back = False
                    redo = False

                # One is worse than P0
                #---------------------
                if redo: 
                    back = True  # moving backwards 
                    guesspar = par0 
                    guessx = x 
                    guessy = y 
                    x = newx 
                    y = newy 


        #==============================
        # ---- CHECKING FORWARD ----
        #==============================
        if ((p3 == False) and (p4 == False)) or (back == False) or noback: 

            # This is the very end 
            if (x1 == 999999.): 
                flag = 1 
                goto, BOMB 

            back = False  # moving forward 

            # Only P1 has been visited before
            #================================
            if p1 and (p2 == False): 
                b1 = gbetter(par0,rms0,noise0,par1,rms1,noise1) 
                redo = True
                # Checking to see if this has been done before 
                #   getting P1 position 
                tnewx1,tnewy1 = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
                redo1 = gredo(tnewx1,tnewy1,x,y,par0) 

                # Moving to P1 (P1 worse than P0) 
                if (b1 == False) and redo1: 
                    newx,newy = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
                    x = tnewx1 
                    y = tnewy1 
                    # getting the guess 
                    guesspar,guessx,guessy = gguess(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
                # Can't redo P1, or P1 better than P0, move another step ahead 
                else: 
                    x = tnewx1 
                    y = tnewy1 
                    redo = False
                    skip = True  # don't fit this one 

            # Only P2 has been visited before, THIS SHOULD NEVER HAPPEN
            #================================                    
            if p2 and (p1 == False): 
                print('This should never happen!!')
                import pdb; pdb.set_trace() 

            # Both have been visited before
            #==============================                  
            if p1 and p2: 
                b1 = gbetter(par0,rms0,noise0,par1,rms1,noise1) 
                b2 = gbetter(par0,rms0,noise0,par2,rms2,noise2) 
                redo = True   # redo unless proven otherwise 
                # Checking to see if this has been done before 
                #   getting P1 position 
                tnewx1,tnewy1 = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
                redo1 = gredo(tnewx1,tnewy1,x,y,par0) 
                # Checking to see if this has been done before 
                #   getting P2 position 
                tnewx2,tnewy2 = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr,p2=True)
                redo2 = gredo(tnewx2,tnewy2,x,y,par0) 
                if (redo1 == False) and (redo2 == False):  # no redo 
                    redo = False

                # P1 worse than P0, and P2 better than P0 (b1=0 and b2=1)
                #----------------------------------------
                if (b1 == False) and b2: 
                    # Can redo, moving to P1 
                    if redo1: 
                        newx,newy = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
                    # can't redo, increment and skip 
                    else: 
                        newx,newy = gincrement(x,y,xsgn=xsgn,  # to P1 
                                                   ysgn=ysgn,xr=xr,yr=yr)
                        redo = False
                        skip = True

                # P2 worse than P0, and P1 better than P0 (b1=1 and b2=0)
                #----------------------------------------
                if b1 and (b2 == False): 
                    # Can redo, moving to P2 
                    if redo2: 
                        newx,newy = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr,p2=True)
                    # can't redo, increment and skip 
                    else: 
                        newx,newy = gincrement(x,y,xsgn=xsgn,   # to P1 
                                                   ysgn=ysgn,xr=xr,yr=yr)
                        redo = False
                        skip = True


                # Both worse than P0
                #-------------------
                if (b1 == False) and (b2 == False):  # both bad, find worst 
                    # Can redo either one 
                    if redo1 and redo1: 
                        b12 = gbetter(par1,rms1,noise1,par2,rms2,noise2) 
                        # moving to P1 (P1 worse than P2) 
                        if b12:  # to P1 
                            newx,newy = gincrement(x,y,xsgn=xsgn,
                                                       ysgn=ysgn,xr=xr,yr=yr)
                        # moving to P2 (P2 worse than P1) 
                        if (b12 == False):  # to P1 
                            newx,newy = gincrement(x,y,xsgn=xsgn,
                                                        ysgn=ysgn,xr=xr,yr=yr,p2=True)

                    # Can't redo P2, goto P1 
                    if redo1 and (redo2 == False): 
                        newx,newy = gincrement(x,y,xsgn=xsgn,   # to P1 
                                                   ysgn=ysgn,xr=xr,yr=yr)
                    # Can't redo P1, goto P2 
                    if (redo1 == False) and redo2: 
                        newx,newy = gincrement(x,y,xsgn=xsgn,   # to P2 
                                                   ysgn=ysgn,xr=xr,yr=yr,p2=True)
                    # Can't do either, increment and skip 
                    if (redo1 == False) and (redo2 == False): 
                        newx,newy = gincrement(x,y,xsgn=xsgn,   # to P1 
                                                   ysgn=ysgn,xr=xr,yr=yr)
                        redo = False
                        skip = True 


                # Both better, increment and skip
                #--------------------------------
                if b1 and b2: 
                    newx,newy = gincrement(x,y,xsgn=xsgn,    # to P1 
                                               ysgn=ysgn,xr=xr,yr=yr)
                    redo = False
                    skip = True 


                x = newx 
                y = newy 

                # Getting the guess 
                if redo:  # redo 
                    # Getting the new guess from backward positions 
                    gguesspar,guessx,guessy = gguess(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)


            # Neither has been visited before, increment
            #===========================================
            if (p1 == False) and (p2 == False): 
                # Increment 
                newx,newy = gincrement(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)
                x = newx 
                y = newy 
                # Getting the guess 
                guesspar,guessx,guessy = gguess(x,y,xsgn=xsgn,ysgn=ysgn,xr=xr,yr=yr)


    return newx,newy,guesspar



def driver(datacube,xstart=0,ystart=0,outfile=None,noprint=False,noplot=False,
           plotxr=None,xsgn=1,ysgn=1,xr=None,yr=None,trackplot=False,
           noback=False,backret=True,gstruc=None,btrack=None,savestep=None,
           gassnum=None,subcube=None,wander=False,clobber=False):
    """
    This program runs the gaussian fitting program 
    on a large part of the HI all sky survey 
     
    This program can be run in three different modes. 
    /BACKRET: The program is allowed to go backwards, but MUST return to 
              the position it came from.  This seems to be the best mode 
              to run the program in for most purposes. 
              This is now the DEFAULT mode. 
    /NOBACK: The program is not allowed to go backwards.  It will always 
              move in a forward direction until it is finished.  Therefore 
              there is essentially no re-decomposition of positions. 
              To use this mode set /NOBACK on the command line. 
    /WANDER: The program is allowed to go backwards and forwards (Haud's 
              original algorithm).  If it goes back in y and only 
              the forward x position has been visited before then 
              it will need to go through all x positions before 
              returning to where it came from.  So either use strips 
              narrow in x or run this only once an initial solution 
              for every position has been found. 
              To use this mode set /WANDER on the command line. 
     
    If you input gstruc and btrack from a previous run, then the program 
    will start where it left off.  This is useful when you ran it with 
    /BACKRET and then want to use the /WANDER mode.  But currently it 
    start where it left off, so if you want it to start at the very 
    beginning again you need to do something (such as adding btrack(0) 
    to the end of btrack: btrack = [btrack,btrack(0)]. 
    This can also be used for debugging purposes. 
     
    Parameters
    ----------
    datacube      Cube object or filename.
    xstart        The x to start with 
    ystart        The y to start with 
    xr            X range 
    yr            Y range 
    xsgn          Direction of x increments (-1 or 1) 
    ysgn          Direction of y increments (-1 or 1) 
    outfile         File to save the structures to 
    plotxr          Plotting xrange 
    /trackplot      Track the progression visually 
    /noplot         Don't plot anything 
    /noprint        Don't print anything 
    /noback         The program is not allowed to go backwards 
    /backret        Any backwards motion must return to the position it 
                      came from 
    /wander         Allow backwards motion. Haud's algorithm 
    gstruc          Structure of gaussians 
    btrack          Tracking structure 
     
    Returns
    -------
    gstruc          This is a structure of all the gaussians found 
    btrack          Structure that keeps track of every move. 
     
    PROGRAMS USED 
    gaussfitter.pro fits an hi spectrum with gaussians, using Haud's method 
    gfit.pro        fits gaussians to spectrum, given initial parameters 
    gdev.pro        returns deviants from gfunc 
    gdev1.pro       returns deviants from gfunc1 
    gfunc.pro       computes gaussian from input parameters 
    gfunc1.pro      computer one gaussian and additional power terms from input pars 
    hinoise.pro     computes the noise in the HI spectrum 
    rdhispec.pro    reads an HI spectrum 
    setlimits.pro   sets default gaussian parameter limits 
    gest.pro        estimates the gaussian parameters at a certain point 
    gplot.pro       plots gaussians 
    gpeak.pro       finds peaks in spectrum 
    parcheck.pro    checks the gaussian parameters for problems 
    printgpar.pro   print the gaussian parameters 
    gremdup.pro     removes duplicate gaussian parameters 
    gremove.pro     removes bad gaussian parameters 
    gfind.pro       finds a position in the gaussian components structure 
    gbetter.pro     determines which of two gaussian fits is better 
    gguess.pro      get the guess for a forward movement 
    gincrement.pro  increment the position forward 
    gredo.pro       can a position be redone with a certain guess 
    gtrackplot.pro  track the progression visually 
    gsort.pro       sorts an array of gaussian parameters with 
                      decreasing area 
    gbtrack.pro     This increased the array size for par,guesspar if necessary 
     
    ANALYSIS TOOLS 
    ghess.pro       makes "Hess"-like diagrams of the gaussian data 
    gimage.pro      makes a 2D total intensity image from the gaussians 
    grecon.pro      reconstructs the datacube from the gaussians 
    greconplot.pro  reconstructs what the gaussian analysis did 
    gclean.pro      cleans out the zero-velocity region and the noise 
                       from the gstruc structure 
    gauss_plots.pro makes 6 diagnostic plots of the gaussian results 
    gauss_plots2-4  various plots of the gaussians 
     
    Created by David Nidever April 2005 
    """
 
    flag = 0 
    count = 0
    tstart = time.time() 
 
    # Setting parameters
    if xr is None:
         xr = [0.,2000.]
    if yr is None:
        yr = [0.,2000.]
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

    # Load the cube
    if type(datacube) is str:
        datacubefile = datacube
        print('Loading '+datacubefile)
        datacube = Cube.read(datacubefile)
        
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
            if type(btrack) is list:
                BTRACK += btrack
            else:
                BTRACK.append(btrack)
            if type(gstruct) is list:
                GSTRUC += gstruc
            else:
                GSTRUC.append(gstruc)
             
            count += 1 
            lastx = x 
            lastlast = y 
         
         
        # FIGURE OUT THE NEXT MOVE 
        #------------------------- 
        if (count > 0):
            import pdb; pdb.set_trace()
            lastx,lasty = x,y
            x,y,guesspar = nextmove(x,y,xsgn,ysgn,xr,yr)


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
 
        # Some bug checking 
        if x is None: 
            import pdb; pdb.set_trace() 
        if (x == lastx) and (y == lasty): 
            import pdb; pdb.set_trace() 
        if count != 0: 
            if (red1+red2+red3+red4 == 0) and redo: 
                import pdb; pdb.set_trace() 
 
 
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
            noise = spec.noise
            npts = spec.n
            #spec,v,lon,lat = loadspec(cubefile,x,y,npts=npts,noise=noise)
 
            # No good spectrum 
            if spec is None or np.sum(spec.flux)==0:
                rms = None
                noise = None
                skip = True
                count += 1
                continue
                #goto,SKIP 
 
            smspec = dln.savgol(spec.flux,21,2) 
            dum,vindcen = dln.closest(spec.vel,0) 
 
            # GETTIING THE VELOCITY RANGE around the zero-velocity MW peak 
            
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
            #par0,sigpar0,rms,noise,v2,spec2,resid2 = fitter.gaussfitter(spec,vmin=vmin,vmax=vmax,initpar=guesspar,
            #                                                            noprint=True,noplot=True)
            results = fitter.gaussfitter(spec,vmin=vmin,vmax=vmax,initpar=guesspar,noprint=True,noplot=True)            
 
            # FIT WITH NO GUESS (if first time and previous fit above with guess) 
            tp0 = gfind(x,y,xr=xr,yr=yr) 
            if (tp0 == 0) and (guesspar is not None):
                tpar0,tsigpar0,trms,noise,v2,spec2,resid2 = fitter.gaussfitter(v,spec,vmin=vmin,vmax=vmax,noprint=True,noplot=True)
                b = gbetter(par0,rms,noise,tpar0,trms,noise) 
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
                inv =v 
 
                
            # RUNNING GAUSSFITTER ON EVERYTHING WITHOUT THE ZERO-VELOCITY REGION, WITH GUESS 
            par0,sigpar0,rms,noise,v3,spec3,resid3 = fitter.gaussfitter(inv,inspec,initpar=guesspar2,noprint=True,noplot=True)
            
 
            # FIT WITH NO GUESS (if first time and previous fit above with guess) 
            if (tp0 == 0) and (len(guesspar) > 1):
                tpar0,tsigpar0,trms,noise,v3,spec3,resid3 = fitter.gaussfitter(inv,inspec,noprint=True,noplot=True)
                b = gbetter(par0,rms,noise,tpar0,trms,noise) 
                # The fit without the guess is better 
                if (b == 1): 
                    par0 = tpar0 
                    sigpar0 = tsigpar0 
                    rms = trms 
 
 
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
                
            print('fitting ',time.time()-t0)
 
            # PLOTTING/PRINTING, IF THERE WAS A FIT 
            if tstr['par'] is not None:
                # Getting the rms of all the components of the whole spectrum 
                th = gfunc(v,tstr['par'])
                rms = np.std(spec-th) 
                tstr['rms'] = rms 
 
                # Printing and plotting
                if noplot == False:
                    utils.gplot(v,spec,tstr.par,xlim=plotxr)
                if noprint == False:
                    utils.printgpar(tstr['par'],tstr['sigpar'],
                                    len(tstr['par'])//3,tstr['rms'],tstr['noise'])
                if trackplot:
                    utils.gtrackplot(x,y,lastx,lasty,redo, count,xr=xr,yr=yr,pstr=pstr,xstr=xstr,ystr=ystr)
            else:
                if noprint == False:
                    print('No gaussians found at this position!')

 
                
            # ADDING SOLUTION TO GSTRUC 
            if count == 0 : 
                tstr = gstruc_add(tstr)
                if count > 0: 
                    old = gfind(x,y,pind=pind1,rms=rms1,par=par1,noise=noise1,xr=xr,yr=yr) 
 
                # This is a re-decomposition 
                if old and redo: 
                    # Checking the two decompositions 
                    par2 = tstr['par']  # new one 
                    rms2 = tstr['rms'][0]
                    b = gbetter(par2,rms2,noise2,par1,rms1,noise1) 
                    # This one's better 
                    if (b == False): 
                        pind1 = gstruc_remove(pind1)  # removing the old 
                        tstr = gstruc_add(tstr)       # putting in the new 
                        t1 = time.time() 
                        print(time.time()-t1)
                        redo_fail = False
                    else: # re-decomposition failed 
                        redo_fail = False 
                        print('REDO FAILED!')
 
                # This is NOT a re-decomposition, add it 
                if (old == False) or (redo == False): 
                    t1 = time.time() 
                    tstr = gstruc_add(tstr)
                    print('gstruc ',time.time()-t1)
                    redo_fail = False
 
        # SKIP FITTING PART
        else: 
            # Creating a dummy structure 
            tstr = gstruc_dict.copy()
            redo_fail = False
            redo = False
            back = False
 
            if trackplot:
                utils.gtrackplot(x,y,lastx,lasty,redo,count,xr=xr,yr=yr,pstr=pstr,xstr=xstr,ystr=ystr)
 
 
        # FINISHING UP THE TRACKING STRUCTURE 
        npar = len(tstr['par'])
 
        track['par'] = tstr['par']
        track['rms'] = rms 
        track['noise'] = noise 
        track['redo_fail'] = redo_fail 
 
        # UPDATING THE TRACKING STRUCTURE
        BTRACK += [track]
 
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


