#!/usr/bin/env python

import os
import time
import numpy as np
from . import utils,fitter

def gincrement(lon,lat,newlon,newlat,lonr=None,latr=None,
               lonsgn=1,latsgn=1,nstep=1,p2=False):
    """
    This program increments the position 
     
    Parameters
    ----------
    lon      Longitude of current position 
    lat      Latitude of current position 
    lonsgn   Sign of longitude increment 
    latsgn   Sign of latitude increment 
    lonr     Two element array of longitude limits 
    latr     Two element array of latitude limits 
    nstep    Number of steps to increment 
    p2      Increment in latitude rather than in longtide 
     
    Returns
    -------
    newlon   Longitude of new position 
    newlat   Latitude of new position 
     
    When this program has problems it returns: 
    newlon = 999999. 
    newlat = 999999. 
     
    Created by David Nidever April 2005 
    """
 
    step = 1.0   # 0.5 
     
    # bad until proven okay 
    newlon = 999999. 
    newlat = 999999. 

    if lonr is None:
        lonr = [0.,2000.] 
    if latr is None: 
        latr = [0.,2000.] 
     
    # Longitude and latitude limits 
    lon0 = float(lonr[0]) 
    lon1 = float(lonr[1]) 
    lat0 = float(latr[0]) 
    lat1 = float(latr[1]) 
     
    # Are we in the box? 
    if (lon < lon0) or (lon > lon1) or (lat < lat0) or (lat > lat1): 
        newlon = 999999. 
        newlat = 999999. 
        goto,BOMB 
     
    # Wrong signs 
    if np.abs(lonsgn) != 1 : 
        goto,BOMB 
    if np.abs(latsgn) != 1 : 
        goto,BOMB 
     
    # figuring out the case 
    #bit = ((latsgn+1.)*0.5) + 2.*((lonsgn+1.)*0.5) 
    #bit = long(bit) 
    # 
    # bit (lon,lat) 
    # 0 - (-1,-1) 
    # 1 - (-1,+1) 
    # 2 - (+1,-1) 
    # 3 - (+1,+1) 
     
    tlon = float(lon) 
    tlat = float(lat) 
     
    # Looping through all the steps 
    for i in range(nstep): 
     
        # p2, incrementing vertically 
        if p2:
            # UP, at the end 
            if (latsgn == 1) and (tlat == lat1): 
                newlon = 999999. 
                newlat = 999999. 
                continue
            # DOWN, at the end 
            if (latsgn == -1) and (tlat == lat0): 
                newlon = 999999. 
                newlat = 999999. 
                continue 
         
            # Not at end, normal increment 
            newlon = tlon 
            newlat = tlat + latsgn * step 
         
        # Incrementing Sideways 
        else: 
            # RIGHT, lonsgn = +1 
            if (lonsgn == 1): 
                # UP, the very end 
                if (latsgn == 1) and (tlon == lon1) and (tlat == lat1): 
                    newlon = 999999. 
                    newlat = 999999. 
                    continue
                # DOWN, the very end 
                if (latsgn == -1) and (tlon == lon1) and (tlat == lat0): 
                    newlon = 999999. 
                    newlat = 999999. 
                    continue 
             
                # At end of longitude, increment latitude 
                if (tlon == lon1): 
                    newlon = lon0 
                    newlat = tlat + latsgn * step 
             
                # Normal increment 
                if (tlon != lon1): 
                    newlon = tlon + lonsgn * step 
                    newlat = tlat 
         
            # LEFT, lonsgn = -1 
            if (lonsgn == -1): 
                # UP, the very end 
                if (latsgn == 1) and (tlon == lon0) and (tlat == lat1): 
                    newlon = 999999. 
                    newlat = 999999. 
                    continue
                # DOWN, the very end 
                if (latsgn == -1) and (tlon == lon0) and (tlat == lat0): 
                    newlon = 999999. 
                    newlat = 999999. 
                    continue 
             
                # At end of longitude, increment latitude 
                if (tlon == lon0): 
                    newlon = lon1 
                    newlat = tlat + latsgn * step 
             
                # Normal increment 
                if (tlon != lon0): 
                    newlon = tlon + lonsgn * step 
                    newlat = tlat 
     
        # In case we're looping 
        tlon = newlon 
        tlat = newlat 
     
    # Final answer 
    newlon = tlon 
    newlat = tlat 

    return newlon,newlat


def gredo(lon,lat,guesslon,guesslat,guesspar):
    """
    This function checks wether we can redo this location again 
    returns 1 if this redo is okay 
    returns 0 if this redo is NOT okay 
    
    Parameters
    ----------
    lon       Longitude of current position 
    lat       Latitude of current position 
    guesslon  Longitude of guess position 
    guesslat  Latitude of guess position 
    guesspar  Guess parameters 
     
    Parameters
    ----------
    flag      Function value is either: 
                  1 - Redo okay. This guess has not been done before 
                  0 - Redo NOT okay. This guess has been done before 
     
    When there are any problems in this program it returns: 
    flag = -1 
     
    Created by David Nidever April 2005 
    """
     
    # Bad to start out with 
    flag = -1 
     
    # Making sure it's the right structure 
    tags = tag_names(*(!btrack.data)) 
    if (len(tags) != 15) : 
        return -1 
    btags = ['COUNT','LON','LAT','RMS','NOISE','PAR','GUESSPAR','GUESSLON','GUESSLAT',
             'BACK','REDO','REDO_FAIL','SKIP','LASTLON','LASTLAT'] 
    comp = (tags == btags) 
    if ((where(comp != 1))(0) != -1) : 
        return -1 
     
    flag = 1 # do it unless proven wrong 
     
    # FROM **ANY** PREVIOUS POSITION 
    prev, = np.where((*(!btrack.data)).lon == lon and (*(!btrack.data)).lat == lat)
    nprev = len(prev)
     
    gd1, = np.where(guesspar != 999999.) 
    if (len(gd1) == 0): 
        return 0 
    tguesspar = guesspar[gd1]
    nguesspar = len(tguesspar) 
    ngg = nguesspar//3 
     
    # FROM **ANY** PREVIOUS POSITION 
    # we have used a guess from this position before 
    #  but have the parameters changed sufficiently 
    if (nprev > 0): 
         
        # Looping through the previous ones 
        for i in range(nprev): 
            guesspar2 = (*(!btrack.data))[prev[i]].guesspar 
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
                        return 0
 
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


def gfind(lon,lat,lonr=None,latr=None):
    """
    This function helps find a latitude and longitude in 
    the gaussian components structure. 
     
    Parameters
    ----------
    lon     Longitude to search for 
    lat     Latitude to search for 
    lonr    Two element array of longitude limits, lonr=[lonmin,lonmax] 
    latr    Two element array of latitude limits, latr=[latmin,latmax] 
     
    Returns
    -------
    flag    The function value is either 0 or 1: 
                 1 - the position exists in the structure 
                 0 - the position does NOT exist in the structure 
                -1 - any problems 
    pind     Index of position in !gstruc. 
    rms     RMS of gaussian fit at the desired position 
    noise   Noise level at desired position 
    par     Parameters of gaussians in gstruc with the desired position 
     
    When there are any problems in this program it returns: 
     flag = -1 
     rms = -1 
     noise = -1 
     par = -1 
     pind = -1 
     
    Created by David Nidever April 2005 
    """

    # Assume bad until proven otherwise 
    flag = -1 
    rms = -1 
    noise = -1 
    par = -1 
    pind = -1 
     
    # Setting the ranges
    if lonr is not None:
        lon0 = lonr[0] 
        lon1 = lonr[1] 
    else: 
        lon0 = 0. 
        lon1 = 359.5
    if latr is not None:
        lat0 = latr[0] 
        lat1 = latr[1] 
    else: 
        lat0 = -90. 
        lat1 = 90. 
     
    if (lon < lon0) or (lon > lon1) or (lat < lat0) or (lat > lat1): 
        flag = -1 
        rms = -1 
        noise = -1 
        par = -1 
        pind = -1 
        return flag,pind,rms,noise,par
     
    # No !gstruc yet, first position 
    DEFSYSV,'!gstruc',exists=gstruc_exists 
    if gstruc_exists == 0: 
        rms = -1 
        noise = -1 
        par = -1 
        pind = -1 
        flag = 0 
        return flag,pind,rms,noise,par 
     
    # Looking for the position 
    t0 = time.time() 
    # LONSTART/LATSTART has a value for each position, faster searching 
    #  use NGAUSS and INDSTART to get the indices into DATA 
    pind, = np.where(*(!gstruc.lonstart) == lon and *(!gstruc.latstart) == lat)
    npind = len(pind)
    print('find ',time.time()-t0)
     
    # Found something, getting the values 
    if npind > 0: 
        ind = l64indgen((*(!gstruc.ngauss))[pind[0]])+(*(!gstruc.indstart))[pind[0]] 
        rms = first_el((*(!gstruc.data))[ind].rms) 
        noise = first_el((*(!gstruc.data))[ind].noise) 
        par = ((*(!gstruc.data))[ind].par)(*) 
        flag = 1 
         
    # Nothing found 
    else: 
        rms = -1 
        noise = -1 
        par = -1 
        flag = 0 
     
    return flag,pind,rms,noise,par


def gguess(lon,lat,lonsgn=1,latsgn=1,lonr=None,latr=None):
    """
    This program finds the best guess for the new profile 
    The one from the backwards positions. 
    If it can't find one then it returns 999999.'s 
     
    Parameters
    ----------
    lon      Longitude of current position 
    lat      Latitude of current position 
    lonsgn   Sign of longitude increment 
    latsgn   Sign of latitude increment 
    lonr     Two element array of longitude limits 
    latr     Two element array of latitude limits 
     
    Returns
    -------
    guesspar  The first guess gaussian parameters 
    guesslon  The longitude of the position where the 
                   guess parameters came from 
    guesslat  The latitude of the position where the 
                   guess parameters came from 
     
    When this program has problems it returns: 
      guesspar = 999999. 
      guesslon = 999999. 
      guesslat = 999999. 
     
    Created by David Nidever April 2005 
    """
     
    # Bad until proven okay 
    guesspar = 999999. 
    guesslon = 999999. 
    guesslat = 999999. 
     
    # Making sure it's the right structure 
    tags = tag_names(*(!gstruc.data)) 
    if (len(tags) != 6) : 
        return guesspar,guesslon,guesslat 
    comp = (tags == ['LON','LAT','RMS','NOISE','PAR','SIGPAR']) 
    if ((where(comp != 1))(0) != -1) :
        return guesspar,guesslon,guesslat         
     
    # Saving the originals 
    orig_lon = lon 
    orig_lat = lat 

    if lonr is None:
        lonr = [0.,359.5]
    if latr is None:
        latr = [-90.,90.] 
     
    # Is the longitude range continuous?? 
    if (lonr[0] == 0.) and (lonr[1] == 359.5): 
        cont = 1 
    else: 
        cont = 0 
     
    # getting the p3 and p4 positions 
    # P3 back in longitude (l-0.5), same latitude 
    # P4 back in latitude (b-0.5), same longitude 
    lon3,lat3 = gincrement(lon,lat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr)
    lon4,lat4 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
     
    # CHECKING OUT THE EDGES 
    # AT THE LEFT EDGE, and continuous, Moving RIGHT, use neighbor on other side 
    # Use it for the guess, but never move to it directly 
    if (lon == lonr[0]) and (lonsgn == 1) and (cont == 1): 
        lat3 = lat 
        lon3 = lonr[1] 
     
    # AT THE RIGHT EDGE, and continuous, Moving LEFT 
    if (lon == lonr[1]) and (lonsgn == -1) and (cont == 1): 
        lat3 = lat 
        lon3 = lonr[0] 
     
    # At the edge, NOT continuous, Moving RIGHT, NO P3 NEIGHBOR 
    if (lon == lonr[0]) and (lonsgn == 1) and (cont == 0): 
        lon3 = 999999. 
        lat3 = 999999. 
     
    # At the edge, NOT continuous, Moving LEFT, NO P3 NEIGHBOR 
    if (lon == lonr[1]) and (lonsgn == -1) and (cont == 0): 
        lon3 = 999999. 
        lat3 = 999999. 
     
    # Have they been visited before? 
    p3 = gfind(lon3,lat3,ind=ind3,rms=rms3,noise=noise3,par=par3) 
    p4 = gfind(lon4,lat4,ind=ind4,rms=rms4,noise=noise4,par=par4) 
     
    # Comparing the solutions 
    b34 = gbetter(par3,rms3,noise3,par4,rms4,noise4) 
     
    # selecting the best guess 
    if (b34 == 0):# using P3 
        guesspar = par3 
        guesslon = lon3 
        guesslat = lat3 
    if (b34 == 1):# using P4 
        guesspar = par4 
        guesslon = lon4 
        guesslat = lat4 
    if (b34 == -1): 
        guesspar = 999999. 
        guesslon = 999999. 
        guesslat = 999999. 
     
    # Putting the originals back 
    lon = orig_lon 
    lat = orig_lat 

    return guesspar,guesslon,guesslat  


def gdriver(lonstart,latstart,cubefile=None,outfile=None,noprint=False,noplot=False,
            plotxr=None,lonsgn=1,latsgn=1,lonr=None,latr=None,trackplot=False,
            noback=False,backret=True,gstruc=None,btrack=None,savestep=None,
            gassnum=None,subcube=None,wander=False):
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
              original algorithm).  If it goes back in latitude and only 
              the forward longitude position has been visited before then 
              it will need to go through all longitude positions before 
              returning to where it came from.  So either use strips 
              narrow in longitude or run this only once an initial solution 
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
    lonstart        The longitude to start with 
    latstart        The latitude to start with 
    =cubefile       The filename of the main datacube. 
    lonr            Longitude range 
    latr            Latitude range 
    lonsgn          Direction of longitude increments (-1 or 1) 
    latsgn          Direction of latitude increments (-1 or 1) 
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
    if lonr is None:
         lonr = [0.,2000.]
    if latr is None:
        latr = [0.,2000.]
    if lonstart is None:
        if lonsgn == 1: 
            lonstart = lonr[0] 
        else: 
            lonstart = lonr[1]
    if latstart is None:
        if latsgn == 1 : 
            latstart = latr[0] 
        else: 
            latstart = latr[1]
    if wander:
        backret = False
    if noback:
        backret = False
 
    # No cube filename input 
    if cubefile is None
        print('Must input CUBEFILE')
        return 
 
    # No mode selected, using default mode (backret) 
    if (backret == False) and (noback == False) and (wander == False): 
        print('' )
        print('!!! WARNING !!!  NO MODE SELECTED  ->  USING DEFAULT (BACKRET) MODE' )
        print('')
        sleep(3) 
        backret=1 
 
    # Restore file 
    restore_file = repstr(outfile,'.fits','_restore.sav') 
 
    # checking the file 
    if not keyword_set(file): 
        date = strsplit(systime(0),/extract) 
        time = strsplit(date(3),':',/extract) 
        # gauss_Apr202005_080136.dat, day, month, year, hour, minute, second 
        outfile = 'gauss_'+date(1)+date(2)+date(4)+'_'+time(0)+time(1)+time(2)+'.dat' 
    dum = findfile(outfile) 
    if dum != '': 
        print('THE FILE ',outfile,' EXISTS ALREADY !!!' )
        print('DO YOU WANT TO CONTINUE?')
        quest='' 
        read,quest 
        if quest != 'y' and quest != 'yes' and quest != 'YES' and     quest != 'yes' and quest != 'Yes' : 
        return
 
    # Printing out the inputs 
    print(' RUNNING GAUSSIAN ANALYSIS WITH THE FOLLOWING PARAMETERS')
    print('-----------------------------------------------------------')
    print(' STARTING POSITION = (',stringize(lonstart,ndec=1),',',stringize(latstart,ndec=1),')')
    print(' LONGITUDE RANGE = [',stringize(lonr(0),ndec=1),',',stringize(lonr(1),ndec=1),']' )
    print(' LATITUDE RANGE = [',stringize(latr(0),ndec=1),',',stringize(latr(1),ndec=1),']' )
    print(' LON DIRECTION = ',stringize(lonsgn))
    print(' LAT DIRECTION = ',stringize(latsgn))
    print(' FILE = ',file)
    print('-----------------------------------------------------------')
    if (backret == 1) : 
        print(' USING (BACKRET) MODE')
    if (noback == 1) : 
        print(' USING (NOBACK) MODE')
    if (wander == 1) : 
        print(' USING (WANDER) MODE')
    print('-----------------------------------------------------------')
    print('')
 
    # Is the longitude range continuous?? 
    if (lonr(0) == 0.) and (lonr(1) == 410.0): 
        cont = 1 
    else: 
        cont = 0 

    # Initializing some parameters 
    redo_fail = False 
    redo = False
    back = False
    lastlon = 999999. 
    lastlat = 999999. 
    p0 = False
    p1 = False
    p2 = False
    p3 = False
    p4 = False
     
    # Where are we starting 
    lon = lonstart 
    lat = latstart 
    #lon = lonr(0) 
    #lat = latr(0) 
     
    np = 99 
    btrack_schema = {count:999999.,lon:999999.,lat:999999.,rms:999999.,noise:999999.,par:fltarr(np)+999999,
                     guesspar:fltarr(np)+999999.,guesslon:999999.,guesslat:999999.,back:999999.,redo:999999.,
                     redo_fail:999999.,skip:999999.,lastlon:999999.,lastlat:999999.} 
     
    gstruc_schema = {lon:999999.,lat:999999.,rms:999999.,noise:999999.,
                     par:fltarr(3)+999999.,sigpar:fltarr(3)+999999.,glon:999999.,glat:999999.} 
     
    # STARTING THE LARGE LOOP 
    while (flag == 0): 
         
        t00 = time.time() 
         
        # P0 is the current position 
        # P1 forward in longitude (l+0.5), same latitude 
        # P2 forward in latitude (b+0.5), same longitude 
        # P3 back in longitude (l-0.5), same latitude 
        # P4 back in latitude (b-0.5), same longitude 
        # 
        # Move forward in longitude if possible 
         
        tstr = None
        tstr1 = None
        tstr2 = None
        skip = False
        guesslon = 999999. 
        guesslat = 999999. 
        guesspar = 999999. 
         
        #stop 
         
        # STARTING WITH BTRACK, RESTORING THE LAST STATE 
        if (count == 0) and (gstruc is not None and btrack is not None):
            nbtrack = len(btrack) 
            count = btrack[nbtrack-1].count 
            lon = btrack[nbtrack-1].lon 
            lat = btrack[nbtrack-1].lat 
            rms = btrack[nbtrack-1].rms 
            noise = btrack[nbtrack-1].noise 
            par = btrack[nbtrack-1].par 
            guesspar = btrack[nbtrack-1].guesspar 
            guesslon = btrack[nbtrack-1].guesslon 
            guesslat = btrack[nbtrack-1].guesslat 
            back = btrack[nbtrack-1].back 
            redo = btrack[nbtrack-1].redo 
            redo_fail = btrack[nbtrack-1].redo_fail 
            skip = btrack[nbtrack-1].skip 
            lastlon = btrack[nbtrack-1].lastlon 
            lastlat = btrack[nbtrack-1].lastlat 
            btrack = btrack_add(brack)
            gstruct = gstruc_add(gstruc)
             
            # REALLY NEED !GSTRUC AND !BTRACK 
             
            count = count+1 
            lastlon = lon 
            lastlast = lat 
         
         
        # FIGURE OUT THE NEXT MOVE 
        #------------------------- 
        if (count > 0): 
            # Get the positions, THIS IS THE PROPER WAY TO DO IT!!!!! 
            lon1,lat1 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
            lon2,lat2 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,p2=True) 
            lon3,lat3 = gincrement(lon,lat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr)
            lon4,lat4 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
             
            # Have they been visited before? 
            p0 = gfind(lon,lat,rms=rms0,noise=noise0,par=par0,lonr=lonr,latr=latr) 
            p1 = gfind(lon1,lat1,rms=rms1,noise=noise1,par=par1,lonr=lonr,latr=latr) 
            p2 = gfind(lon2,lat2,rms=rms2,noise=noise2,par=par2,lonr=lonr,latr=latr) 
            p3 = gfind(lon3,lat3,rms=rms3,noise=noise3,par=par3,lonr=lonr,latr=latr) 
            p4 = gfind(lon4,lat4,rms=rms4,noise=noise4,par=par4,lonr=lonr,latr=latr) 
             
            # PRINTING OUT SOME RELEVANT INFORMATION HERE 
            # comparing them 
            strb1 = gbetter(par0,rms0,noise0,par1,rms1,noise1) 
            strb2 = gbetter(par0,rms0,noise0,par2,rms2,noise2) 
            strb3 = gbetter(par0,rms0,noise0,par3,rms3,noise3) 
            strb4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 
             
            # do we need to redo? 
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

            strlon1 = '%5.1f' % lon1
            strlat1 = '%5.1f' % lat1
            strlon2 = '%5.1f' % lon2
            strlat2 = '%5.1f' % lat2
            strlon3 = '%5.1f' % lon3
            strlat3 = '%5.1f' % lat3
            strlon4 = '%5.1f' % lon4
            strlat4 = '%5.1f' % lat4
            if (lon1 == 999999.): strlon1 = '-----'
            if (lat1 == 999999.): strlat1 = '-----' 
            if (lon2 == 999999.): strlon2 = '-----'
            if (lat2 == 999999.): strlat2 = '-----' 
            if (lon3 == 999999.): strlon3 = '-----'
            if (lat3 == 999999.): strlat3 = '-----' 
            if (lon4 == 999999.): strlon4 = '-----'
            if (lat4 == 999999.): strlat4 = '-----' 
 
            # printing out the info 
            print('Count = %d' % count)
            print('Last/Current Position = (%.1f,%.2f)' %(lon,lat))
            print('Neighbors (position)  visited  better  redo')
            print('P1  (',strlon1,',',strlat1,')  ',p1, strb1, red1) 
            print('P2  (',strlon2,',',strlat2,')  ',p2, strb2, red2)
            print('P3  (',strlon3,',',strlat3,')  ',p3, strb3, red3)
            print('P4  (',strlon4,',',strlat4,')  ',p4, strb4, red4)
            print('')
 
            # if P3 or P4 worse than P0 then move back, to worst decomp 
            # if P3 and P4 better than P0 then move forward, if both 
            #  have been visited before then the worst decomp 
            #  if neither has been visited before then move to P1. 
 
 
            # If back redo and BACKRET=1 then return to pre-redo position 
            # This is done separately from the normal algorithm 
            if backret and back and redo: 
                back = False
                nbtrack = !btrack.count 
                newlon = (*(!btrack.data))[nbtrack-1].lastlon 
                newlat = (*(!btrack.data))[nbtrack-1].lastlat 
                lastlon = (*(!btrack.data))[nbtrack-1].lon 
                lastlat = (*(!btrack.data))[nbtrack-1].lat 
 
                # p0 is the redo position, p5 is the pre-redo position 
                p0 = gfind(lastlon,lastlat,rms=rms0,noise=noise0,par=par0,lonr=lonr,latr=latr) 
                p5 = gfind(newlon,newlat,rms=rms5,noise=noise5,par=par5,lonr=lonr,latr=latr) 
 
                b = gbetter(par0,rms0,noise0,par5,rms5,noise5) 
                redo = gredo(newlon,newlat,lastlon,lastlat,par0) 
                lon = newlon 
                lat = newlat 
 
                # back position better, redo pre-redo position 
                if (b == False) and redo: 
                    # Getting the guess 
                    guesspar,guesslon,guesslat = gguess(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
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
                    nbtrack = btrack.count 
                    lon = btrack['data'][nbtrack-1].lastlon 
                    lat = btrack['data'][nbtrack-1].lastlat 

                    
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
                        tnewlon,tnewlat = gincrement(lon,lat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr)
                        redo = gredo(tnewlon,tnewlat,lon,lat,par0)
                        
                        # P3 worse than P0, moving back
                        #------------------------------
                        if (b3 == False) and redo: 
                            newlon,newlat = gincrement(lon,lat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr)           
                            back = True   # moving backwards 
                            guesspar = par0 
                            guesslon = lon 
                            guesslat = lat 
                            lon = newlon 
                            lat = newlat 
                        else: 
                            back = False
                            redo = False
     
                    # Only P4 visited before
                    #=======================
                    if p4 and (p3 == False): 
                        b4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 
                        # Checking to see if this has been done before 
                        #   getting P4 position 
                        tnewlon,tnewlat = gincrement(lon,lat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
                        redo = gredo(tnewlon,tnewlat,lon,lat,par0)
                        
                        # P4 worse than P0, moving back
                        #------------------------------
                        if (b4 == False) and redo: 
                            newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
                            back = True  # moving backwards 
                            guesspar = par0 
                            guesslon = lon 
                            guesslat = lat 
                            lon = newlon 
                            lat = newlat 
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
                        tnewlon3,tnewlat3 = gincrement(lon,lat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr)
                        redo3 = gredo(tnewlon3,tnewlat3,lon,lat,par0) 
                        # Checking to see if this has been done before 
                        #   getting P4 position 
                        tnewlon4,tnewlat4 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
                        redo4 = gredo(tnewlon4,tnewlat4,lon,lat,par0) 
         
                        # P3 worse than P0, but P4 better than P0 (b3=0 and b4=1)
                        #----------------------------------------
                        if (b3 == False) and b4: 
                            # We can redo it, moving back to P3 
                            if (redo3 == 1): 
                                newlon,newlat = gincrement(lon,lat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr)
                            # Can't redo, move forward 
                            else: 
                                redo = False
                                back = False
         
                        # P4 worse than P0, but P3 better than P0 (b3=1 and b4=0)
                        #----------------------------------------
                        if b3 and b4: 
                            # We can redo it, moving back to P4 
                            if redo4: 
                                newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
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
                                    newlon,newlat = gincrement(lon,lat,lonsgn=-lonsgn,
                                                               latsgn=-latsgn,lonr=lonr,latr=latr)
                                # moving back to P4 (P4 worse than P3) 
                                if (b34 == False) :# to P4 
                                    newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,
                                                               latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
                            # Can't redo P4, goto P3 
                            if redo3 and (redo4 == False): 
                                newlon,newlat = gincrement(lon,lat,lonsgn=-lonsgn,    # to P3 
                                                           latsgn=-latsgn,lonr=lonr,latr=latr) 
                            # Can't redo P3, goto P4 
                            if (redo3 == False) and redo4: 
                                newlon,newlat = gincrement(lon,lat,newlon,newlat,lonsgn=lonsgn,   # to P4 
                                                           latsgn=-latsgn,lonr=lonr,latr=latr,p2=True)
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
                            guesslon = lon 
                            guesslat = lat 
                            lon = newlon 
                            lat = newlat 

                            
                #==============================
                # ---- CHECKING FORWARD ----
                #==============================
                if ((p3 == False) and (p4 == False)) or (back == False) or noback: 
     
                    # This is the very end 
                    if (lon1 == 999999.): 
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
                        tnewlon1,tnewlat1 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
                        redo1 = gredo(tnewlon1,tnewlat1,lon,lat,par0) 
     
                        # Moving to P1 (P1 worse than P0) 
                        if (b1 == False) and redo1: 
                            newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
                            lon = tnewlon1 
                            lat = tnewlat1 
                            # getting the guess 
                            guesspar,guesslon,guesslat = gguess(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
                        # Can't redo P1, or P1 better than P0, move another step ahead 
                        else: 
                            lon = tnewlon1 
                            lat = tnewlat1 
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
                        tnewlon1,tnewlat1 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
                        redo1 = gredo(tnewlon1,tnewlat1,lon,lat,par0) 
                        # Checking to see if this has been done before 
                        #   getting P2 position 
                        tnewlon2,tnewlat2 = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,p2=True)
                        redo2 = gredo(tnewlon2,tnewlat2,lon,lat,par0) 
                        if (redo1 == False) and (redo2 == False):  # no redo 
                            redo = False
     
                        # P1 worse than P0, and P2 better than P0 (b1=0 and b2=1)
                        #----------------------------------------
                        if (b1 == False) and b2: 
                            # Can redo, moving to P1 
                            if redo1: 
                                newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
                            # can't redo, increment and skip 
                            else: 
                                newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,  # to P1 
                                                           latsgn=latsgn,lonr=lonr,latr=latr)
                                redo = False
                                skip = True
     
                        # P2 worse than P0, and P1 better than P0 (b1=1 and b2=0)
                        #----------------------------------------
                        if b1 and (b2 == False): 
                            # Can redo, moving to P2 
                            if redo2: 
                                gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,/p2 
                            # can't redo, increment and skip 
                            else: 
                                newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,   # to P1 
                                                           latsgn=latsgn,lonr=lonr,latr=latr)
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
                                    newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,
                                                               latsgn=latsgn,lonr=lonr,latr=latr)
                                # moving to P2 (P2 worse than P1) 
                                if (b12 == False):  # to P1 
                                    newlon,newlat = gincrement,(on,lat,lonsgn=lonsgn,
                                                                latsgn=latsgn,lonr=lonr,latr=latr,p2=True)
         
                            # Can't redo P2, goto P1 
                            if redo1 and (redo2 == False): 
                                newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,   # to P1 
                                                           latsgn=latsgn,lonr=lonr,latr=latr)
                            # Can't redo P1, goto P2 
                            if (redo1 == False) and redo2: 
                                newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,   # to P2 
                                                           latsgn=latsgn,lonr=lonr,latr=latr,p2=True)
                            # Can't do either, increment and skip 
                            if (redo1 == False) and (redo2 == False): 
                                newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,   # to P1 
                                                           latsgn=latsgn,lonr=lonr,latr=latr)
                                redo = False
                                skip = True 
         

                        # Both better, increment and skip
                        #--------------------------------
                        if b1 and b2: 
                            newlon,newlat = gincrement(lon,lat,,lonsgn=lonsgn,    # to P1 
                                                       latsgn=latsgn,lonr=lonr,latr=latr)
                            redo = False
                            skip = True 

     
                        lon = newlon 
                        lat = newlat 
     
                        # Getting the guess 
                        if redo:  # redo 
                            # Getting the new guess from backward positions 
                            gguesspar,guesslon,guesslat = gguess(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)

 
                    # Neither has been visited before, increment
                    #===========================================
                    if (p1 === False) and (p2 == False): 
                        # Increment 
                        newlon,newlat = gincrement(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
                        lon = newlon 
                        lat = newlat 
                        # Getting the guess 
                        guesspar,guesslon,guesslat = gguess(lon,lat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr)
 
 
 
        # Starting the tracking structure, bad until proven good
        np = 99 #100 ;45 
        DEFSYSV,'!btrack',exists=btrack_exists 
        if btrack_exists == 1 :#100 ;45 
            np = len((*(!btrack.data))[0].par) > 99 
        track = {count:999999.,lon:999999.,lat:999999.,rms:999999.,noise:999999.,par:fltarr(np)+999999,
                 guesspar:fltarr(np)+999999.,guesslon:999999.,guesslat:999999.,back:999999.,redo:999999.,
                 redo_fail:999999.,skip:999999.,lastlon:999999.,lastlat:999999.} 
        nguesspar = len(guesspar) 
        track.count = count 
        track.lon = lon 
        track.lat = lat 
        track.lastlon = lastlon 
        track.lastlat = lastlat 
        track.guesspar(0:nguesspar-1) = guesspar 
        track.guesslon = guesslon 
        track.guesslat = guesslat 
        track.back = back 
        track.redo = redo 
        track.skip = skip 
 
        # Some bug checking 
        if lon == 999999.: 
            import pdb; pdb.set_trace() 
        if (lon == lastlon) and (lat == lastlat): 
            import pdb; pdb.set_trace() 
        if count != 0: 
            if (red1+red2+red3+red4 == 0) and (redo == 1): 
                import pdb; pdb.set_trace() 
 
 
        if skip: 
            print('SKIP')
 
        # FITTING THE SPECTRUM, UNLESS WE'RE SKIPPING IT 
        #------------------------------------------------ 
        if (skip == False: 
            t0 = time.time() 
 
            # Initial Printing 
            print('Fitting Gaussians to the HI spectrum at (',stringize(lon,ndec=1),',',stringize(lat,ndec=1),')')
            strout = ''
            if redo:
                strout = strout+'REDO '
            if back:
                strout = strout+'BACK'
            if back is False:
                strout = strout+'FORWARD' 
            print(strout) 
                          
            # Getting the HI spectrum 
            spec,v,glon,glat = loadspec(cubefile,lon,lat,npts=npts,noise=noise)
 
            # No good spectrum 
            if npts == 0: 
                rms = 999999. 
                noise = 999999. 
                skip = True
                continue
                #goto,SKIP 
 
            smspec = dln.savgol(spec,21,2) 
            dum,vindcen = dln.closest(v,0) 
 
            # GETTIING THE VELOCITY RANGE around the zero-velocity MW peak 
            
            # Finding the vel. low point 
            flag = 0 
            i = vindcen 
            while (flag == 0): 
                if smspec[i] <= noise: 
                    lo = i 
                if smspec[i] <= noise: 
                    flag = 1
                i -= 1 
                if i < 0: 
                    flag = 1 
            if len(lo) == 0:   # never dropped below the noise threshold 
                lo = 0 
            lo = np.maximum(0,(lo-20))
 
            # Finding the vel. high point 
            flag = 0 
            i = vindcen 
            while (flag == 0): 
                if smspec[i] <= noise : 
                    hi = i 
                if smspec[i] <= noise : 
                    flag = 1 
                i += 1 
                if i > npts-1: 
                    flag = 1 
            if len(hi) == 0: 
                hi = npts-1 
            hi = np.minimum((npts-1),(hi+20))
 
            vmin = v[lo] 
            vmax = v[hi] 
 
            # RUNNING GAUSSFITTER ON ZERO VELOCITY REGION, WITH GUESS 
            par0,sigpar0,rms,noise,v2,spec2,resid2 = fitter.gaussfitter(v,spec,vmin=vmin,vmax=vmax,initpar=guesspar,
                                                                        noprint=True,noplot=True)
 
            # FIT WITH NO GUESS (if first time and previous fit above with guess) 
            tp0 = gfind(lon,lat,lonr=lonr,latr=latr) 
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
                tstr1.lon = lon 
                tstr1.lat = lat 
                tstr1.glon = glon 
                tstr1.glat = glat 
                #tstr1.rms = rms 
                tstr1.noise = noise 
 
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
                        guesspar2 = 999999. 
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
                tstr2.lon = lon 
                tstr2.lat = lat 
                tstr2.glon = glon 
                tstr2.glat = glat 
                tstr2.noise = noise 
 
 
            # ADDING THE STRUCTURES TOGETHER, TSTR = [TSTR1,TSTR2] 
            if keyword_set(tstr1) and keyword_set(tstr2) : 
                tstr = [tstr1,tstr2] 
            if keyword_set(tstr1) and not keyword_set(tstr2) : 
                tstr = tstr1 
            if not keyword_set(tstr1) and keyword_set(tstr2) : 
                tstr = tstr2 
            if not keyword_set(tstr1) and not keyword_set(tstr2):# no gaussians 
                tstr = gstruc_schema 
                tstr.lon = lon 
                tstr.lat = lat 
                tstr.glon = glon 
                tstr.glat = glat 
                tstr.rms = rms 
                tstr.noise = noise 
                
            print('fitting ',time.time()-t0)
 
            # PLOTTING/PRINTING, IF THERE WAS A FIT 
            if tstr(0).par(0) != 999999.: 
                # Getting the rms of all the components of the whole spectrum 
                th = gfunc(v,(tstr.par)(*)) 
                rms = np.std(spec-th) 
                tstr.rms = rms 
 
                # Printing and plotting
                if noplot == False:
                    utils.gplot(v,spec,tstr.par,xlim=plotxr)
                if noprint == False:
                    utils.printgpar((tstr.par)(*),(tstr.sigpar)(*),
                                    len((tstr.par)(*))/3,first_el(tstr.rms),first_el(tstr.noise))
                if trackplot:
                    utils.gtrackplot(lon,lat,lastlon,lastlat,redo, count,lonr=lonr,latr=latr,pstr=pstr,xstr=xstr,ystr=ystr)
            else:
                if noprint == False:
                    print('No gaussians found at this position!')

 
                
            # ADDING SOLUTION TO GSTRUC 
            if count == 0 : 
                tstr = gstruc_add(tstr)
                if count > 0: 
                    old = gfind(lon,lat,pind=pind1,rms=rms1,par=par1,noise=noise1,lonr=lonr,latr=latr) 
 
                # This is a re-decomposition 
                if old and redo: 
                    # Checking the two decompositions 
                    par2 = (tstr.par)(*)# new one 
                    rms2 = first_el(tstr.rms) 
                    b = gbetter(par2,rms2,noise2,par1,rms1,noise1) 
                    # This one's better 
                    if (b == False): 
                        gstruc_remove,pind1  # removing the old 
                        gstruc_add,tstr  # putting in the new 
                        t1 = time.time() 
                        print(time.time()-t1)
                        redo_fail = False
                    else: # re-decomposition failed 
                        redo_fail = False 
                        print('REDO FAILED!')
 
                # This is NOT a re-decomposition, add it 
                if (old == False) or (redo == False): 
                    t1 = time.time() 
                    gstruc_add,tstr 
                    print('gstruc ',time.time()-t1)
                    redo_fail = False
 
        # SKIP FITTING PART
        else: 
 
            # creating a dummy structure 
            tstr = gstruc_schema 
            redo_fail = 0 
            redo = 0 
            back = 0 
 
            if trackplot:
                utils.gtrackplot(lon,lat,lastlon,lastlat,redo,count,lonr=lonr,latr=latr,pstr=pstr,xstr=xstr,ystr=ystr)
 
 
        # FINISHING UP THE TRACKING STRUCTURE 
        npar = len((tstr.par)(*)) 
 
        # Increasing btrack parameter arrays 
        if (count > 0): 
            nbpar = len((*(!btrack.data))[0].par) 
            if (npar > nbpar): 
            gbtrack,track,tstr 
 
        track.par[0:npar-1] = (tstr.par)(*) 
        track.rms = rms 
        track.noise = noise 
        track.redo_fail = redo_fail 
 
        # UPDATING THE TRACKING STRUCTURE 
        btrack_add,track 
 
        count += 1 
 
        # Saving the last position 
        lastlon = lon 
        lastlat = lat 
 
        print('this iteration ',time.time()-t00))
 
        # SAVING THE STRUCTURES, periodically
        if savestep == False:
            savestep = 50 
            nsave = savestep 
            if (int(count)//int(nsave) == long(count)/float(nsave)): 
                print('SAVING DATA!')
                #MWRFITS,(*(!gstruc.data))[0:!gstruc.count-1],file,/create 
                #MWRFITS,(*(!btrack.data))[0:!btrack.count-1],file,/silent# append 
                #gstruc = !gstruc & btrack = !btrack 
                #SAVE,gstruc,btrack,file=restore_file 
                #undefine,gstruc,btrack 
 
    # FINAL SAVE 
    print(str(len(gstruc),2),' final Gaussians')
    print('Saving data to ',file)
    #MWRFITS,(*(!gstruc.data))[0:!gstruc.count-1],file,/create 
    #MWRFITS,(*(!btrack.data))[0:!btrack.count-1],file,/silent# append 
    #gstruc = !gstruc & btrack = !btrack 
    #SAVE,gstruc,btrack,file=restore_file 
    #undefine,gstruc,btrack 
 
    print('dt = ',str(time.time()-tstart,2),' sec.')


