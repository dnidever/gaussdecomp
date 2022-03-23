#!/usr/bin/env python

import os
import time
import numpy as np

def gdriver,lonstart,latstart,cubefile=cubefile,file=file,noprint=noprint(noplot=noplot,            plotxr=plotxr,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,            trackplot=trackplot,noback=noback,backret=backret,gstruc=gstruc,            btrack=btrack,savestep=savestep,gassnum=gassnum,subcube=subcube 
     
    # This program runs the gaussian fitting program 
    # on a large part of the HI all sky survey 
    # 
    # This program can be run in three different modes. 
    # /BACKRET: The program is allowed to go backwards, but MUST return to 
    #          the position it came from.  This seems to be the best mode 
    #          to run the program in for most purposes. 
    #          This is now the DEFAULT mode. 
    # /NOBACK: The program is not allowed to go backwards.  It will always 
    #          move in a forward direction until it is finished.  Therefore 
    #          there is essentially no re-decomposition of positions. 
    #          To use this mode set /NOBACK on the command line. 
    # /WANDER: The program is allowed to go backwards and forwards (Haud's 
    #          original algorithm).  If it goes back in latitude and only 
    #          the forward longitude position has been visited before then 
    #          it will need to go through all longitude positions before 
    #          returning to where it came from.  So either use strips 
    #          narrow in longitude or run this only once an initial solution 
    #          for every position has been found. 
    #          To use this mode set /WANDER on the command line. 
    # 
    # If you input gstruc and btrack from a previous run, then the program 
    # will start where it left off.  This is useful when you ran it with 
    # /BACKRET and then want to use the /WANDER mode.  But currently it 
    # start where it left off, so if you want it to start at the very 
    # beginning again you need to do something (such as adding btrack(0) 
    # to the end of btrack: btrack = [btrack,btrack(0)]. 
    # This can also be used for debugging purposes. 
    # 
    # INPUT 
    #  lonstart        The longitude to start with 
    #  latstart        The latitude to start with 
    #  =cubefile       The filename of the main datacube. 
    #  lonr            Longitude range 
    #  latr            Latitude range 
    #  lonsgn          Direction of longitude increments (-1 or 1) 
    #  latsgn          Direction of latitude increments (-1 or 1) 
    #  file            File to save the structures to 
    #  plotxr          Plotting xrange 
    #  /trackplot      Track the progression visually 
    #  /noplot         Don't plot anything 
    #  /noprint        Don't print anything 
    #  /noback         The program is not allowed to go backwards 
    #  /backret        Any backwards motion must return to the position it 
    #                  came from 
    #  /wander         Allow backwards motion. Haud's algorithm 
    #  gstruc          Structure of gaussians 
    #  btrack          Tracking structure 
    # 
    # OUTPUT 
    #  gstruc          This is a structure of all the gaussians found 
    #  btrack          Structure that keeps track of every move. 
    # 
    # PROGRAMS USED 
    #  gaussfitter.pro fits an hi spectrum with gaussians, using Haud's method 
    #  gfit.pro        fits gaussians to spectrum, given initial parameters 
    #  gdev.pro        returns deviants from gfunc 
    #  gdev1.pro       returns deviants from gfunc1 
    #  gfunc.pro       computes gaussian from input parameters 
    #  gfunc1.pro      computer one gaussian and additional power terms from input pars 
    #  hinoise.pro     computes the noise in the HI spectrum 
    #  rdhispec.pro    reads an HI spectrum 
    #  setlimits.pro   sets default gaussian parameter limits 
    #  gest.pro        estimates the gaussian parameters at a certain point 
    #  gplot.pro       plots gaussians 
    #  gpeak.pro       finds peaks in spectrum 
    #  parcheck.pro    checks the gaussian parameters for problems 
    #  printgpar.pro   print the gaussian parameters 
    #  gremdup.pro     removes duplicate gaussian parameters 
    #  gremove.pro     removes bad gaussian parameters 
    #  gfind.pro       finds a position in the gaussian components structure 
    #  gbetter.pro     determines which of two gaussian fits is better 
    #  gguess.pro      get the guess for a forward movement 
    #  gincrement.pro  increment the position forward 
    #  gredo.pro       can a position be redone with a certain guess 
    #  gtrackplot.pro  track the progression visually 
    #  gsort.pro       sorts an array of gaussian parameters with 
    #                  decreasing area 
    #  gbtrack.pro     This increased the array size for par,guesspar if necessary 
    # 
    # ANALYSIS TOOLS 
    #  ghess.pro       makes "Hess"-like diagrams of the gaussian data 
    #  gimage.pro      makes a 2D total intensity image from the gaussians 
    #  grecon.pro      reconstructs the datacube from the gaussians 
    #  greconplot.pro  reconstructs what the gaussian analysis did 
    #  gclean.pro      cleans out the zero-velocity region and the noise 
    #                   from the gstruc structure 
    #  gauss_plots.pro makes 6 diagnostic plots of the gaussian results 
    #  gauss_plots2-4  various plots of the gaussians 
    # 
    # Created by David Nidever April 2005 
     
flag = 0 
count = 0LL 
tstart = time.time() 
 
# Setting parameters 
if len(lonr) == 0 : 
    lonr = [0.,2000.] 
if len(latr) == 0 : 
    latr = [0.,2000.] 
if len(lonsgn) == 0 : 
    lonsgn = 1. 
if len(latsgn) == 0 : 
    latsgn = 1. 
if len(lonstart) == 0 : 
    if lonsgn == 1 : 
        lonstart=lonr(0) 
    else: 
        lonstart=lonr(1) 
if len(latstart) == 0 : 
    if latsgn == 1 : 
        latstart=latr(0) 
    else: 
        latstart=latr(1) 
if len(noplot) == 0 : 
    noplot=0 
if len(noprint) == 0 : 
    noprint=0 
if len(noback) == 0 : 
    noback=0 
if len(backret) == 0 :# USING /BACKRET BY DEFAULT UNLESS NOT DESIRED 
    backret=1 
if len(wander) == 0 : 
    wander=0 
if keyword_set(wander) : 
    backret=0 
if keyword_set(noback) : 
    backret=0 
 
# No cube filename input 
if len(cubefile) == 0: 
    print('Must input CUBEFILE' 
    return 
 
# No mode selected, using default mode (backret) 
if (backret == 0) and (noback == 0) and (wander == 0): 
    print('' 
    print('!!! WARNING !!!  NO MODE SELECTED  ->  USING DEFAULT (BACKRET) MODE' 
    print('' 
    wait,3 
    backret=1 
 
# Restore file 
restore_file = repstr(file,'.fits','_restore.sav') 
 
# checking the file 
if not keyword_set(file): 
    date = strsplit(systime(0),/extract) 
    time = strsplit(date(3),':',/extract) 
    # gauss_Apr202005_080136.dat, day, month, year, hour, minute, second 
    file = 'gauss_'+date(1)+date(2)+date(4)+'_'+time(0)+time(1)+time(2)+'.dat' 
dum = findfile(file) 
if dum != '': 
    print('THE FILE ',file,' EXISTS ALREADY !!!' 
    print('DO YOU WANT TO CONTINUE?' 
    quest='' 
    read,quest 
    if quest != 'y' and quest != 'yes' and quest != 'YES' and     quest != 'yes' and quest != 'Yes' : 
        goto,BOMB 
 
# Printing out the inputs 
print(' RUNNING GAUSSIAN ANALYSIS WITH THE FOLLOWING PARAMETERS' 
print('-----------------------------------------------------------' 
print(' STARTING POSITION = (',stringize(lonstart,ndec=1),',',stringize(latstart,ndec=1),')' 
print(' LONGITUDE RANGE = [',stringize(lonr(0),ndec=1),',',stringize(lonr(1),ndec=1),']' 
print(' LATITUDE RANGE = [',stringize(latr(0),ndec=1),',',stringize(latr(1),ndec=1),']' 
print(' LON DIRECTION = ',stringize(lonsgn) 
print(' LAT DIRECTION = ',stringize(latsgn) 
print(' FILE = ',file 
print('-----------------------------------------------------------' 
if (backret == 1) : 
    print(' USING (BACKRET) MODE' 
if (noback == 1) : 
    print(' USING (NOBACK) MODE' 
if (wander == 1) : 
    print(' USING (WANDER) MODE' 
print('-----------------------------------------------------------' 
print('' 
 
# Is the longitude range continuous?? 
if (lonr(0) == 0.) and (lonr(1) == 410.0): 
    cont=1 
else: 
    cont=0 
#if (lonr(0) eq 0.) and (lonr(1) eq 359.5) then cont=1 else cont=0 
     
    # Initializing some parameters 
    redo_fail = 0 
    redo = 0 
    back = 0 
    lastlon = 999999. 
    lastlat = 999999. 
    p0 = 0 
    p1 = 0 
    p2 = 0 
    p3 = 0 
    p4 = 0 
    #if n_elements(gstruc) eq 0 then gstruc = 0 
    #if n_elements(btrack) eq 0 then btrack = 0 
     
    # Where are we starting 
    lon = lonstart 
    lat = latstart 
    #lon = lonr(0) 
    #lat = latr(0) 
     
    np = 99 
    btrack_schema = {count:999999.,lon:999999.,lat:999999.,rms:999999.,noise:999999.,par:fltarr(np)+999999,                 guesspar:fltarr(np)+999999.,guesslon:999999.,guesslat:999999.,back:999999.,redo:999999.,                 redo_fail:999999.,skip:999999.,lastlon:999999.,lastlat:999999.} 
     
    gstruc_schema = {lon:999999.,lat:999999.,rms:999999.,noise:999999.,                 par:fltarr(3)+999999.,sigpar:fltarr(3)+999999.,glon:999999.,glat:999999.} 
     
    # STARTING THE LARGE LOOP 
    WHILE (flag == 0): 
         
        t00 = time.time() 
         
        # P0 is the current position 
        # P1 forward in longitude (l+0.5), same latitude 
        # P2 forward in latitude (b+0.5), same longitude 
        # P3 back in longitude (l-0.5), same latitude 
        # P4 back in latitude (b-0.5), same longitude 
        # 
        # Move forward in longitude if possible 
         
        tstr = 0 
        tstr1 = 0 
        tstr2 = 0 
        skip = 0 
        guesslon = 999999. 
        guesslat = 999999. 
        guesspar = 999999. 
         
        #stop 
         
        # STARTING WITH BTRACK, RESTORING THE LAST STATE 
        if (count == 0) and (keyword_set(gstruc)) and (keyword_set(btrack)): 
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
            btrack_add,brack 
            gstruc_add,gstruc 
             
            # REALLY NEED !GSTRUC AND !BTRACK 
             
            count = count+1 
            lastlon = lon 
            lastlast = lat 
         
         
        # FIGURE OUT THE NEXT MOVE 
        #------------------------- 
        if (count > 0): 
             
             
            # ; Redo Failed! 
            # ; if we went forward then don't do anything, should continue forward 
            # ; If we went back and backret=1 then return to pre-redo position 
            # If (redo eq 1) and (redo_fail eq 1) and (back eq 1) then begin 
            #   ; Go back to pre-redo position 
            #   lon = lastlon 
            #   lat = lastlat 
            # endif 
             
             
            # Get the positions, THIS IS THE PROPER WAY TO DO IT!!!!! 
            gincrement,lon,lat,lon1,lat1,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr 
            gincrement,lon,lat,lon2,lat2,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,/p2 
            gincrement,lon,lat,lon3,lat3,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr 
            gincrement,lon,lat,lon4,lat4,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
             
            # have they been visited before? 
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
            if (p1 == 1) and (strb1 == 0): 
                red1=1 
            else: 
                red1=0 
            if (p2 == 1) and (strb2 == 0): 
                red2=1 
            else: 
                red2=0 
            if (p3 == 1) and (strb3 == 0): 
                red3=1 
            else: 
                red3=0 
            if (p4 == 1) and (strb4 == 0): 
                red4=1 
            else: 
                red4=0 
             
            strlon1 = stringize(lon1,ndec=1,length=5)  &  strlat1 = stringize(lat1,ndec=1,length=5) 
            strlon2 = stringize(lon2,ndec=1,length=5)  &  strlat2 = stringize(lat2,ndec=1,length=5) 
            strlon3 = stringize(lon3,ndec=1,length=5)  &  strlat3 = stringize(lat3,ndec=1,length=5) 
            strlon4 = stringize(lon4,ndec=1,length=5)  &  strlat4 = stringize(lat4,ndec=1,length=5) 
            if (lon1 == 999999.) : 
                strlon1 = '-----'  &  if (lat1 == 999999.) : 
                strlat1 = '-----' 
        if (lon2 == 999999.) : 
            strlon2 = '-----'  &  if (lat2 == 999999.) : 
            strlat2 = '-----' 
    if (lon3 == 999999.) : 
        strlon3 = '-----'  &  if (lat3 == 999999.) : 
        strlat3 = '-----' 
if (lon4 == 999999.) : 
    strlon4 = '-----'  &  if (lat4 == 999999.) : 
    strlat4 = '-----' 
 
# printing out the info 
print('Count = ',stringize(count) 
print('Last/Current Position = (',stringize(lon,ndec=1),',',stringize(lat,ndec=1),')' 
print('Neighbors (position)  visited  better  redo' 
print('P1  (',strlon1,',',strlat1,')  ',p1, strb1, red1 
print('P2  (',strlon2,',',strlat2,')  ',p2, strb2, red2 
print('P3  (',strlon3,',',strlat3,')  ',p3, strb3, red3 
print('P4  (',strlon4,',',strlat4,')  ',p4, strb4, red4 
print('' 
 
 
#p1 = gfind(gstruc,lon+lonsgn*0.5,lat,ind=ind1,rms=rms1,noise=noise1,par=par1,lonr=lonr,latr=latr) 
#p2 = gfind(gstruc,lon,lat+latsgn*0.5,ind=ind2,rms=rms2,noise=noise2,par=par2,lonr=lonr,latr=latr) 
#p3 = gfind(gstruc,lon-lonsgn*0.5,lat,ind=ind3,rms=rms3,noise=noise3,par=par3,lonr=lonr,latr=latr) 
#p4 = gfind(gstruc,lon,lat-latsgn*0.5,ind=ind4,rms=rms4,noise=noise4,par=par4,lonr=lonr,latr=latr) 
 
# What to do when P3 or P4 are -1 (i.e. they don't exist) ??? 
# Fixed 
 
# if P3 or P4 worse than P0 then move back, to worst decomp 
# if P3 and P4 better than P0 then move forward, if both 
#  have been visited before then the worst decomp 
#  if neither has been visited before then move to P1. 
 
 
# If back redo and BACKRET=1 then return to pre-redo position 
# This is done separately from the normal algorithm 
IF (backret == 1) and (back == 1) and (redo == 1): 
 
back = 0 
 
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
if (b == 0) and (redo == 1): 
     
    # getting the guess 
    gguess,lon,lat,guesspar,guesslon,guesslat,lonsgn=lonsgn,latsgn=latsgn,               lonr=lonr,latr=latr 
     
    redo = 1 
    skip = 0 
     
    # back position worse, or can't redo pre-position, skip 
else: 
     
    redo = 0 
    skip = 1 
 
#stop 
 
# NOT (redo back and /backret) 
# Starting Normal Algorithm 
else: 
 
 
# Redo Failed! 
# if we went forward then don't do anything, should continue forward 
# If we went back and backret=1 then return to pre-redo position 
If (redo == 1) and (redo_fail == 1) and (back == 1): 
    # Go back to pre-redo position 
    nbtrack = !btrack.count 
    lon = (*(!btrack.data))[nbtrack-1].lastlon 
    lat = (*(!btrack.data))[nbtrack-1].lastlat 
 
#---- CHECKING BACKWARDS ---- 
if ((p3 == 1) or (p4 == 1)) and (noback == 0): 
     
    # Only P3 visited before 
    if (p3 == 1) and (p4 != 1): 
        b3 = gbetter(par0,rms0,noise0,par3,rms3,noise3) 
         
        # checking to see if this has been done before 
        # getting P3 position 
        gincrement,lon,lat,tnewlon,tnewlat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr 
        redo = gredo(tnewlon,tnewlat,lon,lat,par0) 
         
        # P3 worse than P0, moving back 
        if (b3 == 0) and (redo == 1): 
            gincrement,lon,lat,newlon,newlat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr 
             
            back = 1# moving backwards 
            guesspar = par0 
            guesslon = lon 
            guesslat = lat 
            lon = newlon 
            lat = newlat 
        else: 
            back = 0 
            redo = 0 
         
        #stop 
     
     
    # Only P4 visited before 
    If (p4 == 1) and (p3 != 1): 
        b4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 
         
        # checking to see if this has been done before 
        # getting P4 position 
        gincrement,lon,lat,tnewlon,tnewlat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
        redo = gredo(tnewlon,tnewlat,lon,lat,par0) 
         
        # P4 worse than P0, moving back 
        if (b4 == 0) and (redo == 1): 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
             
            back = 1# moving backwards 
            guesspar = par0 
            guesslon = lon 
            guesslat = lat 
            lon = newlon 
            lat = newlat 
        else: 
            back = 0 
            redo = 0 
     
     
    # Both visited before 
    If (p3 == 1) and (p4 == 1): 
        b3 = gbetter(par0,rms0,noise0,par3,rms3,noise3) 
        b4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 
        redo = 1# redo unless proven otherwise 
         
        # checking to see if this has been done before 
        # getting P3 position 
        gincrement,lon,lat,tnewlon3,tnewlat3,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr 
        redo3 = gredo(tnewlon3,tnewlat3,lon,lat,par0) 
         
        # checking to see if this has been done before 
        # getting P4 position 
        gincrement,lon,lat,tnewlon4,tnewlat4,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
        redo4 = gredo(tnewlon4,tnewlat4,lon,lat,par0) 
         
        # P3 worse than P0, but P4 better than P0 
        if (b3 == 0) and (b4 == 1): 
             
            # we can redo it, moving back to P3 
            if (redo3 == 1): 
                gincrement,lon,lat,newlon,newlat,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr 
                # can't redo, move forward 
            else: 
                redo=0 
                back=0 
    # b3=0 and b4=1 
         
        # P4 worse than P0, but P3 better than P0 
        if (b3 == 1) and (b4 == 0): 
             
            # we can redo it, moving back to P4 
            if (redo4 == 1): 
                gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
                # can't redo, move forward 
            else: 
                redo=0 
                back=0 
    # b3=1 and b4=0 
         
        # both bad 
        if (b3 == 0) and (b4 == 0): 
             
            # Can redo either one 
            if (redo3 == 1) and (redo4 == 1): 
                b34 = gbetter(par3,rms3,noise3,par4,rms4,noise4) 
                 
                # moving back to P3 (P3 worse than P4) 
                if (b34 == 1) :# to P3 
                    gincrement,lon,lat,newlon,newlat,lonsgn=-lonsgn,$ 
                latsgn=-latsgn,lonr=lonr,latr=latr 
                 
                # moving back to P4 (P4 worse than P3) 
                if (b34 == 0) :# to P4 
                    gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$ 
                latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
             
            # Can't redo P4, goto P3 
            if (redo3 == 1) and (redo4 == 0): 
                gincrement,lon,lat,newlon,newlat,lonsgn=-lonsgn,$# to P3 
                latsgn=-latsgn,lonr=lonr,latr=latr 
             
            # Can't redo P3, goto P4 
            if (redo3 == 0) and (redo4 == 1): 
                gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$# to P4 
                latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
             
            # Can't do either, move forward 
            if (redo3 == 0) and (redo4 == 0): 
                redo = 0 
                back = 0 
             
    # both worse 
         
        # both are better than P0, move forward 
        if (b3 == 1) and (b4 == 1): 
            back = 0 
            redo = 0 
         
        # one is worse than P0 
        if (redo == 1): 
            back = 1# moving backwards 
            guesspar = par0 
            guesslon = lon 
            guesslat = lat 
            lon = newlon 
            lat = newlat 
         
# both p3 and p4 visited before 
     
# checking backwards 
 
 
# ---- CHECKING FORWARD ---- 
IF ((p3 != 1) and (p4 != 1)) or (back == 0) or (noback == 1): 
     
    # This is the very end 
    if (lon1 == 999999.): 
flag = 1 
    goto, BOMB 
 
back = 0# moving forward 
 
# Only P1 has been visited before 
If (p1 == 1) and (p2 != 1): 
    b1 = gbetter(par0,rms0,noise0,par1,rms1,noise1) 
    redo = 1 
     
    # checking to see if this has been done before 
    # getting P1 position 
    gincrement,lon,lat,tnewlon1,tnewlat1,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr 
    redo1 = gredo(tnewlon1,tnewlat1,lon,lat,par0) 
     
    # moving to P1 (P1 worse than P0) 
    if (b1 == 0) and (redo1 == 1): 
        gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr 
         
        lon = tnewlon1 
        lat = tnewlat1 
         
        # getting the guess 
        gguess,lon,lat,guesspar,guesslon,guesslat,lonsgn=lonsgn,latsgn=latsgn,                 lonr=lonr,latr=latr 
         
        # Can't redo P1, or P1 better than P0, move another step ahead 
    else: 
         
        lon = tnewlon1 
        lat = tnewlat1 
         
        redo = 0 
        skip = 1# don't fit this one 
     
# only p1 
 
 
# Only P2 has been visited before, THIS SHOULD NEVER HAPPEN 
If (p2 == 1) and (p1 != 1): 
    print('This should never happen!!' 
    import pdb; pdb.set_trace() 
# only p2 
 
 
# Both have been visited before 
If (p1 == 1) and (p2 == 1): 
    b1 = gbetter(par0,rms0,noise0,par1,rms1,noise1) 
    b2 = gbetter(par0,rms0,noise0,par2,rms2,noise2) 
    redo = 1# redo unless proven otherwise 
     
    # checking to see if this has been done before 
    # getting P1 position 
    gincrement,lon,lat,tnewlon1,tnewlat1,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr 
    redo1 = gredo(tnewlon1,tnewlat1,lon,lat,par0) 
     
    # checking to see if this has been done before 
    # getting P2 position 
    gincrement,lon,lat,tnewlon2,tnewlat2,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,/p2 
    redo2 = gredo(tnewlon2,tnewlat2,lon,lat,par0) 
     
    if (redo1+redo2 == 0) :# no redo 
        redo=0 
     
    # P1 worse than P0, and P2 better than P0 
    if (b1 == 0) and (b2 == 1): 
         
        # can redo, moving to P1 
        if (redo1 == 1): 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr 
            # can't redo, increment and skip 
        else: 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$# to P1 
            latsgn=latsgn,lonr=lonr,latr=latr 
            redo = 0 
            skip = 1 
# b1=0 and b2=1 
     
    # P2 worse than P0, and P1 better than P0 
    if ((b1 == 1) and (b2 == 0)): 
         
        # can redo, moving to P2 
        if (redo2 == 1): 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,/p2 
            # can't redo, increment and skip 
        else: 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$# to P1 
            latsgn=latsgn,lonr=lonr,latr=latr 
            redo = 0 
            skip = 1 
# b1=1 and b2=0 
     
    # both worse than P0 
    if (b1 == 0) and (b2 == 0):# both bad, find worst 
         
        # Can redo either one 
        if (redo1 == 1) and (redo1 == 1): 
            b12 = gbetter(par1,rms1,noise1,par2,rms2,noise2) 
             
            # moving to P1 (P1 worse than P2) 
            if (b12 == 1) :# to P1 
                gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$ 
            latsgn=latsgn,lonr=lonr,latr=latr 
             
            # moving to P2 (P2 worse than P1) 
            if (b12 == 0) :# to P1 
                gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$ 
            latsgn=latsgn,lonr=lonr,latr=latr,/p2 
         
        # Can't redo P2, goto P1 
        if (redo1 == 1) and (redo2 == 0): 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$# to P1 
            latsgn=latsgn,lonr=lonr,latr=latr 
         
        # Can't redo P1, goto P2 
        if (redo1 == 0) and (redo2 == 1): 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$# to P2 
            latsgn=latsgn,lonr=lonr,latr=latr,/p2 
         
        # Can't do either, increment and skip 
        if (redo1 == 0) and (redo2 == 0): 
            gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$# to P1 
            latsgn=latsgn,lonr=lonr,latr=latr 
            redo = 0 
            skip = 1 
         
# both worse 
     
    # both better, increment and skip 
    if (b1 == 1) and (b2 == 1): 
        gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,$# to P1 
        latsgn=latsgn,lonr=lonr,latr=latr 
        redo = 0 
        skip = 1 
# both better 
     
    lon = newlon 
    lat = newlat 
     
    # getting the guess 
    if (redo == 1):# redo 
        # getting the new guess from backward positions 
        gguess,lon,lat,guesspar,guesslon,guesslat,lonsgn=lonsgn,latsgn=latsgn,                 lonr=lonr,latr=latr 
# both 
 
 
# Neither has been visited before, increment 
If (p1 != 1) and (p2 != 1): 
     
    # increment 
    gincrement,lon,lat,newlon,newlat,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr 
     
    lon = newlon 
    lat = newlat 
     
    # getting the guess 
    gguess,lon,lat,guesspar,guesslon,guesslat,lonsgn=lonsgn,latsgn=latsgn,               lonr=lonr,latr=latr 
     
# neither 
 
# checking forwards 
 
# not (redo back and /backret) 
 
 
# Getting the Info to Print 
# 
# get the positions 
#gincrement,lon,lat,lon1,lat1,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr 
#gincrement,lon,lat,lon2,lat2,lonsgn=lonsgn,latsgn=latsgn,lonr=lonr,latr=latr,/p2 
#gincrement,lon,lat,lon3,lat3,lonsgn=-lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr 
#gincrement,lon,lat,lon4,lat4,lonsgn=lonsgn,latsgn=-latsgn,lonr=lonr,latr=latr,/p2 
# 
# have they been visited before? 
#p0 = gfind(gstruc,lon,lat,ind=ind0,rms=rms0,noise=noise0,par=par0,lonr=lonr,latr=latr) 
#p1 = gfind(gstruc,lon1,lat1,ind=ind1,rms=rms1,noise=noise1,par=par1,lonr=lonr,latr=latr) 
#p2 = gfind(gstruc,lon2,lat2,ind=ind2,rms=rms2,noise=noise2,par=par2,lonr=lonr,latr=latr) 
#p3 = gfind(gstruc,lon3,lat3,ind=ind3,rms=rms3,noise=noise3,par=par3,lonr=lonr,latr=latr) 
#p4 = gfind(gstruc,lon4,lat4,ind=ind4,rms=rms4,noise=noise4,par=par4,lonr=lonr,latr=latr) 
# 
# comparing them 
#strb1 = gbetter(par0,rms0,noise0,par1,rms1,noise1) 
#strb2 = gbetter(par0,rms0,noise0,par2,rms2,noise2) 
#strb3 = gbetter(par0,rms0,noise0,par3,rms3,noise3) 
#strb4 = gbetter(par0,rms0,noise0,par4,rms4,noise4) 
# 
# printing out the info 
#print,'Current Position = (',stringize(lon,ndec=1),',',stringize(lat,ndec=1),')' 
#print,'Neighbors 
#print,'P1  (',stringize(lon1,ndec=1),',',stringize(lat1,ndec=1),')  ',strb1 
#print,'P2  (',stringize(lon2,ndec=1),',',stringize(lat2,ndec=1),')  ',strb2 
#print,'P3  (',stringize(lon3,ndec=1),',',stringize(lat3,ndec=1),')  ',strb3 
#print,'P4  (',stringize(lon4,ndec=1),',',stringize(lat4,ndec=1),')  ',strb4 
 
# count ne 0 
 
 
# Starting the tracking structure, bad until proven good 
np = 99#100 ;45 
DEFSYSV,'!btrack',exists=btrack_exists 
if btrack_exists == 1 :#100 ;45 
np = len((*(!btrack.data))[0].par) > 99 
track = {count:999999.,lon:999999.,lat:999999.,rms:999999.,noise:999999.,par:fltarr(np)+999999,           guesspar:fltarr(np)+999999.,guesslon:999999.,guesslat:999999.,back:999999.,redo:999999.,           redo_fail:999999.,skip:999999.,lastlon:999999.,lastlat:999999.} 
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
 
 
# some bug checking 
if lon == 999999. : 
import pdb; pdb.set_trace() 
if (lon == lastlon) and (lat == lastlat) : 
import pdb; pdb.set_trace() 
if count != 0 : 
if (red1+red2+red3+red4 == 0) and (redo == 1) : 
import pdb; pdb.set_trace() 
#if count gt 40 then stop 
 
 
If (skip == 1) : 
print('SKIP' 
 
# FITTING THE SPECTRUM, UNLESS WE'RE SKIPPING IT 
#------------------------------------------------ 
If (skip != 1): 
t0 = time.time() 
 
# Initial Printing 
print('Fitting Gaussians to the HI spectrum at (',stringize(lon,ndec=1),',',stringize(lat,ndec=1),')' 
strout = '' 
if keyword_set(redo) : 
strout = strout+'REDO ' 
if keyword_set(back) : 
strout = strout+'BACK' 
if not keyword_set(back) : 
strout = strout+'FORWARD' 
print(strout 
 
# Getting the HI spectrum 
GLOADSPEC,cubefile,lon,lat,spec,v,glon,glat,npts=npts,noise=noise 
 
# No good spectrum 
if npts == 0: 
rms = 999999. 
noise = 999999. 
skip = 1 
goto,SKIP 
 
smspec = savgolsm(spec,[10,10,2]) 
dum = closest(0,v,ind=vindcen) 
 
# GETTIING THE VELOCITY RANGE around the zero-velocity MW peak 
 
# finding the vel. low point 
flag=0 
i = vindcen 
while (flag == 0): 
if smspec[i] <= noise : 
lo = i 
if smspec[i] <= noise : 
flag = 1 
i = i - 1 
if i < 0 : 
flag=1 
if len(lo) == 0 :# never dropped below the noise threshold 
lo=0 
lo = 0 > (lo-20) 
 
# finding the vel. high point 
flag=0 
i = vindcen 
while (flag == 0): 
if smspec[i] <= noise : 
hi = i 
if smspec[i] <= noise : 
flag = 1 
i = i + 1 
if i > npts-1 : 
flag=1 
if len(hi) == 0 : 
hi=npts-1 
hi = (npts-1) < (hi+20) 
 
vmin = v[lo] 
vmax = v[hi] 
 
# RUNNING GAUSSFITTER ON ZERO VELOCITY REGION, WITH GUESS 
gaussfitter,lon,lat,par0,sigpar0,rms,noise,v2,spec2,resid2,vmin=vmin,vmax=vmax,                /noprint(/noplot,inpar=guesspar,inv=v,inspec=spec 
 
# FIT WITH NO GUESS (if first time and previous fit above with guess) 
tp0 = gfind(lon,lat,lonr=lonr,latr=latr) 
if (tp0 == 0) and (len(guesspar) > 1): 
 
gaussfitter,lon,lat,tpar0,tsigpar0,trms,noise,v2,spec2,resid2,vmin=vmin,vmax=vmax,                  /noprint(/noplot,inv=v,inspec=spec 
 
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
tstr1(i).par = par0(3*i:3*i+2) 
for i in range(ngauss): 
tstr1(i).sigpar = sigpar0(3*i:3*i+2) 
tstr1.lon = lon 
tstr1.lat = lat 
tstr1.glon = glon 
tstr1.glat = glat 
#tstr1.rms = rms 
tstr1.noise = noise 
 
 
# REMOVING ZERO-VELOCITY parameters and spectrum 
if par0(0) != -1: 
th = gfunc(v,par0) 
inspec = spec-th 
inv = v 
npts = len(v) 
if len(guesspar) > 1: 
inpar1 = guesspar 
inpar2 = guesspar 
gremove,inpar1,v(0:lo),spec(0:lo) 
gremove,inpar2,v(hi:npts-1),spec(hi:npts-1) 
guesspar2 = [inpar1,inpar2] 
gd , = np.where(guesspar2 != -1,ngd) 
if (ngd == 0) :# no guess 
    guesspar2 = 999999. 
if (ngd != 0) : 
    guesspar2 = guesspar2(gd) 
else: 
inspec = spec 
inv =v 
 
 
# RUNNING GAUSSFITTER ON EVERYTHING WITHOUT THE ZERO-VELOCITY REGION, WITH GUESS 
gaussfitter,lon,lat,par0,sigpar0,rms,noise,v3,spec3,resid3,inv=inv,inspec=inspec,                /noprint(/noplot,inpar=guesspar2 
 
 
# FIT WITH NO GUESS (if first time and previous fit above with guess) 
if (tp0 == 0) and (len(guesspar) > 1): 
 
gaussfitter,lon,lat,tpar0,tsigpar0,trms,noise,v3,spec3,resid3,inv=inv,inspec=inspec,                  /noprint(/noplot 
 
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
tstr2(i).par = par0(3*i:3*i+2) 
for i in range(ngauss): 
tstr2(i).sigpar = sigpar0(3*i:3*i+2) 
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
 
print('fitting ',time.time()-t0 
 
# PLOTTING/PRINTING, IF THERE WAS A FIT 
if tstr(0).par(0) != 999999.: 
 
# getting the rms of all the components of the whole spectrum 
th = gfunc(v,(tstr.par)(*)) 
rms = np.std(spec-th) 
tstr.rms = rms 
 
# printing and plotting 
if not keyword_set(noplot) : 
gplot,v,spec,(tstr.par)(*),xrange=plotxr 
if not keyword_set(noprint) : 
printgpar,(tstr.par)(*),(tstr.sigpar)(*),         len((tstr.par)(*))/3,first_el(tstr.rms),first_el(tstr.noise) 
if keyword_set(trackplot) : 
gtrackplot,lon,lat,lastlon,lastlat,redo,                                        count,lonr=lonr,latr=latr,pstr=pstr,xstr=xstr,ystr=ystr 
 
else: 
if not keyword_set(noprint) : 
print('No gaussians found at this position!' 
# printing and plotting 
 
 
# ADDING SOLUTION TO GSTRUC 
If count == 0 : 
gstruc_add,tstr 
If count > 0: 
old = gfind(lon,lat,pind=pind1,rms=rms1,par=par1,noise=noise1,lonr=lonr,latr=latr) 
 
# This is a re-decomposition 
if (old == 1) and (redo == 1): 
 
# checking the two decompositions 
par2 = (tstr.par)(*)# new one 
rms2 = first_el(tstr.rms) 
b = gbetter(par2,rms2,noise2,par1,rms1,noise1) 
 
# This one's better 
if (b == 0): 
    gstruc_remove,pind1# removing the old 
    gstruc_add,tstr# putting in the new 
    #remove,ind1,gstruc      ; removing the old 
    t1 = time.time() 
    #gstruc = [gstruc,tstr]  ; putting in the new 
    print(time.time()-t1 
    redo_fail = 0 
else:# re-decomposition failed 
    redo_fail = 1 
    print('REDO FAILED!' 
 
# This is NOT a re-decomposition, add it 
if (old != 1) or (redo != 1): 
t1 = time.time() 
gstruc_add,tstr 
#if not keyword_set(gstruc) then gstruc=tstr else         ;  gstruc = [gstruc,tstr] 
    print('gstruc ',time.time()-t1 
    redo_fail = 0 
 
# count gt 0 
 
 
# SKIP FITTING PART 
else: 
SKIP: 
 
# creating a dummy structure 
tstr = gstruc_schema 
redo_fail = 0 
redo = 0 
back = 0 
 
if keyword_set(trackplot) : 
gtrackplot,lon,lat,lastlon,lastlat,redo,                                     count,lonr=lonr,latr=latr,pstr=pstr,xstr=xstr,ystr=ystr 
 
# skip the fitting part 
 
 
# FINISHING UP THE TRACKING STRUCTURE 
npar = len((tstr.par)(*)) 
 
# increasing btrack parameter arrays 
if (count > 0): 
nbpar = len((*(!btrack.data))[0].par) 
;if (npar > nbpar):gbtrack,btrack,tstr# this is done in btrack_add.pro 
if (npar > nbpar) : 
gbtrack,track,tstr 
 
track.par[0:npar-1] = (tstr.par)(*) 
track.rms = rms 
track.noise = noise 
track.redo_fail = redo_fail 
 
# UPDATING THE TRACKING STRUCTURE 
btrack_add,track 
#if count eq 0 then btrack = track 
#if count gt 0 then btrack = [btrack,track] 
 
count = count + 1 
 
# Saving the last position 
lastlon = lon 
lastlat = lat 
 
print('this iteration ',time.time()-t00 
 
# SAVING THE STRUCTURES, periodically 
if not keyword_set(savestep) : 
savestep=50 
nsave = savestep 
if (long64(count)/long64(nsave) == long64(count)/float(nsave)): 
print('SAVING DATA!' 
MWRFITS,(*(!gstruc.data))[0:!gstruc.count-1],file,/create 
MWRFITS,(*(!btrack.data))[0:!btrack.count-1],file,/silent# append 
gstruc = !gstruc & btrack = !btrack 
SAVE,gstruc,btrack,file=restore_file 
undefine,gstruc,btrack 
 
#ngstruc = n_elements(gstruc) 
#print,gstruc[ngstruc-1].lon,gstruc[ngstruc-1].lat 
 
# Debugging 
#print,strtrim(n_elements(gstruc),2),' Gaussians' 
#dum = where(gstruc.lat le lat,nbelow) 
#print,strtrim(nbelow,2),' Gaussians <= lat=',strtrim(lat,2) 
#print,'Fraction <=  ',strtrim(double(nbelow)/double(n_elements(gstruc)),2) 
#dum = where(gstruc.lon eq 0 and gstruc.lat eq 0,nzero) 
#print,strtrim(nzero,2),' Gaussians at lon=lat=0' 
#print,' ' 
# 
#if lat ge 2 and nzero gt 10 then stop,'too many Gaussians at lon=lat=0' 
# 
#if count gt 20000L and lon eq 0 and lat eq 0 then begin 
#  print,'PROBLEM. lon/lat reset to 0' 
#  retall 
#endif 
 
#if count gt 40 then stop 
#stop 
 
BOMB: 
 
# while endflag eq 0 
 
# FINAL SAVE 
print(str(len(gstruc),2),' final Gaussians' 
print('Saving data to ',file 
MWRFITS,(*(!gstruc.data))[0:!gstruc.count-1],file,/create 
MWRFITS,(*(!btrack.data))[0:!btrack.count-1],file,/silent# append 
gstruc = !gstruc & btrack = !btrack 
SAVE,gstruc,btrack,file=restore_file 
undefine,gstruc,btrack 
 
print('dt = ',str(time.time()-tstart,2),' sec.' 
 
BOMB: 
 
#stop 
 
