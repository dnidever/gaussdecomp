pro gplot_9panel,str,lon,lat,color=color,save=save,file=file,$
          yrange=yrange,xrange=xrange,tit=tit,noresid=noresid,$
          xthick=xthick,ythick=ythick,charsize=charsize,charthick=charthick,$
          thick=thick,xtit=xtit,ytit=ytit,position=position,$
          normal=normal,device=device

; Plots the gaussian fit.  You can plot just the gaussians
; by typing: gplot,v,0,par  (or anything for spec)
;
;  INPUT
;   v        Array of velocities
;   y        Array of spectral points
;   par      Array of gaussian parameters
;   /color   Plot in color
;   /save    Save to postscript file
;   file     Name of postscript file
;   xrange   Range of x-values
;   yrange   Range of y-values
;   tit      Plot title
;   /noresid Don't overplot the residuals
;
;  OUTPUT
;   None
;
; Created by David Nidever April 2005

nstr = n_elements(str)
nlon = n_elements(lon)
nlat = n_elements(lat)

; Bad Input Values
if (n_params() eq 0) or (nstr eq 0) or (nlon eq 0) or (nlat eq 0) then begin
  print,'Syntax - gplot,v,y,par,color=color,save=save,file=file,'
  print,'         yrange=yrange,xrange=xrange,tit=tit'
  return
endif


if not keyword_set(file) then file='gaussfit'
if keyword_set(save) then ps_open,file,color=color

; colors for my screen
red = 250
lred = 210
green = 190000
orange = 310000
yellow = 450000
blue = -25000
lblue = -15000
purple = -20000
white = -1
backgr = 0.
coarr = [green,lblue,yellow,blue,purple,lred,orange]

; setting for postscript
if keyword_set(color) then begin
  loadct,39
  black=0
  purple=30
  blue=60
  aqua=80
  green=155   ;135
  yellow=200
  orange=225
  white=0
  red=250
  lred=240
  ;backg=white
  backg=yellow
  coarr = [green,aqua,yellow,blue,purple,lred,orange]
endif

minlon = min(str.lon)
maxlon = max(str.lon)
minlat = min(str.lat)
maxlat = max(str.lat)

gd = where(str.lon eq lon and str.lat eq lat,ngd)
if ngd eq 0 then begin
  print,'LON=',strtrim(lon,2),' LAT=',strtrim(lat,2),' NOT IN STRUCTURE'
  return
end

!p.multi=[0,3,3]
erase

for i=0,2 do begin
  ilat = lat+i-1
  rdmatthews5,ilat,frame,glon,glat,vel


  for j=0,2 do begin
    ilon = lon+j-1

    count = i*3 + j
    !p.multi=[9-count,3,3]
    if count gt 0 then noerase=1 else noerase=0

    gdi = where(str.lon eq ilon and str.lat eq ilat,ngdi)

    ; We've got some
    if ngdi gt 0 then begin

      y = reform(frame[ilon,*])
      par = str[gdi].par
      gplot,vel,y,par,noerase=noerase,/noresid,yr=[-0.1,1],/noannot
      
      xyouts,min(vel)+5,0.85,'('+strtrim(ilon,2)+','+strtrim(ilat,2)+')',align=0

      ;stop

    endif else begin

    endelse

  end ; for j
end ; for i

if keyword_set(save) then ps_close

if keyword_set(stp) then stop

end
