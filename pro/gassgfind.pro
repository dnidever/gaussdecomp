function gfind,gstruc,lon,lat,ind=ind,rms=rms,noise=noise,$
         par=par,lonr=lonr,latr=latr


; This function helps find a latitude and longitude in
; the gaussian components structure.
;
;  INPUT
;   gstruc  Structure of gaussian parameters
;   lon     Longitude to search for
;   lat     Latitude to search for
;   lonr    Two element array of longitude limits, lonr=[lonmin,lonmax]
;   latr    Two element array of latitude limits, latr=[latmin,latmax]
;
;  OUTPUT
;   flag    The function value is either 0 or 1:
;             1 - the position exists in the structure
;             0 - the position does NOT exist in the structure
;            -1 - any problems
;   ind     Index of gaussians in gstruc with the desired position
;   rms     RMS of gaussian fit at the desired position
;   noise   Noise level at desired position
;   par     Parameters of gaussians in gstruc with the desired position
;
; When there are any problems in this program it returns:
;  flag = -1
;  rms = -1
;  noise = -1
;  par = -1
;  ind = -1
;
; Created by David Nidever April 2005

; assume bad until proven otherwise
flag = -1
rms = -1
noise = -1
par = -1
ind = -1

ngstruc = n_elements(gstruc)
nlon = n_elements(lon)
nlat = n_elements(lat)

; Bad Input Values
if (n_params() eq 0) or (ngstruc eq 0) or (nlon eq 0) or (nlat eq 0) then begin
  print,'Syntax - f = gfind(gstruc,lon,lat,ind=ind,rms=rms,noise=noise,'
  print,'                   par=par,lonr=lonr,latr=latr)'
  return,-1
endif

; Making sure it's the right structure
if (n_tags(gstruc) eq 0) then return,-1
tags = tag_names(gstruc)
if (n_elements(tags) lt 6) then return,-1
comp = (tags eq ['LON','LAT','RMS','NOISE','PAR','SIGPAR'])
if ((where(comp ne 1))(0) ne -1) then return,-1

; setting the ranges
if keyword_set(lonr) then begin
  lon0 = lonr(0)
  lon1 = lonr(1)
endif else begin
  lon0 = 0.
  lon1 = 359.5
endelse
if keyword_set(latr) then begin
  lat0 = latr(0)
  lat1 = latr(1)
endif else begin
  lat0 = -90.
  lat1 = 90.
endelse

if (lon lt lon0) or (lon gt lon1) or (lat lt lat0) or (lat gt lat1) then begin
  flag = -1
  rms = -1
  noise = -1
  par = -1
  ind = -1
  goto,BOMB
endif

; Looking for the position
ind = where(gstruc.lon eq lon and gstruc.lat eq lat,nind)

; Found something, getting the values
if nind gt 0 then rms = first_el(gstruc(ind).rms)
if nind gt 0 then noise = first_el(gstruc(ind).noise)
if nind gt 0 then par = (gstruc(ind).par)(*)

; Nothing found
if nind eq 0 then rms = -1
if nind eq 0 then noise = -1
if nind eq 0 then par = -1

; Setting flag value
if nind eq 0 then flag = 0
if nind gt 0 then flag=1

BOMB:

;stop

return,flag

end
