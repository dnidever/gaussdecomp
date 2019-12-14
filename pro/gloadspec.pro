pro gloadspec,cubefile,lon,lat,spec,v,glon,glat,npts=npts,noise=noise

; This programs loads the spectrum for GDRIVER.PRO
; LON is the longitude index of the datacube
; LAT is the latitude index of the datacube

undefine,spec,v,glon,glat

if n_elements(lon) eq 0 or n_elements(lat) eq 0 then begin
  print,'Syntax - gloadspec,lon,lat,spec,v0,glon,glat,npts=npts,noise=noise'
  return
endif

; The first time
DEFSYSV,'!gauss',exists=gauss_exists
if not gauss_exists then begin
  print,'LOADING DATACUBE'
  FITS_READ,cubefile,cube,head
  FITS_ARRAYS,head,v0,glon0,glat0
  szcube = size(cube)
  ; flip velocity
  ;cube = reverse(cube,3)
  ;v0 = reverse(v0)
  ; Flip longitude
  ;cube = reverse(cube,2)
  ;glon0 = reverse(glon0)
  glon2d = glon0#(fltarr(szcube[3])+1.0)
  glat2d = (fltarr(szcube[2])+1.0)#glat0

  ; missing data are set to roughly -1.52
  ;bd = where(cube lt -1,nbd)
  ;cube[bd] = !values.f_nan

  DEFSYSV,'!gauss',{cube:cube,head:head,v:v0,glon2d:glon2d,glat2d:glat2d}

  undefine,cube,glon2d,glat2d

Endif

; Load the spectrum
glon = !gauss.glon2d[lon,lat]
glat = !gauss.glat2d[lon,lat]
v = !gauss.v
spec = reform(!gauss.cube[*,lon,lat])  ; v,lon,lat

; Only want good points
gd = where(finite(spec) eq 1,ngd)
mingood = 10            ; minimum number of GOOD spectral points (not NAN)

; No spectrum at this position, skip
if ngd lt mingood then begin
  rms = 999999.
  noise = 999999.
  undefine,spec,v
  npts = 0
  return
end else begin
  spec = spec[gd]
  v = v[gd]
  npts = n_elements(v)
  hinoise,v,spec,noise  ; get the noise
endelse

;stop

end
