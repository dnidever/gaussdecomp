pro make_head_cube

; Make a cube for the head of the stream and both MCs

files = 'data/'+['gass_10_0_1111111001_sub3.fits','gass_10_0_1111111001_sub4.fits','gass_10_0_1111111001_sub5.fits',$
                 'gass_-30_0_1111111001_sub1.fits','gass_-30_0_1111111001_sub2.fits','gass_-30_0_1111111001_sub3.fits',$
                 'gass_-30_0_1111111001_sub4.fits','gass_-30_0_1111111001_sub5.fits']
nfiles = n_elements(files)

for i=0,nfiles-1 do begin

  fits_read,files[i],cube1,head
  ;head = headfits(files[i])
  fits_arrays,head,mlon1,mlat1,vel1,/silent
  vel1 /= 1e3
  print,minmax(mlon1)
  
  ; 600-999, 26-356 km/s
  cube1 = cube1[*,*,600:999]
  vel1 = vel1[600:999]

  if i eq 0 then begin
    cube = cube1
    mlon = mlon1
    mlat = mlat1
    vel = vel1
  endif else begin
    g = where(mlon1 lt min(mlon),ngd)
    lo = g[0]
    cube1 = cube1[lo:*,*,*]
    mlon1 = mlon1[lo:*]
    
    ; new array
    sz = size(cube)
    sz1 = size(cube1)
    nx = sz[1]+sz1[1]
    cube0 = cube
    cube = fltarr(nx,sz[2],sz[3])
    cube[0:sz[1]-1,*,*] = cube0
    cube[sz[1]:*,*,*] = cube1
    mlon = [mlon,mlon1]
    undefine,cube0,cube1

  endelse

  ;stop

endfor

;save,cube,mlon,mlat,vel,file='magstream_head_cube,dat'

; This is roughly the same region as Figure 8 of For et al.
loadcol,3
displayc,reverse(transpose(cube[0:700,300:800,186]),2),mlat[300:800],mlon[0:700],/log,min=0.1

; Filament 1
; mlon1b, mlat1
restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil1path_mlonmlat.dat'
; cen1, mlon1
restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil1path_mloncen.dat'
ml1 = scale_vector(findgen(1000),min([mlon1,mlon1b]),max([mlon1,mlon1b]))
interp,reverse(mlon1b),reverse(mlat1),ml1,mb1
interp,reverse(mlon1),reverse(cen1),ml1,v1

; Filament 2
restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil2path_mlonmlat.dat'
restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil2path_mloncen.dat'
ml2 = scale_vector(findgen(1000),min([mlon2,mlon2b]),max([mlon2,mlon2b]))
interp,reverse(mlon2b),reverse(mlat2),ml2,mb2
interp,reverse(mlon2),reverse(cen2),ml2,v2

setdisp
oplot,mlat1,mlon1b,co=250
oplot,mlat2,mlon2b,co=150

stop

end
