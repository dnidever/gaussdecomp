pro fix_gass25b

; construct GSTRUC from BTRACK

restore,'gass2-5b.dat'
undefine,gstruc

dum = {lon:0.0,lat:0.0,rms:0.0,noise:0.0,par:fltarr(3),sigpar:fltarr(3),glon:0.0,glat:0.0}
gstruc = replicate(dum,500000L)

g = where(btrack.par[0] lt 9e5,ng)
count = 0
for i=0LL,ng-1 do begin
  if i mod 1000 eq 0 then print,i
  btrack1 = btrack[g[i]]
  par = btrack1.par
  ind = where(par lt 9e5,nind)
  ngauss = nind/3
  new = replicate(dum,ngauss)
  new.lon = btrack1.lon
  new.lat = btrack1.lat
  new.rms = btrack1.rms
  new.noise = btrack1.noise
  ;new.glon = btrack1.glon
  ;new.glat = btrack1.glat
  for j=0,ngauss-1 do new[j].par=par[j*3:(j+1)*3-1]
  ; NO sigpar/glon/glat in btrack
  gstruc[count:count+ngauss-1] = new
  ;push,gstruc,new
  count += ngauss
  ;stop
endfor

gstruc = gstruc[0:count-1]

stop

end
