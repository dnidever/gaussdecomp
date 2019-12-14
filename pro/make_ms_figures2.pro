pro make_ms_figures2

; Make figure like Figure 8 in For et al.

setdisp

; cube,mlon,mlat,vel
restore,'magstream_head_cube,dat'

; Filament 1
; mlon1b, mlat1
;restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil1path_mlonmlat.dat'
restore,'/net/dl1/users/dnidever/net/halo/dln5q/doradus/papers/hims1/figs/fil1path_mlonmlat.dat'
; cen1, mlon1
;restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil1path_mloncen.dat'
restore,'/net/dl1/users/dnidever/net/halo/dln5q/doradus/papers/hims1/figs/fil1path_mloncen.dat'
ml1 = scale_vector(findgen(1000),min([mlon1,mlon1b]),max([mlon1,mlon1b]))
interp,reverse(mlon1b),reverse(mlat1),ml1,mb1
interp,reverse(mlon1),reverse(cen1),ml1,v1
mag2gal,ml1,mb1,gl1,gb1
;writecol,'ms_filament1.txt',gl1,gb1,v1,ml1,mb1,fmt='(2F10.4,F8.2,2F10.4)'
;writeline,'ms_filament1.txt','#   GLON      GLAT      VLSR     MLON      MLAT',/prepend

; Filament 2
;restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil2path_mlonmlat.dat'
;restore,'/n/Dnidever1/net/halo/dln5q/doradus/papers/hims1/figs/fil2path_mloncen.dat'
restore,'/net/dl1/users/dnidever/net/halo/dln5q/doradus/papers/hims1/figs/fil2path_mlonmlat.dat'
restore,'/net/dl1/users/dnidever/net/halo/dln5q/doradus/papers/hims1/figs/fil2path_mloncen.dat'
ml2 = scale_vector(findgen(1000),min([mlon2,mlon2b]),max([mlon2,mlon2b]))
interp,reverse(mlon2b),reverse(mlat2),ml2,mb2
interp,reverse(mlon2),reverse(cen2),ml2,v2
mag2gal,ml2,mb2,gl2,gb2
;writecol,'ms_filament2.txt',gl2,gb2,v2,ml2,mb2,fmt='(2F10.4,F8.2,2F10.4)'
;writeline,'ms_filament2.txt','#   GLON      GLAT      VLSR     MLON      MLAT',/prepend

velarr = [180.0, 184.0, 188.0, 193.0, 197.0, 201.0, 205.0, 209.0, 213.0, 217.0, 221.0, 225.0, $
          230.0, 234.0, 238.0, 242.0, 246.0, 250.0, 254.0, 259.0, 263.0, 267.0, 271.0, 275.0]
nvelarr = n_elements(velarr)

!p.font = 0
ps_open,'ms_figures1a',/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=12
;!p.multi=[0,3,4]

for i=0,nvelarr-1 do begin

  dum = closest(velarr[i],vel,ind=vind)
  vind = vind[0]

  ; This is roughly the same region as Figure 8 of For et al.
  loadcol,3
  black = fsc_color('black',0)
  ;if i eq 0 or i eq 12 then noerase=0 else noerase=1
  noerase = 1
  row = i / 3
  col = i mod 3
  if row gt 3 then begin
    if row eq 4 and col eq 0 then noerase=0
    row-=4  ; second page
  endif
  row = 3-row  ; flip
  x0 = 0.07
  y0 = 0.05
  x1 = 0.99
  y1 = 0.99
  dx = (x1-x0)/3.0
  dy = (y1-y0)/4.0
  pos = [dx*col+x0, dy*row+y0, dx*(col+1)+x0, dy*(row+1)+y0]
  if row eq 0 then xtit = 'B!dMS!n' else xtit=''
  if row eq 0 then undefine,xtickformat else xtickformat='(A1)'
  if col eq 0 then ytit = 'L!dMS!n' else ytit=''
  if col eq 0 then undefine,ytickformat else ytickformat='(A1)'
  display,reverse(transpose(cube[0:700,300:800,vind]),2),mlat[300:800],mlon[0:700],/log,min=0.1,xtit=xtit,ytit=ytit,$
          xtickformat=xtickformat,ytickformat=ytickformat,noerase=noerase,position=pos
  ;         tit='V!dLSR!n='+stringize(vel[vind],ndec=1)+' km/s',noerase=noerase,position=pos

  loadct,39,/silent

  ; Filament 1
  oplot,mb1,ml1,co=250,linestyle=1,thick=1
  ind1 = where(abs(v1-vel[vind]) lt 10,nind1)
  for j=0,nind1-1 do oplot,[mb1[ind1[j]]],[ml1[ind1[j]]],ps=8,sym=0.2,co=250

  oplot,mb2,ml2,co=150,linestyle=1,thick=1
  ind2 = where(abs(v2-vel[vind]) lt 10,nind2)
  for j=0,nind2-1 do oplot,[mb2[ind2[j]]],[ml2[ind2[j]]],ps=8,sym=0.2,co=150

  xyouts,-24,8,stringize(velarr[i],ndec=1)+' km/s',align=0,charsize=1.1

  if i eq 11 then begin
    ps_close
    !p.font = 0
    ps_open,'ms_figures1b',/color,thick=4,/encap
    device,/inches,xsize=9.5,ysize=12
    ;!p.multi=[0,3,4]
  endif

  ;stop

endfor

ps_close
!p.multi = 0

stop

end
